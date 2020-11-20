import requests
import time
from tenacity import retry, stop_after_attempt, wait_exponential

# GraphQL query code
query_code = '''
query Gene($geneId: String, $geneSymbol: String, $referenceGenome: ReferenceGenomeId!, #$variantId: String!, $dataset: DatasetId = gnomad_r2_1 
){
  gene(gene_id: $geneId, gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
    reference_genome
    gene_id
    gene_version
    symbol
    name
    canonical_transcript_id
    hgnc_id
    omim_id
    canonical_transcript_id
    clinvar_variants {
      variant_id
      reference_genome
      chrom
      pos
      ref
      alt
      clinical_significance
      clinvar_variation_id
      gold_stars
      major_consequence
    }
    gnomad_constraint {
      oe_lof
      oe_lof_lower
      oe_lof_upper
      oe_mis
      oe_mis_lower
      oe_mis_upper
      lof_z
      mis_z
      pLI
    }
  }
}
# query Variant($variantId: String!
# ){
#   variant(dataset: gnomad_r2_1, variantId: $variantId) {
#     colocatedVariants
#     exome {
#       ac_hemi
#       ac_hom
#       ac
#     }
#   }
# }
'''

result_dict = {}
REQS_PER_SEC = 15

class GnomADQuery:
    def __init__(self, identifier, category="gene", assembly="GRCh37"):
        self.category = category
        if self.category == "gene":
            self.query = '''
                query Gene($geneId: String, $geneSymbol: String, $referenceGenome: ReferenceGenomeId!, #$variantId: String!, $dataset: DatasetId = gnomad_r2_1 
                ){
                  gene(gene_id: $geneId, gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
                    reference_genome
                    gene_id
                    gene_version
                    symbol
                    name
                    canonical_transcript_id
                    hgnc_id
                    omim_id
                    canonical_transcript_id
                    clinvar_variants {
                      variant_id
                      reference_genome
                      chrom
                      pos
                      ref
                      alt
                      clinical_significance
                      clinvar_variation_id
                      gold_stars
                      major_consequence
                    }
                    gnomad_constraint {
                      oe_lof
                      oe_lof_lower
                      oe_lof_upper
                      oe_mis
                      oe_mis_lower
                      oe_mis_upper
                      lof_z
                      mis_z
                      pLI
                    }
                  }
                }'''
            self.query_variables = '{"geneId": "' + identifier + '", "referenceGenome": "' + assembly + '"}'
        else:
            self.variant = identifier.replace(":", "-")
            self.query = """query Variant($variantId: String!
                            ){
                              variant(dataset: gnomad_r2_1, variantId: $variantId) {
                                colocatedVariants
                                exome {
                                  ac_hemi
                                  ac_hom
                                  ac
                                }
                                genome {
                                  ac_hemi
                                  ac_hom
                                  ac
                                }
                              }
                            }"""
            self.query_variables = '{"variantId": "' + self.variant + '"}'
        self.req_count = 0
        self.last_req = 0

    @retry(reraise=True,
           stop=stop_after_attempt(10),
           wait=wait_exponential(multiplier=1, min=1, max=5))
    def gnomad_sparql_request(self):
        """General function for handling API communication.
        """
        r = requests.post(url="https://gnomad.broadinstitute.org/api?",
                          json={'query': self.query, 'variables': self.query_variables},
                          headers={'Accept': 'application/vnd.cap-collectif.preview+json'})

        # if status code == 200
        if r.ok:
            return r
        # if error "too many requests" retry after given time
        else:
            # Reraise if some sort of error occurs.
            print(f"Some error occured while requesting VEP data. Retrying...\n"
                  f"{r.status_code}: {r.reason}")
            raise IOError("There has been an issue with a variant while requesting gnomAD.")

    # requests and returns gnomAD gene information
    def get_gnomad_info(self):
        time.sleep(0.2)
        # check if we need to rate limit ourselves
        if self.req_count >= REQS_PER_SEC:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        self.req_count += 1

        try:
            r = self.gnomad_sparql_request()
            status_code = r.status_code
        except IOError:
            print("Some problem with gnomAD API request!")
            status_code = 400

        if status_code == 200:
            status_code = 200
            decoded = r.json()
            result_dict = self.recursion(decoded)
            if self.category == "gene":
                if result_dict is not None:
                    if result_dict.get("oe_mis") and result_dict.get("oe_mis_lower") and result_dict.get("oe_mis_upper"):
                        result_dict["oe_mis_interval"] = "{value}  [{lower} - {upper}]".format(
                            value=round(result_dict.get("oe_mis"), 2),
                            lower=round(result_dict.get("oe_mis_lower"), 2),
                            upper=round(result_dict.get("oe_mis_upper"), 2))
                    if result_dict.get("oe_lof") and result_dict.get("oe_lof_lower") and result_dict.get("oe_lof_upper"):
                        result_dict["oe_lof_interval"] = "{value}  [{lower} - {upper}]".format(
                            value=round(result_dict.get("oe_lof"), 2),
                            lower=round(result_dict.get(
                                "oe_lof_lower"), 2),
                            upper=round(result_dict.get(
                                "oe_lof_upper"), 2))
                else:
                    status_code = 495
            else:
                if "ac_hemi_exome" in result_dict.keys():
                    result_dict["ac_hemi"] = result_dict.get("ac_hemi_exome") + result_dict.get("ac_hemi_genome")
                    result_dict["ac_hom"] = result_dict.get("ac_hom_exome") + result_dict.get("ac_hom_genome")
                    result_dict["ac"] = result_dict.get("ac_exome") + result_dict.get("ac_genome")

                if self.variant[0] in ["X", "x"]:
                    result_dict["male_count"] = result_dict.get("ac_hemi") or 0
                    total_allele_count = result_dict.get("ac") or 0
                    result_dict["female_count"] = total_allele_count - result_dict.get("male_count")

                if "errors" in result_dict.keys():
                    for i in range(len(result_dict.get("errors"))):
                        try:
                            if result_dict.get("errors")[i].get("message") == "Variant not found":
                                not_in_gnomad_dict = {}
                                not_in_gnomad_dict["ac_hemi"] = 0
                                not_in_gnomad_dict["ac_hom"] = 0
                                not_in_gnomad_dict["ac"] = 0
                                return not_in_gnomad_dict, status_code
                        except (AttributeError, TypeError, IndexError):
                            print("Could not retrieve Error Code from gnomAD Server response!")
            return result_dict, status_code

        else:
            return None, status_code


    # formats data into a one dimensional dictionary
    def recursion(self, dict, suffix="", sub_dirs=["exome", "genome"]):
        for key, value in dict.items():
            if type(value) == type(dict):
                if key in sub_dirs:
                    suffix = "_" + key
                self.recursion(value, suffix)
            else:
                if isinstance(value, float):
                    result_dict[key + suffix] = round(value, 2)
                else:
                    result_dict[key + suffix] = value
        return result_dict

# print(GnomADQuery("X:153032470:G:A", "variant").get_gnomad_info())
# print(GnomADQuery("1-55516888-G-GA", "variant").get_gnomad_info())
# gnomad_variant_result = gnomad_variant_instance.get_gnomad_info()
# gnomad_instance = GnomADQuery("ENSG00000198753")
# gnomad_info = gnomad_instance.get_gnomad_info()
# print(gnomad_info)
# print(gnomad_variant_result)