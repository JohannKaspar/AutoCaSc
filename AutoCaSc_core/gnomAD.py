import requests
import time

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
                  }
                }"""
            self.query_variables = '{"variantId": "' + self.variant + '"}'
        self.req_count = 0
        self.last_req = 0

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

        r = requests.post(url="https://gnomad.broadinstitute.org/api?",
                                 json={'query': self.query, 'variables': self.query_variables},
                                 headers={'Accept': 'application/vnd.cap-collectif.preview+json'})

        if r.status_code == 200:
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
                if result_dict.get("ac_hemi") is not None and result_dict.get("ac") is not None:
                    result_dict["male_count"] = result_dict.get("ac_hemi")
                    result_dict["female_count"] = result_dict.get("ac") - result_dict.get("ac_hemi")
                else:
                    result_dict["male_count"] = 0
                    result_dict["female_count"] = 0

            return result_dict, status_code


        elif r.status_code == 429:
            if 'Retry-After' in r.headers:
                retry = r.headers['Retry-After']
                time.sleep(float(retry))
                self.get_gnomad_info()
        else:
            return None, r.status_code


    # formats data into a one dimensional dictionary
    def recursion(self, dict):
        for key, value in dict.items():
            if type(value) == type(dict):
                self.recursion(value)
            else:
                if isinstance(value, float):
                    result_dict[key] = round(value, 2)
                else:
                    result_dict[key] = value
        return result_dict

# print(GnomADQuery("X:153032470:G:A", "variant").get_gnomad_info())
# gnomad_variant_result = gnomad_variant_instance.get_gnomad_info()
# gnomad_instance = GnomADQuery("ENSG00000198753")
# gnomad_info = gnomad_instance.get_gnomad_info()
# print(gnomad_info)
# print(gnomad_variant_result)