import os
import requests
import pickle
from tenacity import retry, stop_after_attempt, wait_exponential, wait_random
from tools import write_new_api_request

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

# result_dict = {}
REQS_PER_SEC = 15

class GnomADQuery:
    def __init__(self, identifier, category="gene", assembly="GRCh37", path_to_request_cache_dir=None):
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
        self.assembly = assembly
        self.variant = identifier
        self.result_dict = {}
        self.path_to_request_cache_dir = path_to_request_cache_dir


    @retry(reraise=True,
           stop=stop_after_attempt(10),
           wait=wait_exponential(multiplier=2, min=1, max=30))
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
            print(f"GNOMAD ERROR '{r.status_code}: {r.reason}' for {self.query_variables}. Retrying...")
            raise IOError("There has been an issue with a variant while requesting gnomAD.")

    @retry(stop=stop_after_attempt(5),
           wait=wait_random(0.1, 1))
    def open_pickle_file(self):
        self.gnomad_requests = {}
        with open(f"{self.path_to_request_cache_dir}gnomad_requests_{self.assembly}", "rb") as gnomad_requests_file:
            self.gnomad_requests = pickle.load(gnomad_requests_file)


    # requests and returns gnomAD gene information
    def get_gnomad_info(self):
        r = None
        if self.path_to_request_cache_dir:
            self.gnomad_requests = {}
            # trying to load stored gnomad requests
            try:
                self.open_pickle_file()
                if self.gnomad_requests.get(self.variant):
                    r = self.gnomad_requests.get(self.variant)
            except pickle.UnpicklingError:
                print("could not open gnomad pickle")

            try:
                if r.status_code == 200:
                    pass
                else:
                    r = None
            except AttributeError:
                r = None

        if r is None:
            try:
                r = self.gnomad_sparql_request()
            except IOError:
                print("Some problem with gnomAD API request!")
                status_code = 400
            if r.status_code == 200 and self.path_to_request_cache_dir:
                if self.gnomad_requests != {}:
                    new_gnomad_request = {self.variant: r}
                    write_new_api_request(f"{self.path_to_request_cache_dir}tmp/gnomad", new_gnomad_request)
        if r is not None:
            status_code = r.status_code

        if status_code == 200:
            decoded = r.json()
            self.recursion(decoded)
            if self.category == "gene":
                if self.result_dict is not None:
                    if self.result_dict.get("oe_mis") and self.result_dict.get("oe_mis_lower") and self.result_dict.get("oe_mis_upper"):
                        self.result_dict["oe_mis_interval"] = "{value}  [{lower} - {upper}]".format(
                            value=round(self.result_dict.get("oe_mis"), 2),
                            lower=round(self.result_dict.get("oe_mis_lower"), 2),
                            upper=round(self.result_dict.get("oe_mis_upper"), 2))
                    if self.result_dict.get("oe_lof") and self.result_dict.get("oe_lof_lower") and self.result_dict.get("oe_lof_upper"):
                        self.result_dict["oe_lof_interval"] = "{value}  [{lower} - {upper}]".format(
                            value=round(self.result_dict.get("oe_lof"), 2),
                            lower=round(self.result_dict.get(
                                "oe_lof_lower"), 2),
                            upper=round(self.result_dict.get(
                                "oe_lof_upper"), 2))
                else:
                    status_code = 495
            else:
                self.result_dict["ac_hemi"] = (self.result_dict.get("ac_hemi_exome") or 0) \
                                         + (self.result_dict.get("ac_hemi_genome") or 0)
                self.result_dict["ac_hom"] = (self.result_dict.get("ac_hom_exome") or 0) \
                                        + (self.result_dict.get("ac_hom_genome") or 0)
                self.result_dict["ac"] = (self.result_dict.get("ac_exome") or 0) \
                                    + (self.result_dict.get("ac_genome") or 0)
                if self.variant[0] in ["X", "x"]:
                    self.result_dict["male_count"] = self.result_dict.get("ac_hemi") or 0
                    total_allele_count = self.result_dict.get("ac") or 0
                    self.result_dict["female_count"] = total_allele_count - self.result_dict.get("male_count")
            return self.result_dict, status_code
        else:
            return None, status_code


    # formats data into a one dimensional dictionary
    def recursion(self, dict, suffix=""):
        sub_dirs = ["exome", "genome"]
        for key, value in dict.items():
            if type(value) == type(dict):
                if key in sub_dirs:
                    suffix = "_" + key
                self.recursion(value, suffix)
            else:
                if isinstance(value, float):
                    self.result_dict[key + suffix] = round(value, 2)
                else:
                    self.result_dict[key + suffix] = value

#x = GnomADQuery("1:55516888:G:GA", "variant").get_gnomad_info()
#print("")
# print(GnomADQuery("1-55516888-G-GA", "variant").get_gnomad_info())
# gnomad_variant_result = gnomad_variant_instance.get_gnomad_info()
# gnomad_instance = GnomADQuery("ENSG00000198753")
# gnomad_info = gnomad_instance.get_gnomad_info()
# print(gnomad_info)
# print(gnomad_variant_result)