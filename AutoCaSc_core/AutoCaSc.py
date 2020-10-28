import math
import re
import time
from pathlib import Path
from statistics import mean

import click
import pandas as pd
import requests
from numpy import isnan

from AutoCaSc_core.gnomAD import GnomADQuery
from AutoCaSc_core.tools import safe_get, filterTheDict

AUTOCASC_VERSION = 0.95
ROOT_DIR = str(Path(__file__).parent) + "/data/"

gene_scores = pd.read_csv(ROOT_DIR + "all_gene_data.csv")
# This loads the calculated gene scores and gnomad constraints.

class AutoCaSc:
    """AutoCaSc_core is a tool for quantifying plausibility of gene variants as a cause for neurodevelopmental delay.
    This class is the core program. A variant is annotated by initializing the class.
    Class attributes correspond to variant parameters ranging from position to candidate score.
    """

    def __init__(self, variant, inheritance="other",
                 family_history=False, other_impact="unknown", assembly="GRCh37", CADD_based=False):
        """This function assigns basic parameters and initializes the scoring process.

        :param variant: the variant including position and alternative sequence
        :param inheritance: segregation/zygosity --> de_novo, homo, comphet, other, ad_inherited, x_linked
        :param family_history: are there multiple affected family members? if x_linked affected male family members?
        :param other_impact: impact of the corresponding variant if it is compound heterozygous
        :param assembly: either GRCh37 or GRCh38
        """

        self.mutationtaster_converted_rankscore = None
        self.version = AUTOCASC_VERSION
        self.variant = variant
        self.inheritance = inheritance
        self.comphet_id = None
        self.family_history = family_history
        self.other_impact = other_impact
        self.assembly = assembly
        self.status_code = 200  # initial value is set to 200 = all good

        if self.assembly == "GRCh37":
            self.server = "http://grch37.rest.ensembl.org"  # API endpoint for GRCh37
        else:
            self.server = "http://rest.ensembl.org"  # API endpoint for GRCh38

        # ToDo das hier hübscher machen?

        # assign initial "None" to all parameters
        self.ada_score = None
        self.rf_score = None
        self.maxentscan_ref = None
        self.maxentscan_alt = None
        self.maxentscan_decrease = None
        self.cadd_phred = 0
        self.vcf_string = None
        self.ref_seq_vep = None
        self.gnomad_frequency = None
        self.maf = None
        self.impact = None
        self.gerp_rs_rankscore = None
        self.transcript = None
        self.consequence = None
        self.hgvsc_change = None
        self.protein = None
        self.hgvsp_change = None
        self.polyphen_prediction = None
        self.ada_consequence = None
        self.maxentscan_consequence = None
        self.hgvsc_change = None
        self.pLI = None
        self.mis_z = None
        self.oe_lof = None
        self.oe_lof_lower = None
        self.oe_lof_upper = None
        self.oe_mis = None
        self.oe_mis_lower = None
        self.oe_mis_upper = None
        self.oe_mis_interval = None
        self.oe_lof_interval = None
        self.explanation_dict = {}
        self.inheritance_score = 0
        self.variant_score = 0
        self.literature_score = 0
        self.gene_attribute_score = 0
        self.candidate_score = 0
        self.sift_converted_rankscore = None
        self.mutationtaster_converted_rankscore = None
        self.mutationassessor_rankscore = None
        self.mgi_score = None

        self.check_variant_format()  # this function is called to check if the entered variant is valid

        if self.status_code != 401:  # if variant format is valid (not 401) continue
            if self.variant_format == "vcf":
                # definition of instance variables
                self.chromosome = self.variant.split(":")[0]
                self.pos_start = self.variant.split(":")[1]
                self.ref_seq = self.variant.split(":")[2]
                self.pos_end = str(int(self.pos_start) + len(self.ref_seq) - 1)
                # calculates end of sequence using length of reference sequence to start position
                self.alt_seq = self.variant.split(":")[3]
                self.ext_url = "/vep/human/region/" + ":".join([str(self.chromosome), str(self.pos_start), str(
                    self.pos_end)]) + "/" + self.alt_seq + "?&hgvs=1&vcf_string=1&numbers=1&dbscSNV=1&CADD=1&GeneSplicer=1&SpliceRegion=1&MaxEntScan=1&canonical=1&dbNSFP=SIFT_converted_rankscore,GERP%2B%2B_RS_rankscore,MutationAssessor_rankscore,MutationTaster_converted_rankscore&content-type=application/json"

            elif self.variant_format == "hgvs":
                # Assigning sequences in case they have been entered. This enables the program to check if
                # the entered reference sequence matches the VEP reference sequence at the given position.
                if ">" in self.variant:
                    self.ref_seq = re.search(r"(?<=\d)[CTGA]+(?=>)", variant).group()
                    self.alt_seq = re.search(r"(?<=>)[CTGA]+", variant).group()
                else:
                    self.ref_seq = None
                    self.alt_seq = None
                self.ext_url = "/vep/human/hgvs/" + self.variant + "?&hgvs=1&vcf_string=1&numbers=1&dbscSNV=1&CADD=1&GeneSplicer=1&SpliceRegion=1&MaxEntScan=1&canonical=1&dbNSFP=SIFT_converted_rankscore,GERP%2B%2B_RS_rankscore,MutationAssessor_rankscore,MutationTaster_converted_rankscore&content-type=application/json"

            self.annotate()  # this method call initiates the annotation of the given variant

        try:
            # gene contraint gnomad
            if self.status_code != 498:
                # 498 would be if no matching transcript index has been found, usually occuring when entering an
                # intergenic variant.
                if self.gene_id in gene_scores.ensemble_id.to_list():
                    # checks if there are contraint values in gnomad for the given gene
                    for parameter in ["pLI", "mis_z", "oe_lof", "oe_lof_lower", "oe_lof_upper", "oe_mis",
                                      "oe_mis_lower",
                                      "oe_mis_upper"]:
                        self.__dict__[parameter] = round(
                            gene_scores.loc[gene_scores.ensemble_id == self.gene_id, parameter].values[0], 2)
                        if isnan(self.__dict__[parameter]):
                            self.__dict__[parameter] = None

                    self.oe_mis_interval = f"{self.oe_mis}  [{self.oe_mis_lower} - {self.oe_mis_upper}]"
                    self.oe_lof_interval = f"{self.oe_lof}  [{self.oe_lof_lower} - {self.oe_lof_upper}]"

            self.get_scores(CADD_based)  # this calls the function for scoring the variant

        except AttributeError:
            pass
            # print(f"ATTRIBUTE ERROR IN AUTOCASC: {self.status_code}")
        except TypeError:
            pass
            # print(f"TYPE ERROR IN AUTOCASC: {self.status_code}")

    # master calculation function, calls sub calculations as described in Büttner et al and returns dict of category results
    def get_scores(self, CADD_based=False):
        """This method calls all the scoring functions and assigns their results to class attributes.
        """
        self.explanation_dict = {}

        if CADD_based:
            if self.gene_id in gene_scores.ensemble_id.to_list():
                # If the gene_id is in the computed gene score table, its results are assigned to the class attributes.
                self.literature_score = round(
                    gene_scores.loc[gene_scores.ensemble_id == self.gene_id, "weighted_score"].values[0], 2)
            else:
                # If not, the values are set to 0.
                for score in ["literature_score", "pubtator_score", "gtex_score", "denovo_rank_score", "disgenet_score",
                              "mgi_score", "string_score"]:
                    self.__dict__[score] = 0.

            if self.cadd_phred == 0.0 or self.cadd_phred == 0:
                self.candidate_score = 0
            else:
                self.candidate_score = (self.literature_score * math.log10(self.cadd_phred))

        else:
            self.rate_inheritance()
            self.rate_pli_z()

            self.rate_impact()
            self.rate_in_silico()
            self.rate_conservation()
            self.rate_frequency()
            self.variant_score = round(sum(filter(None, [self.impact_score,
                                                         self.in_silico_score,
                                                         self.conservation_score,
                                                         self.frequency_score
                                                         ])), 2)

            if self.gene_id in gene_scores.ensemble_id.to_list():
                # If the gene_id is in the computed gene score table, its results are assigned to the class attributes.
                self.literature_score = round(
                    6.0 * gene_scores.loc[gene_scores.ensemble_id == self.gene_id, "weighted_score"].values[0], 2)
                # Assigning the plausibility subscores for being able to call them individually.
                for score in ["pubtator_score", "gtex_score", "denovo_rank_score", "disgenet_score",
                              "mgi_score", "string_score"]:
                    self.__dict__[score] = round(gene_scores.loc[gene_scores.ensemble_id == self.gene_id, score].values[0],
                                                 2)
            else:
                # If not, the values are set to 0.
                for score in ["literature_score", "pubtator_score", "gtex_score", "denovo_rank_score", "disgenet_score",
                              "mgi_score", "string_score"]:
                    self.__dict__[score] = 0.

            self.calculate_canidate_score()

    def calculate_canidate_score(self):
        candidate_score_list = [i if i is not None else 0 for i in [
            self.inheritance_score,
            self.gene_attribute_score,
            self.variant_score,
            self.literature_score]]
        self.candidate_score = round(sum(candidate_score_list), 2)

    def rate_inheritance(self):
        """This function scores zygosity/segregation.
        """
        if self.inheritance == "de_novo":
            # ToDo weniger auf de_novo geben wenn family_history positiv?
            self.inheritance_score, self.explanation_dict["inheritance"] = 2, "de novo    -->    2"
        if self.inheritance == "other":
            self.inheritance_score, self.explanation_dict["inheritance"] = 0, "other    -->    0"
        if self.inheritance == "ad_inherited":
            self.inheritance_score, self.explanation_dict["inheritance"] = 0, "inherited autosomal dominant    -->    0"
        if self.inheritance == "comphet":
            if self.family_history in [False, "no"]:
                self.inheritance_score, self.explanation_dict["inheritance"] = 1,\
                    "compound heterozygous, one affected child    -->    1"
            elif self.family_history in [True, "yes"]:
                self.inheritance_score, self.explanation_dict["inheritance"] = 3,\
                    "compound heterozygous, multiple affected children    -->    3"
        if self.inheritance == "homo":
            if self.family_history is True:
                self.inheritance_score, self.explanation_dict["inheritance"] = 3,\
                    "homo, multiple affected children    -->    3"
            else:
                self.inheritance_score, self.explanation_dict["inheritance"] = 2, "homo, one affected child    -->    2"

        if self.inheritance == "x_linked":
            if self.family_history == True:
                # ToDo was wenn weibliche Verwandte betroffen, extra differenzieren?
                self.inheritance_score, self.explanation_dict["inheritance"] = 2,\
                    "x_linked and another maternal male relative    -->    2"
            else:
                self.inheritance_score, self.explanation_dict["inheritance"] = 1, "x_linked and a boy    -->    1"

    def rate_pli_z(self):
        """This function scores gnomAD constraints (pLI/Z).
        """
        if self.inheritance in ["homo", "comphet"]:
            self.gene_attribute_score, self.explanation_dict["pli_z"] = 0, "recessive    -->    0"
        else:
            # if impact is not None:
            if self.impact == "high":
                if self.pLI is not None:
                    if self.pLI < 0.5:
                        if self.inheritance == "de_novo":
                            self.gene_attribute_score, self.explanation_dict[
                                "pli_z"] = -2, "de novo LoF & pLI < 0.5    -->    -2"
                        else:
                            self.gene_attribute_score, self.explanation_dict["pli_z"] = 0, "LoF & pLI < 0.5    -->    0"
                    elif self.pLI < 0.9:
                        if self.inheritance == "de_novo":
                            self.gene_attribute_score, self.explanation_dict[
                                "pli_z"] = 0, "de novo LoF & 0.5 <= pLI < 0.9    -->    0"
                        else:
                            self.gene_attribute_score, self.explanation_dict[
                                "pli_z"] = 0.5, "LoF & pLI <= 0.5 pLI < 0.9    -->    0.5"
                    else:
                        self.gene_attribute_score, self.explanation_dict["pli_z"] = 1, "LoF & pLI >= 0.9    -->    1"
                else:
                    self.gene_attribute_score, self.explanation_dict["pli_z"] = 0, "no data on pLI    --> 0"
            elif self.impact == "moderate":
                if self.mis_z is not None:
                    if self.mis_z < 0:
                        self.gene_attribute_score, self.explanation_dict["pli_z"] = 0, "missense & Z < 0    -->    0"
                    elif self.mis_z >= 0 and self.mis_z < 2.5:
                        self.gene_attribute_score, self.explanation_dict[
                            "pli_z"] = 0.5, "missense & 0 <= Z < 2.5    -->    0.5"
                    elif self.mis_z >= 2.5:
                        self.gene_attribute_score, self.explanation_dict["pli_z"] = 1, "missense & Z >= 2.5    -->    1"
                else:
                    self.gene_attribute_score, self.explanation_dict["pli_z"] = 0, "no data on Z score    -->    0"
            else:
                self.gene_attribute_score, self.explanation_dict["pli_z"] = 0, "low impact or no data on impact --> 0"

    def rate_impact(self):
        """This function scores the variants predicted impact.
        """
        self.impact_score, self.explanation_dict["impact"] = 0, f"impact {self.impact}    -->    0"

        if self.impact == "moderate":
            self.impact_score, self.explanation_dict["impact"] = 0, "impact moderate    -->    0"
        if self.impact == "high":
            self.impact_score, self.explanation_dict["impact"] = 2, "impact high, heterozygous    -->    2"
            if self.inheritance == "homo":
                self.impact_score, self.explanation_dict["impact"] = 3, "impact high, biallelic    -->    3"
            if self.other_impact == "high":
                self.impact_score, self.explanation_dict["impact"] = 3, "compound heterozygous, both impacts high    -->    3"

    def rate_in_silico(self):
        """This function scores in silico predictions.
        """
        # ToDo ranked scores benutzt? checken!
        score_list = []
        affected_splice_counter = []

        self.in_silico_score, self.explanation_dict["in_silico"] = 0, f"{self.impact} and no prediction values --> 0"

        if self.sift_converted_rankscore is not None:
            score_list.append(float(self.sift_converted_rankscore))
        if self.mutationtaster_converted_rankscore is not None:
            score_list.append(float(self.mutationtaster_converted_rankscore))
        if self.mutationassessor_rankscore is not None:
            score_list.append(float(self.mutationassessor_rankscore))

        # "average of in silico"
        if score_list:
            self.in_silico_score, self.explanation_dict["in_silico"] = round(mean(score_list), 2),\
                    f"mean in silico score    -->    {round(mean(score_list), 2)}"
        elif self.impact == "high":
            # If there are no in silico scores for the variant and its impact is high, we assume a score of 1 point.
            self.in_silico_score, self.explanation_dict["in_silico"] = 1, "LoF --> 1"

        # In case of absent in silico values and a moderate or low impact: consider splice site predictions, if existing
        else:
            if self.rf_score is not None:
                if self.rf_score >= 0.6:  # cutoff from "Leipzig guidelines for use of Alamut splice site prediction tools"
                    affected_splice_counter.append(1)
            if self.ada_score is not None:
                if self.ada_score >= 0.6:
                    affected_splice_counter.append(1)
            if self.maxentscan_decrease is not None:
                if self.maxentscan_decrease <= -0.15:
                    affected_splice_counter.append(1)

            if not affected_splice_counter and self.impact == "moderate":
                self.in_silico_score = 0.5
                self.explanation_dict["in_silico"] = "missense, no available in silico values    -->    0.5"
            if sum(affected_splice_counter) == 1:  # "splicing affected in one program"
                self.in_silico_score = 0.5
                self.explanation_dict["in_silico"] = "splicing affected in one program    -->    0.5"
            if sum(affected_splice_counter) >= 2:  # "splicing affected in two or more programs"
                self.in_silico_score = 1
                self.explanation_dict["in_silico"] = "splicing affected in two or more programs    -->    1"

    def rate_conservation(self):
        """This function scores the variants conservation.
        """
        if self.gerp_rs_rankscore is not None:
            self.conservation_score, self.explanation_dict["conservation"] = self.gerp_rs_rankscore,\
                   f"GERP++ rankscore {self.gerp_rs_rankscore}    -->    {self.gerp_rs_rankscore}"

        # "estimation in case of absent GERP++"
        else:
            if self.impact == "high":
                self.conservation_score, self.explanation_dict["conservation"] = 1, "LoF    -->    1"
            else:
                self.conservation_score, self.explanation_dict["conservation"] = 0, "no data!    -->    0"
        # TODO "use UCSC and give Numbers between 0 - 1"?

    def rate_frequency(self):
        """This function scores the variants frequency.
        """
        self.frequency_score, self.explanation_dict["frequency"] = 0, "other    -->    0"

        if self.inheritance in ["de_novo", "ad_inherited"]:
            if self.maf < 0.000005:
                self.frequency_score, self.explanation_dict["frequency"] = 1,\
                    f"{self.inheritance} & MAF < 0.000005    -->    1"
            elif self.maf < 0.00002:
                self.frequency_score, self.explanation_dict["frequency"] = 0.5,\
                    f"{self.inheritance} & MAF < 0.00002    -->    0.5"
            else:
                self.frequency_score, self.explanation_dict["frequency"] = 0, "de novo & MAF > 0.00002    -->    0"

        if self.inheritance in ["homo", "comphet"]:
            if self.maf < 0.00005:
                self.frequency_score, self.explanation_dict["frequency"] = 1,\
                    "autosomal recessive & MAF < 0.00005    -->    1"
            elif self.maf < 0.0005:
                self.frequency_score, self.explanation_dict["frequency"] = 0.5,\
                    "autosomal recessive & MAF < 0.0005    -->    0.5"

        if self.inheritance == "x_linked":
            gnomad_variant_result, _ = GnomADQuery(self.vcf_string, "variant").get_gnomad_info()
            self.male_count = gnomad_variant_result.get("male_count") or 0
            self.female_count = gnomad_variant_result.get("female_count") or 0
            if self.male_count <= 1 and self.female_count >= 5:
                self.frequency_score, self.explanation_dict["frequency"] = 2,\
                    "X linked and discrepancy of MAF in gnomAD between males and females (max.1/min.5)    -->    2"
                #TODO DAS HIER DRUNTER WURDE ADAPTIERT, MUSS MIT RAMI BESPROCHEN WERDEN
            elif self.male_count == 0:
                self.frequency_score = 1
                self.explanation_dict["frequency"] = "X linked no hemizygote in gnomAD  -->    2"
            else:
                self.frequency_score = 0
                self.explanation_dict["frequency"] = "X linked and at least one hemizygote in gnomad"

    def check_variant_format(self):
        """This function checks if the entered variant matches either HGVS or VCF format after doing some formatting.
        """
        input = self.variant
        input = re.sub(r"^[\W]", "", input)  # deletes characters that are not a word character
        input = re.sub(r"Chr|chr", "", input)

        input_vcf = re.sub(r"\s+", " ", input)  # substitutes whitespaces and tabs with one whitespace
        input_vcf = re.sub(r"[^A-Za-z0-9](?!$)", ":", input_vcf)  # substitutes special characters by a colon
        input_vcf = re.sub(r":{2,}", ":", input_vcf)

        input_hgvs = re.sub(r"\s+", "", input)

        if re.fullmatch(r"(\d{1,2}|X):\d+:[a-zA-Z]+:([a-zA-Z]*|-)+", input_vcf) is not None:
            self.variant_format = "vcf"
            self.variant = input_vcf
        elif re.fullmatch(
                r"^.+:[cgp]\.\d+([_+-]\d+)?(_\d+)?([+-]\d+)?([CTGA]+>[CTGA]+|del[CTGA]*|dup[CTGA]*)?(ins[CTGA]+)?",
                input_hgvs) is not None:
            self.variant_format = "hgvs"
            self.variant = input_hgvs
        else:
            self.status_code = 401

    # core function for annotation
    def annotate(self):
        """This function requests the VEP API in order to annotate a given variant and selects relevant parameters.
        :return:
        """
        response_decoded = self.retrieve_data()  # calls function to request VEP data on the given variant
        if self.status_code == 200:
            transcript_index = self.get_transcript_index(response_decoded)
            # get the index of the transcript to consider for further annotations
            if self.status_code == 200:
                self.assign_results(response_decoded, transcript_index)

    # function that makes the server request
    # @retry(tries=10, delay=5, backoff=1.2, jitter=(1, 3), logger=logger)
    def retrieve_data(self, retries=0):
        """General function for handling API communication.
        It takes care of certain server responses, for example if too many requests have been sent.
        Return obtained data as a dict.
        """
        REQS_PER_SEC = 15  # limit requests to 15 per seconds
        req_count = 0
        last_req = 0

        # check if we need to rate limit ourselves
        if req_count >= REQS_PER_SEC:
            delta = time.time() - last_req
            if delta < 1:
                time.sleep(1 - delta)
            last_req = time.time()
            req_count = 0
        req_count += 1

        # try to pull the requested data and check if request was successful
        r = requests.get(self.server + self.ext_url, headers={"content-type": "application/json"})

        # if content is not empty return content
        if r.status_code == 200:
            return safe_get(r.json(), 0)

        # if error "too many requests" retry after given time
        elif r.status_code == 429:
            if 'Retry-After' in r.headers:
                retry = r.headers['Retry-After']
                time.sleep(float(retry))
                self.retrieve_data()
        elif r.status_code == 503 and retries < 5:
            time.sleep(5)
            self.retrieve_data(retries + 1)
        else:
            # if some sort of error occurs
            if "matches reference" in str(r.content, "utf-8"):
                # Return None if the given alternative sequence matches the GRCh reference sequence.
                self.status_code = 496
                return None
            else:
                # Return None if some sort of different error occurs.
                self.status_code = r.status_code
                return None

    # stores and returns relevant parameters in a dictionary
    def assign_results(self, variant_anno_data, transcript_index):
        """This function filters and formats the received annotation results.

        :param variant_anno_data: received annotation data
        :param transcript_index: index of the transcript to work with from transcript_consequences
        """
        selected_transcript_consequences = variant_anno_data["transcript_consequences"][transcript_index]

        if variant_anno_data.get("vcf_string") is not None:
            self.vcf_string = re.sub(r"-", ":", variant_anno_data.get("vcf_string"))
        if self.variant_format == "vcf":
            # formats given variant into VCF format
            self.id = re.sub(r"[^a-zA-Z^0-9-]", ":", variant_anno_data.get("id"))
            # extracts reference sequence from variant id in vcf format
            self.ref_seq_vep = safe_get(self.id.split(":"), 2)
        else:
            self.id = variant_anno_data.get("id")
        if self.ref_seq is None:
            self.ref_seq = variant_anno_data.get("vcf_string").split("-")[2]
            self.alt_seq = variant_anno_data.get("vcf_string").split("-")[3]

        # gets varint gnomad exome frequency and MAF
        self.gnomad_frequency, self.maf = self.get_frequencies(variant_anno_data.get("colocated_variants"))

        # This list of parameters is assigned to corresponding class attributes.
        key_list = ["gene_symbol",
                    "amino_acids",
                    "gene_id",
                    "cadd_phred",
                    "sift_converted_rankscore",
                    "mutationtaster_converted_rankscore",
                    "mutationassessor_rankscore",
                    "ada_score",
                    "rf_score",
                    "maxentscan_ref",
                    "maxentscan_alt"]

        for key in key_list:
            if selected_transcript_consequences.get(key) is not None:
                value = selected_transcript_consequences.get(key)
                if isinstance(value, float):
                    self.__dict__[key] = round(selected_transcript_consequences.get(key), 2)
                else:
                    self.__dict__[key] = selected_transcript_consequences.get(key)

        if selected_transcript_consequences.get("impact") is not None:
            self.impact = selected_transcript_consequences.get("impact").lower()

        if selected_transcript_consequences.get("gerp++_rs_rankscore") is not None:
            self.gerp_rs_rankscore = round(selected_transcript_consequences.get("gerp++_rs_rankscore"), 2)

        if selected_transcript_consequences.get("consequence_terms") is not None:
            self.consequence = safe_get(selected_transcript_consequences.get("consequence_terms"), 0)
        if selected_transcript_consequences.get("hgvsc") is not None:
            self.transcript = safe_get(selected_transcript_consequences.get("hgvsc").split(":"), 0) or \
                              selected_transcript_consequences.get("transcript_id")
            self.hgvsc_change = safe_get(selected_transcript_consequences.get("hgvsc").split(":"), 1)
        if selected_transcript_consequences.get("consequence_terms") is not None:
            self.consequence = re.sub(r"_", " ", safe_get(selected_transcript_consequences.get("consequence_terms"), 0))

        if selected_transcript_consequences.get("hgvsp") is not None:
            self.protein = safe_get(selected_transcript_consequences.get("hgvsp").split(":"), 0)
            self.hgvsp_change = safe_get(selected_transcript_consequences.get("hgvsp").split(":"), 1)

        if selected_transcript_consequences.get("polyphen_prediction") is not None:
            self.polyphen_prediction = re.sub(r"_", " ", selected_transcript_consequences.get("polyphen_prediction"))

        # cutoffs correspond to Leipzig guidelines (Alamut)
        if self.ada_score is not None:
            if self.ada_score >= 0.6:
                self.ada_consequence = "splicing affected"
            else:
                self.ada_consequence = "splicing not affected"

        if self.rf_score is not None:
            if self.rf_score >= 0.6:
                self.rf_consequence = "splicing affected"
            else:
                self.rf_consequence = "splicing not affected"

        if self.maxentscan_ref is not None and self.maxentscan_alt is not None:
            self.maxentscan_decrease = (self.maxentscan_alt - self.maxentscan_ref) / self.maxentscan_ref
            if self.maxentscan_decrease <= -0.15:
                self.maxentscan_consequence = "splicing affected"
            else:
                self.maxentscan_consequence = "splicing not affected"

        # checks if reference sequences match and returns result dictionary and status code accordingly
        if self.ref_seq is not None:
            if self.ref_seq_vep is not None:
                if self.ref_seq_vep != self.ref_seq:
                    self.status_code = 201

    # gets variant frequency in gnomAD exomes and minor allele frequency
    def get_frequencies(self, colocated_variants):
        """This function checks the annotation data for the variants frequency in gnomad and its
        minor allele frequency. If it doesn't find a MAF it sets it to the variants frequency.


        :param colocated_variants: slice of variant annotation results
        :return: variant frequency in gnomad, minor allele frequency
        """
        gnomad_frequency = 0
        maf = 0

        if colocated_variants is None:
            return gnomad_frequency, maf
        else:
            # getting variant frequency
            for colocated_list in colocated_variants:
                frequencies_dict = colocated_list.get("frequencies")
                if frequencies_dict is not None:
                    if frequencies_dict.get(self.alt_seq) is not None:
                        gnomad_frequency = frequencies_dict.get(self.alt_seq).get("gnomad") or 0
                        # TODO 'NM_002642.3:c.422C>T' freq for "aa" but not "gnomad" --> calculate total frequency?

            # getting MAF
            for colocated_dict in colocated_variants:
                if colocated_dict.get("minor_allele_freq") is not None:
                    maf = colocated_dict.get("minor_allele_freq")
            if maf == 0:
                maf = gnomad_frequency

            return gnomad_frequency, maf

    # selects relevant transcript from given transcripts
    def get_transcript_index(self, response_dec):
        """This function selects the transcript index to consider for further calculations, based on
        1) the transcript with the most severe consequence
        2) the canonical transcript (if existing if existing and equally severe consequence
        :param response_dec: VEP response dictionary with a list of transcripts
        :return: transcript index to consider for further calculations
        """

        impact_severity_dict = {"MODIFIER": 0, "LOW": 1, "MODERATE": 2, "HIGH": 3}

        transcript_df = pd.DataFrame()
        # selection prioritisation: 1) most severe consequence    2) canonical
        try:
            for i, transcript in enumerate(response_dec["transcript_consequences"]):
                transcript_df.loc[i, "transcript_id"] = i
                transcript_df.loc[i, "impact_level"] = impact_severity_dict.get(transcript.get("impact"))
                transcript_df.loc[i, "biotype"] = transcript.get("biotype")
                transcript_df.loc[i, "canonical"] = transcript.get("canonical") or 0
            transcript_df = transcript_df.sort_values(by=["impact_level", "canonical"],
                                                      ascending=[False, False])
            transcript_df = transcript_df.reset_index(drop=True)
            if len(transcript_df) > 1 and transcript_df.loc[0, "impact_level"] == transcript_df.loc[1, "impact_level"] \
                    and transcript_df.loc[0, "canonical"] == transcript_df.loc[1, "canonical"]:
                if "protein_coding" in transcript_df.biotype.unique():
                    transcript_df = transcript_df.loc[transcript_df.biotype == "protein_coding"]
                    if len(transcript_df) > 1:
                        print("CAVE! Two protein_coding transcripts are affected!")
                transcript_df.reset_index(inplace=True, drop=True)
            return int(transcript_df.loc[0, "transcript_id"])

        except AttributeError:
            self.status_code = 499
            return None

        except KeyError:
            self.status_code = 498  # intergenic variant
            return None

        except TypeError:
            self.status_code = 497
            return None


# print(AutoCaSc_core("2:46803383:G:A", "de_novo").candidate_score)  # benign
# print(AutoCaSc_core("1:1737948:A:G", "de_novo").candidate_score)  # pathogenic
# print(AutoCaSc_core("X:153040798:C:A", "de_novo").candidate_score)  # benign
# print(AutoCaSc_core("1:3732936:A:G", "de_novo").candidate_score)  # pathogenic
# print(AutoCaSc_core("1:7725246:G:A", "de_novo").candidate_score)  # pathogenic
# print(AutoCaSc_core("1:2234850:G:A", "de_novo").candidate_score)  # pathogenic
# print(AutoCaSc_core("1:3755660:A:C", "de_novo").candidate_score)  # pathogenic
# print(AutoCaSc_core("1:1959699:G:A", "other").candidate_score)  # benign, high expressesed, GABA receptor.......
# print(AutoCaSc_core("CBLB:c.1822C>T", "other").candidate_score)
# print(AutoCaSc_core("NM_001113498.2:c.794T>A", "de_novo").candidate_score)
# print(AutoCaSc_core("chr4:175812275:T:C", inheritance="de_novo", assembly="GRCh38").candidate_score)
# print(AutoCaSc_core("3:149619809:xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx:-", "de_novo").candidate_score)  # benign, high expressesed, GABA receptor.......


@click.group(invoke_without_command=True)  # Allow users to call our app without a command
@click.pass_context
@click.option('--verbose', '-v',
              is_flag=True,
              help="Increase output verbosity level")
@click.option("--output_path", "-o",
              help="Output path for csv-file.")
def main(ctx, verbose, output_path):
    group_commands = ['single', 'batch']
    """
            AutoCaSc_core is a command line tool that helps scoring
            variants of NDD patients for their plausible pathogenicity.\n
            example: python AutoCaSc_core.py batch -i input/file/path
    """

    if ctx.invoked_subcommand is None:
        # No command supplied
        # Inform user on the available commands when running the app
        click.echo("Specify one of the commands below")
        print(*group_commands, sep='\n')

    ctx.obj['VERBOSE'] = verbose
    ctx.obj['OUTPUT_PATH'] = output_path

@click.pass_context
def score_variants(ctx, variants, inheritances, fam_histories):
    # verbose = ctx.obj.get('VERBOSE')
    verbose = True
    results_df = pd.DataFrame()
    variant_dict = {}
    for _variant, _inheritance, _fam_history in zip(variants, inheritances, fam_histories):
        variant_dict[_variant] = AutoCaSc(_variant, _inheritance, _fam_history, assembly="GRCh37")

    if "comphet" in inheritances:
        gene_dict = {}
        for _variant, _inheritance in zip(variants, inheritances):
            if _inheritance == "comphet":
                gene_id = variant_dict.get(_variant).gene_id
                if gene_id:
                    gene_dict[_variant] = gene_id
        # comphet_df = pd.DataFrame.from_dict(gene_dict, orient='index').reset_index()
        # comphet_df.columns = ['variant', 'gene_id']

        for _variant in variants:
            try:
                variant_gene = gene_dict.get(_variant)
                other_variant = list(filterTheDict(gene_dict, variant_gene, _variant).keys())[0]
                other_impact = variant_dict.get(other_variant).impact
                variant_dict[_variant].other_impact = other_impact
                variant_dict[_variant].get_scores()
            except IndexError:
                pass

        comphet_id = 0  # unique identifier for clear identification of corresponding compound heterozygous variants
        for _variant in variants:
            if variant_dict.get(_variant).comphet_id:
                continue
            variant_gene = gene_dict.get(_variant)
            other_variant = list(filterTheDict(gene_dict, variant_gene, _variant).keys())[0]

            combined_variant_score = mean([variant_dict.get(_variant).variant_score,
                                           variant_dict.get(other_variant).variant_score])
            variant_dict.get(_variant).variant_score = combined_variant_score
            variant_dict.get(other_variant).variant_score = combined_variant_score
            variant_dict.get(_variant).comphet_id = comphet_id
            variant_dict.get(other_variant).comphet_id = comphet_id
            variant_dict.get(_variant).calculate_canidate_score()
            variant_dict.get(other_variant).calculate_canidate_score()
            comphet_id += 1


    for _variant, _inheritance, _fam_history in zip(variants, inheritances, fam_histories):
        variant_instance = variant_dict.get(_variant)
        results_df.loc[_variant, "candidate_score"] = variant_instance.candidate_score
        if verbose:
            results_df.loc[_variant, "inheritance_score"] = variant_instance.inheritance_score
            results_df.loc[_variant, "gene_score"] = variant_instance.gene_attribute_score
            results_df.loc[_variant, "variant_score"] = variant_instance.variant_score
            results_df.loc[_variant, "literature_score"] = variant_instance.literature_score
        results_df.loc[_variant, "inheritance_mode"] = _inheritance
        results_df.loc[_variant, "comphet_id"] = variant_instance.comphet_id
        results_df.loc[_variant, "status_code"] = str(int(variant_instance.status_code))

    results_df.index.names = ["variant"]
    # if ctx.obj.get('OUTPUT_PATH'):
    #     results_df.to_csv(ctx.obj.get('OUTPUT_PATH'))
    return results_df

@main.command("batch")
@click.option("--input_file", "-i",
              type=click.Path(exists=True),
              help="Path to file containing variants.")
def score_batch(input_file):
    lines = open(input_file).readlines()
    lines = [line_.strip("\n ") for line_ in lines]
    try:
        variants = [line_.split(",")[0].strip() for line_ in lines]
        inheritances = [line_.split(",")[1].strip() for line_ in lines]
        fam_histories = [line_.split(",")[2].strip() for line_ in lines]
    except IndexError:
        #ToDo error handling
        variants = [line_.strip("\n ") for line_ in lines]
        inheritances = ["de_novo" for _ in lines]

    results_df = score_variants(variants, inheritances, fam_histories)

    click.echo(results_df.head(len(results_df)))


@main.command("single")
@click.option("--variant", "-vcf",
              required=True,
              help="Variant in either VCF or HGVS format.")
@click.option("--inheritance", "-ih",
              default="de_novo",
              help="Inheritance mode of variant.")
@click.option("--family_history", "-f",
              default="no",
              help="Are there any relatives that are affected too?")
def score_single(variant, inheritance, family_history):
    results_df = score_variants([variant], [inheritance], [family_history])
    click.echo(results_df.head(len(results_df)))

# @main.command("single")
# @click.option("--variant", help="Single variant to score.")
# @click.option("--inheritance", default="de_novo", help="Enter inheritance mode. Default is de novo.")


if __name__ == "__main__":
    score_single(["--variant", "1:7725246:G:A",
                 "-ih", "de_novo",
                 "-f", "yes"])
    # score_batch(["--input_file", "/home/johann/variant_test_file.txt"])
    main(obj={})
