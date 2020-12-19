import math
import re
import time
from pathlib import Path
from statistics import mean
from tenacity import retry, stop_after_attempt, wait_exponential
import click
import pandas as pd
import requests
from numpy import isnan, poly1d, product
from AutoCaSc_core.gnomAD import GnomADQuery
from AutoCaSc_core.tools import safe_get, filterTheDict

# from gnomAD import GnomADQuery
# from tools import safe_get, filterTheDict

AUTOCASC_VERSION = 0.96
ROOT_DIR = str(Path(__file__).parent) + "/data/"

gene_scores = pd.read_csv(ROOT_DIR + "all_gene_data.csv")
# This loads the calculated gene scores and gnomad constraints.

class AutoCaSc:
    """AutoCaSc_core is a tool for quantifying plausibility of gene variants as a cause for neurodevelopmental delay.
    This class is the core program. A variant is annotated by initializing the class.
    Class attributes correspond to variant parameters ranging from position to candidate score.
    """
    def __init__(self, variant,
                 inheritance="other",
                 parent_affected=False,
                 has_sibling=False,
                 cosegregating=False,
                 sex="XY",
                 other_autocasc_obj=None,
                 other_impact="unknown",
                 other_variant=None,
                 assembly="GRCh37",
                 transcript_num=None):
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
        self.parent_affected = parent_affected
        self.has_sibling = has_sibling
        self.cosegregating = cosegregating  # this is meant to be True if a sibling is affected and has the same variant
        self.sex = sex
        self.other_autocasc_obj = other_autocasc_obj
        self.other_impact = other_impact
        self.other_variant = other_variant
        self.assembly = assembly
        self.status_code = 200  # initial value is set to 200 = all good

        self.multiple_transcripts = False
        self.transcript_num = transcript_num

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
        if not self.variant_format == "incorrect":  # if variant format is valid (not 401) continue
            self.retrieve_data()
        else:
            self.status_code = 401

    def retrieve_data(self):
        self.get_vep_data()  # this method call initiates the annotation of the given variant
        if self.status_code == 200:
            try:
                # 498 = no matching transcript index has been found (e.g. variant is intergenic)
                self.get_gnomad_constraint()
                self.get_gnomad_counts()  # gets allele counts in gnomad
                if self.other_variant:
                    self.other_autocasc_obj = AutoCaSc(variant=self.other_variant,
                                                       inheritance=self.inheritance,
                                                       parent_affected=self.parent_affected,
                                                       has_sibling=self.has_sibling,
                                                       cosegregating=self.cosegregating,
                                                       sex=self.sex,
                                                       assembly=self.assembly,
                                                       transcript_num=self.transcript_num
                                                       )

            except (AttributeError, TypeError) as e:
                raise e
                # print(f"variant: {self.variant}\n"
                #       f"status_code: {self.status_code}\n"
                #       f"{e}\n")

    def create_url(self):
        if self.variant_format == "vcf":
            # definition of instance variables
            self.chromosome, self.pos_start, self.ref_seq, self.alt_seq = self.variant.split(":")
            self.pos_end = str(int(self.pos_start) + len(self.ref_seq) - 1)
            # calculates end of sequence using length of reference sequence to start position
            self.ext_url = "/vep/human/region/" + ":".join([str(self.chromosome), str(self.pos_start), str(
                self.pos_end)]) + "/" + self.alt_seq + "?&hgvs=1&vcf_string=1&numbers=1&dbscSNV=1&CADD=1&GeneSplicer=1&SpliceRegion=1&MaxEntScan=1&canonical=1&dbNSFP=SIFT_converted_rankscore,GERP%2B%2B_RS_rankscore,MutationAssessor_rankscore,MutationTaster_converted_rankscore&content-type=application/json"

        elif self.variant_format == "hgvs":
            # Assigning sequences in case they have been entered. This enables the program to check if
            # the entered reference sequence matches the VEP reference sequence at the given position.
            if ">" in self.variant:
                self.ref_seq = re.search(r"(?<=\d)[CTGA]+(?=>)", self.variant).group()
                self.alt_seq = re.search(r"(?<=>)[CTGA]+", self.variant).group()
            else:
                self.ref_seq = None
                self.alt_seq = None
            self.ext_url = "/vep/human/hgvs/" + self.variant + "?&hgvs=1&vcf_string=1&numbers=1&dbscSNV=1&CADD=1&GeneSplicer=1&SpliceRegion=1&MaxEntScan=1&canonical=1&dbNSFP=SIFT_converted_rankscore,GERP%2B%2B_RS_rankscore,MutationAssessor_rankscore,MutationTaster_converted_rankscore&content-type=application/json"

    def get_gnomad_constraint(self):
        if self.gene_id in gene_scores.ensemble_id.to_list():
            # checks if there are contraint values in gnomad for the given gene
            for parameter in ["pLI",
                              "mis_z",
                              "oe_lof",
                              "oe_lof_lower",
                              "oe_lof_upper",
                              "oe_mis",
                              "oe_mis_lower",
                              "oe_mis_upper",
                              "obs_hom_lof",
                              "virlof_ar_enrichment"]:
                self.__dict__[parameter] = round(
                    gene_scores.loc[gene_scores.ensemble_id == self.gene_id, parameter].values[0], 2)
                if isnan(self.__dict__[parameter]):
                    self.__dict__[parameter] = None

            self.oe_mis_interval = f"{self.oe_mis}  [{self.oe_mis_lower} - {self.oe_mis_upper}]"
            self.oe_lof_interval = f"{self.oe_lof}  [{self.oe_lof_lower} - {self.oe_lof_upper}]"

    def get_gnomad_counts(self):
        gnomad_variant_result, self.status_code = GnomADQuery(self.vcf_string, "variant").get_gnomad_info()
        if self.status_code == 200:
            self.ac_hom = gnomad_variant_result.get("ac_hom")
            self.allele_count = gnomad_variant_result.get("ac")
            self.n_hemi = gnomad_variant_result.get("ac_hemi")

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
            self.variant_format = "incorrect"
            self.status_code = 401

    @retry(reraise=True,
           stop=stop_after_attempt(15),
           wait=wait_exponential(multiplier=1.3, min=0.1, max=3))
    def vep_api_request(self):
        """General function for handling API communication.
        Return obtained data as a dict.
        """
        url = self.server + self.ext_url

        # try to pull the requested data and check if request was successful
        r = requests.get(url, headers={"content-type": "application/json"})

        # if status code == 200
        if r.ok:
            return safe_get(r.json(), 0), 200
        # if error "too many requests" retry after given time
        else:
            # if some sort of error occurs
            if "matches reference" in str(r.content, "utf-8"):
                # Return None if the given alternative sequence matches the GRCh reference sequence.
                return None, 496
            else:
                # Return None if some sort of different error occurs.
                print(f"VEP ERROR '{r.status_code}: {r.reason}' occured for {self.variant}. Retrying...")
                raise IOError("There has been an issue with a variant.")

    # core function for annotation
    def get_vep_data(self):
        """This function requests the VEP API in order to annotate a given variant and selects relevant parameters.
        :return:
        """
        self.create_url()
        try:
            response_decoded, status_code = self.vep_api_request()
        except IOError:
            response_decoded = None
            self.status_code = 400

        if self.status_code == 200:
            transcript_index = self.get_transcript_index(response_decoded)
            # get the index of the transcript to consider for further annotations
            if self.status_code == 200:
                self.assign_results(response_decoded, transcript_index)

    # function that makes the server request
    # @retry(tries=10, delay=5, backoff=1.2, jitter=(1, 3), logger=logger)

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
                transcript_df = transcript_df.loc[transcript_df.impact_level == transcript_df.loc[0, "impact_level"]]
                transcript_df = transcript_df.loc[transcript_df.canonical == transcript_df.loc[0, "canonical"]]
                if "protein_coding" in transcript_df.biotype.unique():
                    transcript_df = transcript_df.loc[transcript_df.biotype == "protein_coding"]
                    if len(transcript_df) > 1:
                        if not self.transcript_num:
                            print("CAVE! Two protein_coding transcripts are affected!")
                            self.multiple_transcripts = len(transcript_df)
                        else:
                            return int(transcript_df.loc[self.transcript_num, "transcript_id"])
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


    # master calculation function, calls subcalculations
    def calculate_candidate_score(self):
        """This method calls all the scoring functions and assigns their results to class attributes.
        """
        if self.status_code == 200:
            self.explanation_dict = {}
            self.factors = []

            if self.inheritance in ["de_novo", "ad_inherited", "unknown"]:  # autosomal dominant branch
                self.score_dominant()

            elif self.inheritance in ["comphet", "homo"]:
                self.score_recessive()

            elif self.inheritance == "x_linked":
                if self.sex == "XX":
                    self.score_recessive()
                elif self.sex == "XY":
                    self.score_x_hemi()

            if self.gene_id in gene_scores.ensemble_id.to_list():
                # If the gene_id is in the computed gene score table, its results are assigned to the class attributes.
                self.literature_score = round(gene_scores.loc[gene_scores.ensemble_id == self.gene_id,
                                                              "weighted_score"].values[0], 2)
                self.factors.append((self.literature_score, "precalculated plausibility score"))
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
                mean_weighted = round(gene_scores["weighted_score"].mean(), 2)
                self.factors.append((mean_weighted, "mean of plausibility scores, not available for this gene"))

            factor_list, explanation_list = [[factor for factor, _ in self.factors],
                                             [explanation for _, explanation in self.factors]]
            self.candidate_score = round(product(factor_list) * 10, 2)

        # for debugging purpose
        # for factor, explanation in self.factors:
        #     print(f"{factor}  because {explanation}")

    def score_dominant(self):
        if self.inheritance == "de_novo":
            if self.allele_count == 0:
                self.factors.append((1, "denovo and not in gnomad"))
                # ToDo add explanations for factors
            elif self.allele_count == 1:
                self.factors.append((0.9, "denovo and once in gnomad"))
            else:
                self.factors.append((0, "denovo and multiple times in gnomad"))
        elif self.inheritance == "unknown":
            if self.allele_count == 0:
                self.factors.append((0.9, "inheritance unknown and not in gnomad"))
            else:
                self.factors.append((0, "inheritance unknown, present in gnoamd"))
        elif self.inheritance == "ad_inherited":
            # hier wird davon ausgegangen, dass ad_inherited nur möglich wenn self.parent_affected == True
            inheritance_factor = 1 - ((self.allele_count + 1) * 0.1)
            self.factors.append((inheritance_factor,
                                 f"inherited autosomal dominant, {self.allele_count}x in gnomad"))

        if self.impact == "high":
            if self.oe_lof_upper:
                self.factors.append((self.get_loeuf_factor(), f"impact high and LOEUF {self.oe_lof_upper}"))
            else:
                self.factors.append((1.0, "LOEUF not available"))
        elif self.impact == "moderate":
            self.factors.append((self.get_Z_factor(), f"impact moderate and Z {self.mis_z}"))
        else:
            self.factors.append((0, "low impact"))  # if impact not moderate or high --> variant == 0

        if not self.cadd_phred == 0:
            self.factors.append((self.get_cadd_factor(), f"CADD {self.cadd_phred}"))
        else:
            self.factors.append((1, f"no data on CADD available"))

    def score_recessive(self):
        if self.impact not in ["high", "moderate"]:
            self.factors.append((0, "low impact"))

        if self.has_sibling:
            if self.cosegregating:
                self.factors.append((1, "variant cosegregating"))
            else:
                self.factors.append((0, "sibling affected but variant not present or "
                                        "not affected but variant present"))
        else:
            self.factors.append((0.95, "no sibling available"))

        if self.inheritance == "homo":
            if self.ac_hom == 0:
                self.factors.append((1, "moderate impact, variant not homozygous in gnomad"))
            elif self.ac_hom == 1:
                self.factors.append((0.5, "moderate impact, variant only once homozygous in gnomad"))
            else:
                self.factors.append((0, f"moderate imapct, variant {self.ac_hom}x in gnomad"))

        if self.inheritance == "comphet":
            if self.other_autocasc_obj:
                if self.other_autocasc_obj.status_code == 200:
                    for autocasc_object in [self, self.other_autocasc_obj]:
                        # check prevalence of homozygous variant in gnomad
                        if self.ac_hom == 0:
                            self.factors.append((1, f"moderate impact, variant {autocasc_object.ac_hom} not homozygous in gnomad"))
                        elif self.ac_hom == 1:
                            self.factors.append((0.5, f"moderate impact, variant {autocasc_object.ac_hom} only once homozygous in gnomad"))
                        else:
                            self.factors.append((0, f"moderate imapct, variant {autocasc_object.ac_hom}x in gnomad"))
                else:
                    self.status_code = self.other_autocasc_obj.status_code
            else:
                self.status_code = 402

        self.factors.append((self.get_gevir_score(), f"GEVIR virlof ar enrichment {self.virlof_ar_enrichment}"))

        if (self.inheritance == "homo") or ((self.inheritance == "x_linked" and self.sex == "XX")):
            if not self.cadd_phred == 0:
                self.factors.append((self.get_cadd_factor(), f"CADD {self.cadd_phred}"))
            else:
                self.factors.append((1, f"no data on CADD available"))
        elif self.inheritance == "comphet":
            if (self.cadd_phred != 0) and (self.other_autocasc_obj.cadd_phred != 0):
                mean_cadd = mean([self.cadd_phred, self.other_autocasc_obj.cadd_phred])
                self.factors.append((self.get_cadd_factor(mean_cadd), f"comphet, one CADD {self.cadd_phred} "
                                                                      f"other CADD {self.other_autocasc_obj.cadd_phred}"))
            else:
                self.factors.append((1, f"no data on CADD available"))

    def score_x_hemi(self):
        if self.n_hemi == 0:
            self.factors.append((self.get_cadd_factor(),
                                 f"x-linked hemizygous, CADD {self.cadd_phred}"))
        else:
            self.factors.append((0,
                                f"x-linked hemizygous, variant {self.n_hemi}x hemizygous in gnomad"))

    def get_loeuf_factor(self):
        if self.oe_lof_upper >= 0.6:
            loeuf_score = 0.
        elif self.oe_lof_upper <= 0.2:
            loeuf_score = 1.
        else:
            loeuf_model = poly1d([218.18181818, -341.01010101,175.33333333,-37.33665224,3.82452381])
            loeuf_score = loeuf_model(self.oe_lof_upper)
        if loeuf_score >= 1.:
            loeuf_score = 1.
        if loeuf_score < 0:
            loeuf_score = 0
        self.loeuf_score = loeuf_score
        return loeuf_score

    def get_Z_factor(self):
        if self.mis_z <= 0.:
            z_score = 0.
        elif self.mis_z >= 3.09:
            z_score = 1.
        else:
            z_model = poly1d([0.00542974, -0.09183755, 0.55584498, 0.00952595])
            z_score = z_model(self.mis_z)
        if z_score >= 1.:
            z_score = 1.
        self.z_score = z_score
        return z_score

    def get_cadd_factor(self, cadd_phred=None):
        if not cadd_phred:
            cadd_phred = self.cadd_phred
        if cadd_phred >= 30:
            cadd_score = 1.
        else:
            cadd_model = poly1d([3.83333333e-05, -3.85000000e-03, 1.28666667e-01, -4.40000000e-01])
            cadd_score = cadd_model(cadd_phred)

        if cadd_score < 0.:
            cadd_score = 0.
        elif cadd_score > 1.:
            cadd_score = 1.
        self.cadd_score = cadd_score
        return cadd_score

    def get_gevir_score(self):
        return self.virlof_ar_enrichment / 2.0


@click.group(invoke_without_command=True)  # Allow users to call our app without a command
@click.pass_context
@click.option('--verbose', '-v',
              is_flag=True,
              help="Increase output verbosity level")
@click.option("--assembly", "-a",
              default="GRCh37",
              help="Reference assembly to use.")
@click.option("--parent_affected", "-pa",
              default="n",
              help="Is at least one of the parents affected either y or n.")
@click.option("--cosegregating", "-c",
              default="no_sibling",
              help="This option can either be no_sibling, cosegregating or not_cosegregating.")
@click.option("--sex", "-s",
              default="XY",
              help="Sex of index patient. Either XX or XY.")
@click.option("--output_path", "-o",
              help="Output path for csv-file.")
def main(ctx, verbose, assembly, parent_affected, cosegregating, sex, output_path):
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

    ctx.obj["VERBOSE"] = verbose
    ctx.obj["ASSEMBLY"] = assembly
    ctx.obj["PARENT_AFFECTED"] = parent_affected
    ctx.obj["SEX"] = sex
    ctx.obj["OUTPUT_PATH"] = output_path

    if cosegregating == "no_sibling":
        ctx.obj["HAS_SIBLING"] = False
        ctx.obj["COSEGREGATING"] = False
    elif cosegregating == "cosegregating":
        ctx.obj["HAS_SIBLING"] = True
        ctx.obj["COSEGREGATING"] = True
    elif cosegregating == "not_cosegregating":
        ctx.obj["HAS_SIBLING"] = True
        ctx.obj["COSEGREGATING"] = False

@click.pass_context
def score_variants(ctx, variants, inheritances):
    if ctx.obj:
        verbose = ctx.obj.get('VERBOSE')
        assembly = ctx.obj.get('ASSEMBLY')
        parent_affected = ctx.obj.get('PARENT_AFFECTED')
        sex = ctx.obj.get('SEX')
        output_path = ctx.obj.get('OUTPUT_PATH')
        has_sibling = ctx.obj.get('HAS_SIBLING')
        cosegregating = ctx.obj.get('COSEGREGATING')
    else:
        verbose = True
        assembly = "GRCh37"
        parent_affected = False
        sex = "XY"
        output_path = None
        has_sibling = False
        cosegregating = False

    results_df = pd.DataFrame()
    variant_dict = {}
    for _variant, _inheritance in zip(variants, inheritances):
        variant_dict[_variant] = AutoCaSc(variant=_variant,
                              inheritance=_inheritance,
                              parent_affected=parent_affected,
                              has_sibling=has_sibling,
                              cosegregating=cosegregating,
                              sex=sex,
                              assembly=assembly)

    #ToDo: check this again, maybe write one function for both vcf and batch CLI

    comphet_genes_dict = {}
    for _variant, _inheritance in zip(variants, inheritances):
        if _inheritance == "comphet":
            gene_id = variant_dict.get(_variant).gene_id
            if gene_id:
                comphet_genes_dict[_variant] = gene_id
        else:
            variant_dict.get(_variant).calculate_candidate_score()
    # comphet_df = pd.DataFrame.from_dict(gene_dict, orient='index').reset_index()
    # comphet_df.columns = ['variant', 'gene_id']

    for _variant in comphet_genes_dict.keys():
        try:
            variant_gene = comphet_genes_dict.get(_variant)
            other_variant = list(filterTheDict(comphet_genes_dict, variant_gene, _variant).keys())[0]
            variant_dict.get(_variant).other_autocasc_obj = variant_dict.get(other_variant)
            variant_dict.get(_variant).calculate_candidate_score()
        except IndexError:
            pass

    # create the final result_df containing scoring results
    for _variant, _inheritance in zip(variants, inheritances):
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
        print(variant_instance.factors)

    results_df.index.names = ["variant"]
    if output_path:
        results_df.to_csv(output_path, index=False)
    return results_df


@main.command("batch")
@click.option("--input_file", "-i",
              type=click.Path(exists=True),
              help="Path to file containing variants.")
def batch(input_file):
    lines = open(input_file).readlines()
    lines = [line_.strip("\n ") for line_ in lines]
    try:
        variants = [line_.split(",")[0].strip() for line_ in lines]
        inheritances = [line_.split(",")[1].strip() for line_ in lines]
    except IndexError:
        #ToDo error handling
        variants = [line_.strip("\n ") for line_ in lines]
        inheritances = ["de_novo" for _ in lines]

    results_df = score_variants(variants, inheritances)

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
def single(variant, inheritance, family_history):
    results_df = score_variants([variant], [inheritance])
    click.echo(results_df.head(len(results_df)))


if __name__ == "__main__":
    single(["--variant", "NM_006651.3:c.250dup ",
                 "-ih", "de_novo",
                 "-f", "yes"])
    # batch(["--input_file",
    #        "/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/AutoCaSc_core/data/CLI_batch_test_variants.txt"])
    # main(obj={})

# print(AutoCaSc_core("2:46803383:G:A", "de_novo").candidate_score)  # benign
# print(AutoCaSc_core("1:1737948:A:G", "de_novo").candidate_score)  # pathogenic
# print(AutoCaSc_core("X:153040798:C:A", "de_novo").candidate_score)  # benign
# print(AutoCaSc_core("1:3732936:A:G", "de_novo").candidate_score)  # pathogenic --> impact high LOEUF 0.74
# print(AutoCaSc_core("1:7725246:G:A", "de_novo").candidate_score)  # pathogenic
# print(AutoCaSc_core("1:2234850:G:A", "de_novo").candidate_score)  # pathogenic
# print(AutoCaSc_core("1:3755660:A:C", "de_novo").candidate_score)  # pathogenic
# print(AutoCaSc_core("1:1959699:G:A", "other").candidate_score)  # benign, high expressesed, GABA receptor.......
# print(AutoCaSc_core("CBLB:c.1822C>T", "other").candidate_score)
# print(AutoCaSc_core("NM_001113498.2:c.794T>A", "de_novo").candidate_score)
# print(AutoCaSc_core("chr4:175812275:T:C", inheritance="de_novo", assembly="GRCh38").candidate_score)
# print(AutoCaSc_core("3:149619809:xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx:-", "de_novo").candidate_score)  # benign, high expressesed, GABA receptor.......


# homo: NM_001199266.1:c.132_134del