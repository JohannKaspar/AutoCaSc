import os
import pickle
import sys
from pathlib import Path
import random
import tenacity

sys.path.insert(0, str(Path(__file__).parent))

import time
import re
from statistics import mean
from tenacity import retry, stop_after_attempt, wait_exponential, wait_random
import click
import pandas as pd
import requests
from numpy import isnan
from gnomAD import GnomADQuery
from tools import safe_get, get_seq_difference, write_new_api_request

VERSION = 1.0
ROOT_DIR = str(Path(__file__).parent) + "/data/"

gene_scores = pd.read_csv(ROOT_DIR + "all_gene_data.csv")
# This loads the calculated gene scores and gnomad constraints.

class AutoCaSc:
    """AutoCaSc is a tool for quantifying plausibility of gene variants as a cause for neurodevelopmental delay.
    This class is the core. A variant is annotated by initializing the class.
    Class attributes correspond to variant parameters ranging from position to candidate score.
    """
    def __init__(self, variant,
                 inheritance="unknown",
                 other_autocasc_obj=None,
                 other_impact="unknown",
                 other_variant=None,
                 assembly="GRCh37",
                 transcript=None,
                 family_history=False,
                 mode="default",
                 path_to_request_cache_dir=None):
        """This function assigns basic parameters and initializes the scoring process.

        :param variant: the variant including position and alternative sequence
        :param inheritance: segregation/zygosity: de_novo, homo, comphet, other, ad_inherited, x_linked
        :param family_history: are there multiple affected family members? if x_linked affected male family members?
        :param other_impact: impact of the corresponding variant if it is compound heterozygous
        :param assembly: either GRCh37 or GRCh38
        """
        default_none_attributes = ["mutationtaster_converted_rankscore", "comphet_id", "candidate_score",
                                   "ada_score", "rf_score", "maxentscan_ref", "maxentscan_alt", "maxentscan_decrease",
                                   "vcf_string", "ref_seq_vep", "gnomad_frequency", "maf", "impact",
                                   "gerp_rs_rankscore", "transcript", "consequence", "hgvsc_change", "protein",
                                   "hgvsp_change", "polyphen_prediction", "ada_consequence", "maxentscan_consequence",
                                   "hgvsc_change", "pLI", "mis_z", "oe_lof", "oe_lof_lower", "oe_lof_upper", "oe_mis",
                                   "oe_mis_lower", "oe_mis_upper", "oe_mis_interval", "oe_lof_interval",
                                   "sift_converted_rankscore", "mutationtaster_converted_rankscore",
                                   "mutationassessor_rankscore", "mgi_score", "virlof_ar_enrichment", "ref_seq",
                                    "filter_fail_explanation", "factors", "canonical_transcripts", "response_decoded",
                                   "affected_transcripts", "cadd_phred"
                                   ]
        for attribute in default_none_attributes:
            self.__dict__[attribute] = None

        self.variant = variant
        self.inheritance = inheritance
        self.other_autocasc_obj = other_autocasc_obj
        self.other_impact = other_impact
        self.other_variant = other_variant
        self.assembly = assembly
        self.status_code = 200  # initial value is set to 200 = all good
        self.family_history = family_history
        self.num_transcripts = False
        self.transcript = transcript
        self.mode = mode
        self.data_retrieved = False
        self.path_to_request_cache_dir = path_to_request_cache_dir
        self.filter_pass = True
        self.transcript_instances = {}

        if self.assembly == "GRCh37":
            self.server = "http://grch37.rest.ensembl.org"  # API endpoint for GRCh37
        else:
            self.server = "http://rest.ensembl.org"  # API endpoint for GRCh38

        self.explanation_dict = {}
        self.inheritance_score = 0
        self.frequency_score = 0
        self.variant_score = 0
        self.literature_score = 0
        self.gene_constraint_score = 0
        self.candidate_score = 0

        self.check_for_other_variant()
        self.check_variant_format()  # this function is called to check if the entered variant is valid
        if self.variant_format == "incorrect":  # if variant format is valid (not 401) continue
            self.status_code = 401
        elif not self.mode == "web":
            self.retrieve_data()

    def get(self, attribute):
        return self.__dict__.get(attribute)

    def check_for_other_variant(self):
        """Check if a corresponding autocasc objet exists. If so set other_variant to its variant.
        """
        if self.other_autocasc_obj and self.other_variant is None:
            if self.other_autocasc_obj.vcf_string:
                self.other_variant = self.other_autocasc_obj.vcf_string
            else:
                self.other_variant = self.other_autocasc_obj.variant

    def update_inheritance(self, inheritance):
        self.inheritance = inheritance
        for transcript, instance in self.transcript_instances.items():
            instance["inheritance"] = inheritance

    def retrieve_data(self, gnomad=True):
        """Method for collecting data from VEP and gnomad. If used on a local machine it is possible to check,
        if the variant of interest has already been requested and saved before in order to increase analysis speed.
        """
        self.get_vep_data()  # this method call initiates the annotation of the given variant
        self.data_retrieved = True
        if self.status_code == 200 and gnomad is True:
            # 498 = no matching transcript index has been found (e.g. variant is intergenic)
            self.get_gnomad_constraint()
            start = time.time()
            self.get_gnomad_counts()  # gets allele counts in gnomad
            end = time.time()
            print(f"gnomad execution time {round(end - start, 2)}")
            if self.other_variant is not None and self.other_autocasc_obj is None:
                other_instance = AutoCaSc(variant=self.other_variant,
                                          inheritance=self.inheritance,
                                          assembly=self.assembly,
                                          transcript=self.transcript,
                                          mode=self.mode)
                self.other_autocasc_obj = other_instance
                self.status_code = other_instance.status_code

        if self.status_code == 201 and self.variant_format == "hgvs":
            self.hgvs_strand_shift(gnomad=gnomad)


    def hgvs_strand_shift(self, gnomad=True):
        """It sometimes occurs, that VEP returns error codes if a variant is entered in HGVS format. Nevertheless, it
        returns the variant in the correct VCF format. This can then be used to resend the request, this time resulting
        in the correct results.
        """
        test_instance = AutoCaSc(self.vcf_string,
                                 inheritance=self.inheritance,
                                 assembly=self.assembly,
                                 mode=self.mode,
                                 other_variant=self.other_variant,
                                 other_autocasc_obj=self.other_autocasc_obj,
                                 transcript=self.transcript,
                                 family_history=self.family_history,
                                 path_to_request_cache_dir=self.path_to_request_cache_dir
                                 )
        if not test_instance.data_retrieved:
            test_instance.retrieve_data(gnomad=gnomad)
        variant_pre_standshift = self.variant
        if test_instance.status_code == 200:
            self.__dict__ = test_instance.__dict__
            self.variant = variant_pre_standshift


    def create_url(self):
        """This method concatenates all paprts for the request URL to VEP.
        """
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
        """This method pulls gene constraint parameters which are stored together with gene plausibility scores in
        all_gene_data.csv
        """
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
        """This method initiates a GnomADQuery instance in order to retrieve allele counts from gnomad.
        """
        gnomad_variant_result, self.status_code = GnomADQuery(self.vcf_string,
                                                              "variant",
                                                              path_to_request_cache_dir=self.path_to_request_cache_dir).get_gnomad_info()
        if self.status_code == 200:
            self.ac_hom = gnomad_variant_result.get("ac_hom")
            self.allele_count = gnomad_variant_result.get("ac")
            self.n_hemi = gnomad_variant_result.get("ac_hemi")
            self.male_count = gnomad_variant_result.get("male_count") or 0
            self.female_count = gnomad_variant_result.get("female_count") or 0

    def check_variant_format(self):
        """This function checks if the entered variant matches either HGVS or VCF format after doing some formatting.
        """
        if not type(self.variant) == str:
            input = ""
        else:
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
           stop=stop_after_attempt(5),
           wait=wait_exponential(multiplier=1.3, min=0.1, max=3))
    def vep_api_request(self):
        """General function for handling API communication. If there is some error with the returned data, the request
        will be retried a couple of times. Return obtained data as a dict.
        """
        url = self.server + self.ext_url

        # try to pull the requested data and check if request was successful
        r = requests.get(url, headers={"content-type": "application/json"})

        if r.ok:
            return safe_get(r.json(), 0), 200
        else:
            # if some sort of error occurs
            if "matches reference" in str(r.content, "utf-8"):
                # Return None if the given alternative sequence matches the GRCh reference sequence.
                return None, 496
            else:
                # Return None if some sort of different error occurs.
                print(f"VEP ERROR '{r.status_code}: {r.reason}' occured for {self.variant}. Retrying...")
                raise IOError("There has been an issue with a variant.")

    @retry(stop=stop_after_attempt(5),
           wait=wait_random(0.1, 1))
    def open_pickle_file(self):
        """Method for loading stored VEP requests. Useful if the same variant is analysed multiple times.
        """
        self.vep_requests = {}
        with open(f"{self.path_to_request_cache_dir}vep_requests_{self.assembly}", "rb") as vep_requests_file:
            self.vep_requests = pickle.load(vep_requests_file)

    def get_vep_data(self):
        """This function retrieves VEP data either by loading from stored VEP requests or addressing the REST API of
        VEP. It then calls assign_results to select relevant parameters and do some formatting.
        """
        if not self.data_retrieved:
            start = time.time()
            if self.path_to_request_cache_dir is not None:
                try:
                    self.open_pickle_file()
                    if self.vep_requests.get(self.variant):
                        self.response_decoded = self.vep_requests.get(self.variant)
                except (pickle.UnpicklingError, EOFError, tenacity.RetryError):
                    print("could not open vep pickle")

            if self.response_decoded is None:
                self.create_url()
                try:
                    self.response_decoded, self.status_code = self.vep_api_request()
                except IOError:
                    self.response_decoded = None
                    self.status_code = 400

                if self.status_code == 200:
                    if self.path_to_request_cache_dir is not None:
                        new_vep_request = {self.variant: self.response_decoded}
                        write_new_api_request(f"{self.path_to_request_cache_dir}tmp/vep", new_vep_request)

            end = time.time()
            print(f"VEP execution time {round(end - start, 2)}")

        if self.status_code == 200:
            self.get_transcript_index()
            # get the index of the transcript to consider for further annotations
            if self.status_code == 200:
                self.assign_results()


    def assign_results(self, transcript_id=None):
        """This function filters and formats the received annotation results.

        :param self.response_decoded: received annotation data
        :param transcript_index: index of the transcript to work with from transcript_consequences
        """
        if not transcript_id:
            transcript_id = self.transcript
        selected_transcript_consequences = None
        for result in self.response_decoded["transcript_consequences"]:
            if result.get("transcript_id") == transcript_id:
                selected_transcript_consequences = result
                break

        if self.response_decoded.get("vcf_string") is not None:
            self.vcf_string = re.sub(r"-", ":", self.response_decoded.get("vcf_string"))
            # extracts reference sequence from variant id in vcf format
            self.ref_seq_vep = safe_get(self.vcf_string.split(":"), 2)
            self.alt_seq_vep = safe_get(self.vcf_string.split(":"), 3)
            if self.ref_seq is None:
                self.ref_seq = self.response_decoded.get("vcf_string").split("-")[2]
                self.alt_seq = self.response_decoded.get("vcf_string").split("-")[3]
        self.id = self.response_decoded.get("id")
        if self.variant_format == "vcf":
            self.id = re.sub(r"[^a-zA-Z^0-9-]", ":", self.id)

        # checks if reference sequences match and returns result dictionary and status code accordingly
        if any([x is None for x in [self.ref_seq, self.ref_seq_vep, self.alt_seq, self.alt_seq_vep]]):
            pass
        else:
            if not get_seq_difference(self.ref_seq, self.alt_seq) == get_seq_difference(self.ref_seq_vep,
                                                                                            self.alt_seq_vep):
                self.status_code = 201
            else:
                # gets variant gnomad exome frequency and MAF
                self.gnomad_frequency, self.maf = self.get_frequencies(self.response_decoded.get("colocated_variants"))
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
                    self.hgvsc_transcript, self.hgvsc_change = selected_transcript_consequences.get("hgvsc").split(":")

                if selected_transcript_consequences.get("consequence_terms") is not None:
                    self.consequence = re.sub(r"_", " ", safe_get(selected_transcript_consequences.get("consequence_terms"), 0))

                if selected_transcript_consequences.get("hgvsp") is not None:
                    self.protein, self.hgvsp_change = selected_transcript_consequences.get("hgvsp").split(":")

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
            for _variant in colocated_variants:
                frequencies_dict = _variant.get("frequencies")
                if frequencies_dict is not None:
                    if frequencies_dict.get(self.alt_seq) is not None:
                        gnomad_frequency = frequencies_dict.get(self.alt_seq).get("gnomad") or 0
            # getting MAF
            for colocated_dict in colocated_variants:
                if colocated_dict.get("minor_allele_freq") is not None:
                    maf = colocated_dict.get("minor_allele_freq")
            if maf == 0:
                maf = gnomad_frequency
            return gnomad_frequency, maf

    def get_transcript_index(self):
        """This function selects the transcript index to consider for further calculations, based on
        1) the transcript with the most severe consequence
        2) the canonical transcript (if existing if existing and equally severe consequence)
        3) protein coding transcripts
        In case multiple transcripts seem to be equally relevant the first of the list is returned and the number of
        equally important transcripts stored to the class. Depending on the usecase, the other transcripts can be
        scored afterwards.

        :param self.response_decoded: VEP response dictionary with a list of transcripts
        :return: transcript index to consider for further calculations
        """
        impact_severity_dict = {"MODIFIER": 1, "LOW": 1, "MODERATE": 2, "HIGH": 3}
        transcript_df = pd.DataFrame()
        if not self.transcript:
            try:
                for i, transcript in enumerate(self.response_decoded["transcript_consequences"]):
                    # transcript_df.loc[i, "transcript_id"] = i
                    transcript_df.loc[i, "transcript_id"] = transcript.get("transcript_id")
                    transcript_df.loc[i, "impact_level"] = impact_severity_dict.get(transcript.get("impact"))
                    transcript_df.loc[i, "biotype"] = transcript.get("biotype")
                    transcript_df.loc[i, "canonical"] = transcript.get("canonical") or 0

                if "protein_coding" in transcript_df.biotype.unique() and len(transcript_df.biotype.unique()) > 1:
                    transcript_df = transcript_df.loc[transcript_df.biotype == "protein_coding"].reset_index()
                transcript_df = transcript_df.sort_values(by=["impact_level", "canonical"],
                                                          ascending=[False, False],
                                                          ignore_index=True)

                self.affected_transcripts = transcript_df.transcript_id.to_list()
                self.num_transcripts = len(transcript_df)
                transcript_df = transcript_df.loc[
                    transcript_df.impact_level == transcript_df.loc[0, "impact_level"]]
                transcript_df = transcript_df.loc[
                    transcript_df.canonical == transcript_df.loc[0, "canonical"]]
                self.high_priority_transcripts = transcript_df.transcript_id.to_list()
                self.transcript = self.high_priority_transcripts[0]

                if 1 in transcript_df.canonical.to_list():
                    self.canonical_transcripts = transcript_df.loc[transcript_df.canonical == 1, "transcript_id"].to_list()

            except AttributeError:
                self.status_code = 499
            except KeyError:
                self.status_code = 498  # intergenic variant
            except TypeError:
                self.status_code = 497


    def calculate_candidate_score(self, recursively=True):
        """This method calls all the scoring functions and assigns their results to class attributes.
        """
        self.check_for_other_variant()
        self.explanation_dict = {}

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

        if not self.filter_pass:
            self.candidate_score = 0
        else:
            candidate_score_list = [i if i is not None else 0 for i in [
                self.inheritance_score,
                self.gene_constraint_score,
                self.variant_score,
                self.literature_score]]
            self.candidate_score = round(sum(candidate_score_list), 2)

        if self.inheritance == "comphet":
            if self.other_autocasc_obj and recursively:
                self.other_autocasc_obj.calculate_candidate_score(recursively=False)

                self.candidate_score = round(mean([self.candidate_score,
                                                   self.other_autocasc_obj.__dict__.get("candidate_score")]), 2)
                if self.impact == "high" and self.other_autocasc_obj.__dict__.get("impact") == "high":
                    self.candidate_score = round(self.candidate_score + 1., 2)
                    self.impact_score, self.explanation_dict["impact"] = 3, "impact high, biallelic: 3"
                    self.other_autocasc_obj.__dict__["impact_score"], self.other_autocasc_obj.__dict__["explanation_dict"]["impact"] = 3, "impact high, biallelic: 3"

    def rate_inheritance(self):
        """This function scores zygosity/segregation.
        """
        if self.inheritance == "de_novo":
            self.inheritance_score, self.explanation_dict["inheritance"] = 2, "de novo: 2"
        if self.inheritance == "unknown":
            self.inheritance_score, self.explanation_dict["inheritance"] = 0, "other: 0"
        if self.inheritance == "ad_inherited":
            self.inheritance_score, self.explanation_dict["inheritance"] = 0, "inherited autosomal dominant: 0"
        if self.inheritance == "comphet":
            if self.family_history in [True, "yes"]:
                self.inheritance_score, self.explanation_dict["inheritance"] = 3, \
                    "compound heterozygous, multiple affected children: 3"
            else:
                self.inheritance_score, self.explanation_dict["inheritance"] = 1, \
                    "compound heterozygous, one affected child: 1"
        if self.inheritance == "homo":
            if self.family_history in [True, "yes"]:
                self.inheritance_score, self.explanation_dict["inheritance"] = 3, \
                    "homo, multiple affected children: 3"
            else:
                self.inheritance_score, self.explanation_dict["inheritance"] = 2, "homo, one affected child: 2"
        if self.inheritance == "x_linked":
            if self.family_history in [True, "yes"]:
                self.inheritance_score, self.explanation_dict["inheritance"] = 2, \
                    "x_linked and another maternal male relative: 2"
            else:
                self.inheritance_score, self.explanation_dict["inheritance"] = 1, "x_linked and a boy: 1"

    def rate_pli_z(self):
        """This function scores gnomAD constraints (pLI/Z).
        """
        if self.inheritance in ["homo", "comphet"]:
            self.gene_constraint_score, self.explanation_dict["pli_z"] = 0, "recessive: 0"
        else:
            # if impact is not None:
            if self.impact == "high":
                if self.pLI is not None:
                    if self.pLI < 0.5:
                        if self.inheritance == "de_novo":
                            self.gene_constraint_score, self.explanation_dict[
                                "pli_z"] = -2, "de novo LoF & pLI < 0.5: -2"
                        else:
                            self.gene_constraint_score, self.explanation_dict["pli_z"] = 0, "LoF & pLI < 0.5: 0"
                    elif self.pLI < 0.9:
                        if self.inheritance == "de_novo":
                            self.gene_constraint_score, self.explanation_dict[
                                "pli_z"] = 0, "de novo LoF & 0.5 <= pLI < 0.9: 0"
                        else:
                            self.gene_constraint_score, self.explanation_dict[
                                "pli_z"] = 0.5, "LoF & pLI <= 0.5 pLI < 0.9: 0.5"
                    else:
                        self.gene_constraint_score, self.explanation_dict["pli_z"] = 1, "LoF & pLI >= 0.9: 1"
                else:
                    self.gene_constraint_score, self.explanation_dict["pli_z"] = 0, "no data on pLI: 0"
            elif self.impact == "moderate":
                if self.mis_z is not None:
                    if self.mis_z < 0:
                        self.gene_constraint_score, self.explanation_dict["pli_z"] = 0, "missense & Z < 0: 0"
                    elif self.mis_z >= 0 and self.mis_z < 2.5:
                        self.gene_constraint_score, self.explanation_dict[
                            "pli_z"] = 0.5, "missense & 0 <= Z < 2.5: 0.5"
                    elif self.mis_z >= 2.5:
                        self.gene_constraint_score, self.explanation_dict["pli_z"] = 1, "missense & Z >= 2.5: 1"
                else:
                    self.gene_constraint_score, self.explanation_dict["pli_z"] = 0, "no data on Z score: 0"
            else:
                self.gene_constraint_score, self.explanation_dict["pli_z"] = 0, "low impact or no data on impact: 0"

    def rate_impact(self):
        """This function scores the variants predicted impact.
        """
        self.impact_score, self.explanation_dict["impact"] = 0, f"impact {self.impact}: 0"

        if self.impact == "moderate":
            self.impact_score, self.explanation_dict["impact"] = 0, "impact moderate: 0"
        elif self.impact == "high":
            self.impact_score, self.explanation_dict["impact"] = 2, "impact high, heterozygous: 2"
            if self.inheritance == "homo":
                self.impact_score, self.explanation_dict["impact"] = 3, "impact high, biallelic: 3"
            # for comphets, this is checked seperately later

    def rate_in_silico(self):
        """This function scores in silico predictions.
        """
        score_list = []

        self.in_silico_score, self.explanation_dict["in_silico"] = 0, f"{self.impact} and no prediction values: 0"

        if self.sift_converted_rankscore is not None:
            score_list.append(float(self.sift_converted_rankscore))
        if self.mutationtaster_converted_rankscore is not None:
            score_list.append(float(self.mutationtaster_converted_rankscore))
        if self.mutationassessor_rankscore is not None:
            score_list.append(float(self.mutationassessor_rankscore))

        # "average of in silico"
        if score_list:
            self.in_silico_score, self.explanation_dict["in_silico"] = round(mean(score_list), 2),\
                    f"mean in silico score: {round(mean(score_list), 2)}"
        elif self.impact == "high":
            # If there are no in silico scores for the variant and its impact is high, we assume a score of 1 point.
            self.in_silico_score, self.explanation_dict["in_silico"] = 1, "LoF: 1"

        # In case of absent in silico values and a moderate or low impact: consider splice site predictions, if existing
        else:
            affected_splice_counter = 0
            if self.rf_score is not None:
                if self.rf_score >= 0.6:  # cutoff from "Leipzig guidelines for use of Alamut splice site prediction tools"
                    affected_splice_counter += 1
            if self.ada_score is not None:
                if self.ada_score >= 0.6:
                    affected_splice_counter += 1
            if self.maxentscan_decrease is not None:
                if self.maxentscan_decrease <= -0.15:
                    affected_splice_counter += 1

            if (affected_splice_counter == 0) and (self.impact == "moderate"):  # todo talk to Rami about absent insilico in synonymous variant
                self.in_silico_score = 0.5
                self.explanation_dict["in_silico"] = "missense, no available in silico values: 0.5"
            if affected_splice_counter == 1:
                self.in_silico_score = 0.5
                self.explanation_dict["in_silico"] = "splicing affected in one program: 0.5"
            if affected_splice_counter >= 2:
                self.in_silico_score = 1
                self.explanation_dict["in_silico"] = "splicing affected in two or more programs: 1"

    def rate_conservation(self):
        """This function scores the variants conservation.
        """
        if self.gerp_rs_rankscore is not None:
            self.conservation_score, self.explanation_dict["conservation"] = self.gerp_rs_rankscore,\
                   f"GERP++ rankscore {self.gerp_rs_rankscore}: {self.gerp_rs_rankscore}"
        else:  # "estimation in case of absent GERP++"
            if self.impact == "high":
                self.conservation_score, self.explanation_dict["conservation"] = 1, "LoF: 1"
            else:
                self.conservation_score, self.explanation_dict["conservation"] = 0, "no data!: 0"

    def rate_frequency(self):
        """This function scores the variants frequency.
        """
        self.frequency_score, self.explanation_dict["frequency"] = 0, "other: 0"

        if self.inheritance in ["de_novo", "ad_inherited"]:
            if self.maf < 0.000005:
                self.frequency_score, self.explanation_dict["frequency"] = 1,\
                    f"{self.inheritance} & MAF < 0.000005: 1"
            elif self.maf < 0.00002:
                self.frequency_score, self.explanation_dict["frequency"] = 0.5,\
                    f"{self.inheritance} & MAF < 0.00002: 0.5"
            else:
                self.frequency_score, self.explanation_dict["frequency"] = 0, "de novo & MAF > 0.00002: 0"

        if self.inheritance in ["homo", "comphet"]:
            if self.maf < 0.00005:
                self.frequency_score, self.explanation_dict["frequency"] = 1,\
                    "autosomal recessive & MAF < 0.00005: 1"
            elif self.maf < 0.0005:
                self.frequency_score, self.explanation_dict["frequency"] = 0.5,\
                    "autosomal recessive & MAF < 0.0005: 0.5"

        if self.inheritance == "x_linked":
            gnomad_variant_result, _ = GnomADQuery(self.vcf_string, "variant").get_gnomad_info()
            if self.n_hemi <= 1 and self.female_count >= 5:
                self.frequency_score, self.explanation_dict["frequency"] = 2, \
                    "X linked and discrepancy of MAF in gnomAD between males and females (max.1/min.5): 2"
            elif self.n_hemi == 0:
                self.frequency_score = 1
                self.explanation_dict["frequency"] = "X linked 0x hemizygous in gnomAD: 1"
            else:
                self.frequency_score = 0
                self.explanation_dict["frequency"] = f"X linked and {self.n_hemi}x hemizygous in gnomad!"


@click.group(invoke_without_command=True)  # Allow users to call our app without a command
@click.pass_context
@click.option('--verbose', '-v',
              is_flag=True,
              help="Increase output verbosity level")
@click.option("--assembly", "-a",
              default="GRCh37",
              help="Reference assembly to use.")
@click.option("--output_path", "-o",
              help="Output path for csv-file.")
def main(ctx, verbose, assembly, output_path):
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
    ctx.obj["OUTPUT_PATH"] = output_path

@click.pass_context
def score_variants(ctx, variants, inheritances, corresponding_variants, family_histories):
    if ctx.obj:
        verbose = ctx.obj.get('VERBOSE')
        assembly = ctx.obj.get('ASSEMBLY')
        output_path = ctx.obj.get('OUTPUT_PATH')
    else:
        verbose = True
        assembly = "GRCh37"
        output_path = None

    results_df = pd.DataFrame()
    variant_dict = {}
    for _variant, _inheritance, _corresponding_variant, _family_history in \
            zip(variants, inheritances, corresponding_variants, family_histories):
        _instance = AutoCaSc(variant=_variant,
                             inheritance=_inheritance,
                             other_variant=_corresponding_variant,
                             family_history=_family_history,
                             assembly=assembly)
        if _instance.status_code == 200:
            if _instance.other_autocasc_obj:
                if _instance.other_autocasc_obj.status_code == 200:
                    _instance.calculate_candidate_score()
            else:
                _instance.calculate_candidate_score()
        variant_dict[_variant] = _instance

    # create the final result_df containing scoring results
    for _variant, _inheritance in zip(variants, inheritances):
        variant_instance = variant_dict.get(_variant)
        results_df.loc[_variant, "candidate_score"] = variant_instance.candidate_score
        results_df.loc[_variant, "gene_symbol"] = variant_instance.__dict__.get("gene_symbol")
        results_df.loc[_variant, "other_variant"] = variant_instance.other_variant
        if verbose:
            results_df.loc[_variant, "inheritance_score"] = variant_instance.inheritance_score
            results_df.loc[_variant, "gene_score"] = variant_instance.gene_constraint_score
            results_df.loc[_variant, "variant_score"] = variant_instance.variant_score
            results_df.loc[_variant, "literature_score"] = variant_instance.literature_score
        results_df.loc[_variant, "inheritance_mode"] = variant_instance.inheritance
        results_df.loc[_variant, "status_code"] = str(int(variant_instance.status_code))

    results_df.index.names = ["variant"]
    if output_path:
        results_df.to_csv(output_path)
    return results_df


def check_inheritance_input(inheritances):
    for _inheritance in inheritances:
        if not _inheritance in ["comphet", "de_novo", "homo", "ad_inherited", "x_linked", "unknown"]:
            raise IOError(f"Could not interpret inheritance '{_inheritance}'")


@main.command("batch")
@click.option("--input_file", "-i",
              type=click.Path(exists=True),
              help="Path to file with comma-separated columns containing "
                      "Variant | Inheritance | Corresponding comphet variant | Family history. "
                   "For family history use: "
                   "'de_novo' for de novo variants,  'homo' for homozygous variants,  'comphet' for compound "
                   "heterozygous variants, 'ad_inherited' for inherited variants, 'x_linked' for X-linked variants, "
                   "'unknown' for variants where inheritance pattern is unknown")
def batch(input_file):
    lines = open(input_file).readlines()
    lines = [_line.strip("\n ") for _line in lines]
    try:
        variants = [line_.split(",")[0].strip() for line_ in lines]
        inheritances = [line_.split(",")[1].strip() for line_ in lines]
    except IndexError:
        raise IOError("Could not interpret input. Please provide a file with comma-separated columns containing "
                      "Variant | Inheritance | Corresponding variant | Family history")
    check_inheritance_input(inheritances)
    try:
        corresponding_variants = [line_.split(",")[2].strip() for line_ in lines]
        corresponding_variants = [None if _cv == "" else _cv for _cv in corresponding_variants]
    except IndexError:
        corresponding_variants = [None] * len(variants)

    try:
        family_histories = [interpret_family_history(line_.split(",")[3].strip()) for line_ in lines]
    except IndexError:
        family_histories = [False] * len(variants)

    results_df = score_variants(variants=variants,
                                inheritances=inheritances,
                                corresponding_variants=corresponding_variants,
                                family_histories=family_histories)

    pd.set_option('display.max_columns', None)
    click.echo(results_df)

def interpret_family_history(input):
    if isinstance(input, str):
        if "true" in input.lower() or "yes" in input.lower():
            return True
    if input == True:
        return True
    return False

@main.command("single")
@click.option("--variant", "-v",
              required=True,
              help="Variant in either VCF or HGVS format.")
@click.option("--corresponding_variant", "-cv",
              default=None,
              help="Corresponding compound heterozygous variant.")
@click.option("--inheritance", "-ih",
              default="de_novo",
              help="Inheritance mode of variant.")
@click.option("--family_history", "-f",
              default=False,
              help="Is there another family member who could have this variant as well? True / False")
def single(variant, corresponding_variant, inheritance, family_history):
    family_history = interpret_family_history(family_history)
    results_df = score_variants([variant], [inheritance], [corresponding_variant], [family_history])
    pd.set_option('display.max_columns', None)
    click.echo(results_df)


if __name__ == "__main__":
    # single(["-v", "22:45255644:G:T"])
    single(["-v", "9:134736022:G:A", "-ih", "homo"])
    # batch(["-i", "/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/CLI_batch_test_variants.txt"])
    main(obj={})