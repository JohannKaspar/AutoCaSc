import pickle

import pandas as pd
from AutoCaSc_core.tools import add_categories, rank_genes, lin_rank, negative_product
from scipy.stats import spearmanr, mannwhitneyu
import psutil  # used for counting CPUs
import random
from concurrent.futures import ProcessPoolExecutor
import os
import itertools
import time
import wget
import shutil
import gzip
import requests
from io import StringIO
import ray

# ray.init(object_store_memory=30000000000)
random.seed(42)
# ToDo automatisches Downloaden der Dateien wäre cool,
# ToDo pubtator central rohdaten prozessierung

ROOT_DIR = "/home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_maintenance/data/"
# ROOT_DIR = "/Users/johannkaspar/OneDrive/Promotion/AutoCaSc_project_folder/AutoCaSc_maintenance/data/"

sysid_primary = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv",
                            usecols=["Entrez id", "Ensembl id"])
sysid_primary.columns = ["entrez_id", "ensemble_id"]
sysid_candidates = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv",
                               usecols=["Entrez id", "Ensembl id"])
sysid_candidates.columns = ["entrez_id", "ensemble_id"]
princeton_negative = pd.read_csv(ROOT_DIR + "ASD_translated_to_ensembl.csv")["entrez_id"].to_list()
sysid_primary_gene_symbols = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv", usecols=["Gene symbol"])[
    "Gene symbol"].to_list()
sysid_candidates_gene_symbols = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv", usecols=["Gene symbol"])[
    "Gene symbol"].to_list()

morbid_gene_symbols_list = list(set(pd.read_csv("/home/johann/AutoCaSc/data/pubtator_central/MorbidGenes-Panel"
                                                "-v5_2020-08-26_for_varvis.csv", header=None).iloc[:, 0].to_list()) -
                                set(sysid_primary_gene_symbols + sysid_candidates_gene_symbols))
all_genes_df = pd.read_csv(ROOT_DIR + "hgnc_protein_coding.tsv", index_col=False,
                           usecols=["entrez_id", "gene_symbol"], sep="\t",
                           dtype={"entrez_id": "Int32", "gene_symbol": str})
morbid_genes = all_genes_df.loc[all_genes_df.gene_symbol.isin(morbid_gene_symbols_list)]["entrez_id"].to_list()
# negative_gene_list = princeton_negative
negative_gene_list = random.sample(morbid_genes, round(0.8 * len(morbid_genes)))
del all_genes_df, morbid_gene_symbols_list

def HGNC():
    r = requests.get(
        "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_pub_refseq_ids&col=gd_pub_eg_id&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit")
    if r.ok:
        hgnc_table = pd.read_csv(StringIO(r.text),
                                 sep="\t",
                                 header=0,
                                 names=["hgnc_id", "gene_symbol", "refseq_id", "entrez_id", "ensemble_id"])
        if os.path.isfile("/home/johann/AutoCaSc_core/data/hgnc_protein_coding.tsv"):
            os.rename("/home/johann/AutoCaSc_core/data/hgnc_protein_coding.tsv",
                      "/home/johann/AutoCaSc_core/data/hgnc_protein_coding.tsv_old")
        hgnc_table.to_csv("/home/johann/AutoCaSc_core/data/hgnc_protein_coding.tsv",
                          sep="\t",
                          index=False)
    else:
        print("Could not load HGNC data! Using old data instead.")

def fuse_data():
    all_genes = pd.read_csv("/home/johann/AutoCaSc_core/data/hgnc_protein_coding.tsv",
                            sep="\t",
                            usecols=["entrez_id", "ensemble_id"])

    gtex = pd.read_csv("/home/johann/AutoCaSc_core/data/gtex/gene_scores.csv", usecols=["ensemble_id", "rank_score"])
    gtex.columns = ["ensemble_id", "gtex_score"]

    denovo = pd.read_csv("/home/johann/AutoCaSc_core/data/psymukb/all_genes_denovo_NDD.csv",
                         usecols=["entrez_id", "rank_score"])
    denovo.columns = ["entrez_id", "denovo_rank_score"]

    disgenet = pd.read_csv("/home/johann/AutoCaSc_core/data/disgenet/disgenet_gene_scores.csv",
                           usecols=["entrez_id", "gene_score"])
    disgenet.columns = ["entrez_id", "disgenet_score"]

    # mgi = pd.read_csv("/home/johann/AutoCaSc_core/data/mgi/gene_annoations.csv")
    # mgi.columns = ["entrez_id", "mgi_phenotypes", "mgi_neuro_behavioral"]
    # mgi = mgi[["entrez_id", "mgi_neuro_behavioral"]]

    mgi = pd.read_csv("/home/johann/AutoCaSc_core/data/mgi/mgi_gene_scores.csv", usecols=["entrez_id", "rank_score"])
    mgi.columns = ["entrez_id", "mgi_score"]

    pubtator = pd.read_csv("/home/johann/AutoCaSc_core/data/pubtator_central/gene_scores.csv",
                           usecols=["entrez_id", "gene_score"])
    pubtator.columns = ["entrez_id", "pubtator_score"]

    string = pd.read_csv("/home/johann/AutoCaSc_core/data/string/gene_scores.csv", usecols=["gene_id", "gene_score"])
    string.columns = ["ensemble_id", "string_score"]

    gnomad = pd.read_csv("/home/johann/AutoCaSc_core/data/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.csv", sep="\t")
    gnomad = gnomad[["gene_id", "pLI", "oe_lof",
                     "oe_lof_lower", "oe_lof_upper",
                     "oe_mis", "oe_mis_lower",
                     "oe_mis_upper", "mis_z"]]
    gnomad.columns = ["ensemble_id", "pLI", "oe_lof",
                      "oe_lof_lower", "oe_lof_upper",
                      "oe_mis", "oe_mis_lower",
                      "oe_mis_upper", "mis_z"]

    all_data = all_genes.merge(gtex, on="ensemble_id", how="outer")
    all_data = all_data.merge(denovo, on="entrez_id", how="outer")
    all_data = all_data.merge(disgenet, on="entrez_id", how="outer")
    all_data = all_data.merge(mgi, on="entrez_id", how="outer")
    all_data = all_data.merge(pubtator, on="entrez_id", how="outer")
    all_data = all_data.merge(string, on="ensemble_id", how="outer")

    all_data = all_data.fillna(0)

    # careful with filling nans, o/e = 0 would bias results
    all_data = all_data.merge(gnomad, on="ensemble_id", how="outer")

    all_data = add_categories(all_data, "entrez_id", "entrez")

    # all_data = all_data.dropna(subset=["mgi_score", "string_score", "pubtator_score", "disgenet_score", "gtex_score", "denovo_rank_score"])

    all_data[["mgi_score", "string_score", "pubtator_score", "disgenet_score", "gtex_score", "denovo_rank_score"]] = \
        all_data[
            ["mgi_score", "string_score", "pubtator_score", "disgenet_score", "gtex_score",
             "denovo_rank_score"]].fillna(0)

    all_data[
        "gene_sum"] = all_data.gtex_score + all_data.denovo_rank_score + all_data.disgenet_score + all_data.mgi_score + all_data.pubtator_score + all_data.string_score
    all_data["weighted_score"] = sum(
        [spearmanr(all_data[parameter], all_data.sys_primary)[0] * all_data[parameter] for parameter in
         ["mgi_score", "string_score", "pubtator_score", "disgenet_score", "gtex_score", "denovo_rank_score"]])
    all_data.weighted_score = all_data.weighted_score / all_data.weighted_score.max()

    all_data = all_data.sort_values(by="weighted_score", ascending=False).drop_duplicates(
        subset=["entrez_id", "ensemble_id"], keep="first")
    all_data.to_csv("/home/johann/AutoCaSc_core/data/all_gene_data.csv", index=False)


def update_gnomad_data():
    all_genes = pd.read_csv("/home/johann/AutoCaSc_core/data/protein_gene_table.csv")
    gnomad_data = pd.read_csv("/home/johann/AutoCaSc_core/data/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.csv",
                              sep="\t",
                              usecols=["transcript", "pLI", "oe_lof", "oe_lof_lower", "oe_lof_upper",
                                       "oe_mis", "oe_mis_lower", "oe_mis_upper", "mis_z"])

    df = all_genes.merge(gnomad_data, left_on="transcript_id", right_on="transcript", how="left").sort_values(by="pLI",
                                                                                                              ascending=False)
    df = df.drop_duplicates(subset=["entrez_id"], keep="first")
    df = df.loc[~df.entrez_id.isnull()]

    df.to_csv("/home/johann/AutoCaSc_core/data/gnomad/gnomad_data.csv", index=False)


class MGI:
    """
    This class takes in two files from MGI, containing a translation from entrez to MGI identifiers and all
    phenotypes attributed to a given genotype.
    """

    def __init__(self, n_cores=4, download=False):
        if download:
            if os.path.isfile("/home/johann/AutoCaSc_core/data/mgi/HMD_HumanPhenotype.rpt"):
                os.rename("/home/johann/AutoCaSc_core/data/mgi/HMD_HumanPhenotype.rpt",
                          "/home/johann/AutoCaSc_core/data/mgi/HMD_HumanPhenotype_old.rpt")
            wget.download("http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt",
                               "/home/johann/AutoCaSc_core/data/mgi/HMD_HumanPhenotype.rpt")
        self.entrez_mgi = pd.read_csv("/home/johann/AutoCaSc_core/data/mgi/HMD_HumanPhenotype.rpt", sep="\t",
                                      usecols=[0, 1, 5])
        self.entrez_mgi.columns = ["gene_symbol", "entrez_id", "mgi_id"]
        self.entrez_mgi.mgi_id = self.entrez_mgi.mgi_id.str.strip()
        self.mgi_mp = pd.read_csv("/home/johann/AutoCaSc_core/data/mgi/MGI_PhenoGenoMP.rpt.txt", sep="\t", usecols=[3, 5])
        self.mgi_mp.columns = ["phenotype", "mgi_id"]
        self.df = self.entrez_mgi.merge(self.mgi_mp, on="mgi_id")
        self.n_cores = n_cores
        self.pval_df = None

        self.generate_empricial_p_values()
        self.calculate_gene_scores()

    def generate_empricial_p_values(self):
        self.observed_df = self.df.loc[self.df.entrez_id.isin(sysid_primary.entrez_id)].groupby(
            "phenotype").size().reset_index(
            name="count").sort_values(by="count", ascending=False)

        # term_count_pos.to_csv("/home/johann/AutoCaSc_core/data/mgi/mp_count_observed.csv", index=False)

        """
        This part performs bootstraps on the data to identify phenotypes specific to NDD-genes.
        """

        # sufficient for bonferroni correction:
        n_experiments = 20 * len(self.observed_df) + 1

        processes_per_chunk = round(n_experiments / self.n_cores)
        process_id_lists = [list(range(n * processes_per_chunk, (n + 1) * processes_per_chunk)) for n in
                            range(self.n_cores)]

        global mgi_df
        mgi_df = self.df
        with ProcessPoolExecutor() as executor:
            self.bootstrap_df = pd.DataFrame(self.observed_df, columns=["phenotype"])
            results = executor.map(self.mgi_bootstrap_helper, process_id_lists)
            for chunk in results:
                self.bootstrap_df = self.bootstrap_df.merge(chunk, on="phenotype")
            self.bootstrap_df.columns = ["phenotype"] + [f"max_count_{i + 1}" for i in range(self.n_cores)]

        """
        This part calculates empirical p-values by dividing the observed experiments with an equal or higher count
        by the number of all experiments done.
        """
        row_chunks = [list(range(round(n * len(self.bootstrap_df) / self.n_cores),
                                 round((n + 1) * len(self.bootstrap_df) / self.n_cores)))
                      for n in range(self.n_cores)]
        with ProcessPoolExecutor() as executor:
            self.pval_df = pd.DataFrame()
            results = executor.map(self.empirical_p_helper, row_chunks)
            for chunk in results:
                self.pval_df = pd.concat([self.pval_df, chunk])

        self.pval_df.to_csv(ROOT_DIR + "mgi/phenotype_p_values.csv", index=False)

    def calculate_gene_scores(self, correction_exp=0.5, rank_exp=2):
        """
        This part calculates gene scores by counting the number of times a gene has been associated with a significantly
        relevant phenotype and divides it by the square root of the total number of times a gene has been associated with
        some phenotype.
        """
        if self.pval_df is None:
            self.pval_df = pd.read_csv(ROOT_DIR + "mgi/phenotype_p_values.csv", index_col=False)
        relevant_terms = self.pval_df.loc[self.pval_df.p_empirical < .05 / (len(self.pval_df) - 1)][
            "phenotype"].to_list()
        df_filtered = self.df.loc[self.df.phenotype.isin(relevant_terms)]

        gene_counts = df_filtered.groupby("entrez_id").size().reset_index(name="relevant_count")
        gene_df = self.df.groupby("entrez_id").size().reset_index(name="all_count").merge(gene_counts, on="entrez_id")
        gene_df["gene_score"] = gene_df.relevant_count / (gene_df.all_count ** correction_exp)
        gene_df = gene_df[["entrez_id", "gene_score"]]

        """
        Finally gene scores are calculated by ranking genes and taking their rank score to the power of two.
        """
        self.gene_scores = rank_genes(gene_df, "gene_score", rank_exp)

        self.gene_scores.to_csv(ROOT_DIR + "mgi/mgi_gene_scores.csv", index=False)

        return self.gene_scores

    def mgi_bootstrap_helper(self, process_ids):
        df_mgi_subset = mgi_df.loc[mgi_df.phenotype.isin(self.observed_df.phenotype)][["entrez_id", "phenotype"]]
        df_mgi_subset.entrez_id = pd.to_numeric(df_mgi_subset.entrez_id, downcast="unsigned")
        chunk = pd.DataFrame(self.observed_df, columns=["phenotype"])
        for n, process in enumerate(process_ids):
            print(f"{n + 1} of {len(process_ids)}")
            random.seed(process)
            pos_list = random.sample(set(self.df.entrez_id.unique()) - set(
                sysid_candidates.entrez_id.to_list() + sysid_primary.entrez_id.to_list()),
                                     len(sysid_primary.entrez_id))

            df_pos = df_mgi_subset.loc[df_mgi_subset.entrez_id.isin(pos_list)]
            term_count_pos = df_pos.groupby("phenotype").size().reset_index(name=f"count_{process + 1}")

            chunk = chunk.merge(term_count_pos[["phenotype", f"count_{process + 1}"]], on="phenotype", how="left")

        highest_count = pd.DataFrame(chunk.phenotype, columns=["phenotype", "highest_count"])
        chunk = chunk.fillna(0)
        # for phenotype in chunk.phenotype:
        # highest_count.loc[highest_count.phenotype == phenotype, "highest_count"] = chunk.iloc[:, 1:].values.max()
        highest_count["highest_count"] = chunk.max(axis=1)
        return highest_count

    def empirical_p_helper(self, rows):
        bootstraps = [f"max_count_{i + 1}" for i in range(len(self.bootstrap_df.columns) - 1)]
        chunk = self.bootstrap_df.iloc[rows, :].reset_index(drop=True)
        observed_chunk = self.observed_df.loc[self.observed_df.phenotype.isin(chunk.phenotype)]

        for x, phenotype in enumerate(chunk.phenotype):
            print(f"{x} of {len(chunk)}")
            chunk_twisted = chunk.loc[chunk.phenotype == phenotype, bootstraps].T
            chunk_twisted.columns = ["records"]
            count_equal_or_greater = len(chunk_twisted.loc[chunk_twisted.records >= observed_chunk.loc[
                observed_chunk.phenotype == phenotype, "count"].values[0]])
            chunk.loc[chunk.phenotype == phenotype, "p_empirical"] = 1. * count_equal_or_greater / (
                    20 * len(self.observed_df) + 1)

        return chunk[["phenotype", "p_empirical"]]


def evaluate_parameters(df, nomenclature, parameter, version):
    if nomenclature == "entrez_id":
        df["sys_primary"] = df.entrez_id.isin(sysid_primary.entrez_id).astype(int)
        df["sys_candidate"] = df.entrez_id.isin(sysid_candidates.entrez_id).astype(int)
        df["sys"] = df.entrez_id.isin(sysid_primary.entrez_id.to_list() + sysid_candidates.entrez_id.to_list()).astype(
            int)
    if nomenclature == "ensemble_id":
        df["sys_primary"] = df.ensemble_id.isin(sysid_primary.ensemble_id).astype(int)
        df["sys_candidate"] = df.ensemble_id.isin(sysid_candidates.ensemble_id).astype(int)
        df["sys"] = df.ensemble_id.isin(sysid_primary.ensemble_id + sysid_candidates.ensemble_id).astype(int)
    print(
        f"{version}: {mannwhitneyu(df.loc[df.sys_candidate == 1][parameter].to_list(), df.loc[df.sys == 0][parameter].to_list(), alternative='greater')[1]}")


class StringDB:
    """This class takes in a file from StrinDB containing the interactions strengths between protein pairs,
    maps the proteins to genes and calculates a gene score, that is dependent on the amount of interactions with
    known NDD genes and the overall number of interactions.
    """

    def __init__(self, n_cores=4, download=False):
        self.n_cores = n_cores
        # if download:
        #     r = requests.get("https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz", stream=True)
        #     if r.status_code == 200:
        #         if os.path.isfile("/home/johann/AutoCaSc_core/data/string/9606.protein.links.detailed.v11.0.txt"):
        #             os.rename("/home/johann/AutoCaSc_core/data/string/9606.protein.links.detailed.v11.0.txt",
        #                       "/home/johann/AutoCaSc_core/data/string/9606.protein.links.detailed.v11.0.txt_old")
        #         with open("/home/johann/AutoCaSc_core/data/string/9606.protein.links.detailed.v11.0.txt", 'wb') as f:
        #             r.raw.decode_content = True  # just in case transport encoding was applied
        #             gzip_file = gzip.GzipFile(fileobj=r.raw)
        #             shutil.copyfileobj(gzip_file, f)
        #     else:
        #         print("Could not load StrinDB data! Using old data instead.")
        self.string_data = pd.read_csv("/home/johann/AutoCaSc_core/data/string/9606.protein.links.detailed.v11.0.txt",
                                       sep=" ",
                                       usecols=["protein1", "protein2", "combined_score"])

        r = requests.get(
            'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_peptide_id" /></Dataset></Query>')
        if r.ok:
            self.pg_table = pd.read_csv(StringIO(r.text), sep="\t", names=["ensemble_id", "protein_id"])
        else:
            print("Could not load biomart information from Ensebml BioMart! Using old data instead.")
            self.pg_table = pd.read_csv("/home/johann/AutoCaSc_core/data/protein_gene_table.csv",
                                        usecols=["ensemble_id", "protein_id"])

        # cleaning the protein IDs by deleting organism identifier
        self.string_data["protein1"] = self.string_data.protein1.apply(lambda x: x.split(".")[1])
        self.string_data["protein2"] = self.string_data.protein2.apply(lambda x: x.split(".")[1])

        # mapping gene IDs to protein IDs
        self.string_data = self.string_data.merge(self.pg_table, left_on="protein1", right_on="protein_id")
        self.string_data = self.string_data.merge(self.pg_table, left_on="protein2", right_on="protein_id",
                                                  suffixes=["_1", "_2"])

        self.string_data = self.string_data[["ensemble_id_1", "ensemble_id_2", "combined_score",
                                             "protein1", "protein2"]]
        self.string_data.columns = ["gene_1", "gene_2", "score", "protein_1", "protein_2"]

        # self.string_data["score"] = self.string_data["score"] / 999.
        # in case of two genes being mapped to each other take the greatest interaction score in case are differences
        self.string_data = self.string_data.sort_values(by="score", ascending=False)
        self.string_data = self.string_data.drop_duplicates(["gene_1", "gene_2"], "first")

        # make columns to identify if the second gene is a known NDD gene or not
        # self.string_data["sys_primary_2"] = self.string_data.gene_2.isin(sysid_primary.ensemble_id).astype(int)
        self.string_data["sys_2"] = self.string_data.gene_2.isin(sysid_primary.ensemble_id.to_list() +
                                                                 sysid_candidates.ensemble_id.to_list()).astype(int)

        # split data into chunks for parallel processing
        gene_chunks = [self.string_data.gene_1.unique()[round(i * self.string_data.gene_1.nunique() / self.n_cores):
                                                        round((i + 1) * self.string_data.gene_1.nunique() / \
                                                              self.n_cores)] for i in range(self.n_cores)]

        with ProcessPoolExecutor() as executor:
            results = executor.map(self.gene_score_helper, gene_chunks)
            gene_scores_df = pd.DataFrame()
            for chunk in results:
                gene_scores_df = pd.concat([gene_scores_df, chunk], ignore_index=True)

        # sum_corrected_exp_07 was the best correction parameter of all investigated ones
        gene_scores_df = lin_rank(gene_scores_df, "sum_corrected_exp_07")
        gene_scores_df.to_csv("/home/johann/AutoCaSc_core/data/string/gene_scores.csv", index=False)

    def gene_score_helper(self, gene_chunk):
        """This is the helper function that is executed by the parallel processes.
        """
        df = pd.DataFrame()
        for n, gene_1 in enumerate(gene_chunk):
            print(f"processing {n} of {len(gene_chunk)}")

            df_chunk = self.string_data.loc[self.string_data.gene_1 == gene_1]
            sys_sum = df_chunk.loc[df_chunk["sys_2"] == True]["score"].sum(axis=0)
            other_sum = df_chunk.loc[df_chunk.sys_2 == False]["score"].sum(axis=0)
            sum_corrected_exp_07 = sys_sum / (len(df_chunk) ** 0.7)

            is_sys_primary = gene_1 in sysid_primary.ensemble_id.to_list()
            is_sys = is_sys_primary or gene_1 in sysid_candidates.ensemble_id.to_list()

            df = pd.concat([df, pd.DataFrame([[gene_1, sys_sum, other_sum, sum_corrected_exp_07,
                                               is_sys_primary, is_sys]],
                                             columns=["gene_id", "sys_sum", "other_sum", "sum_corrected_exp_07",
                                                      "sys_primary", "sys"])])
        return df


class PsyMuKB:
    def __init__(self):
        df = pd.DataFrame()
        # if download:
        #     r = requests.get("http://47.97.198.109:8080/PsyMuKB%20v1.5/DNMs/MasterFile_allDNMs-codingLoc_v1.5.csv")
        #     if r.ok:
        #         df = pd.read_csv(StringIO(r.text),
        #                          usecols=["EntrezID", "PrimaryPhenotype", "PubmedID", "ExonicFunc.refGene"])
        #     else:
        #         print("Could not load StrinDB data! Using old data instead.")
        if df.empty:
            df = pd.read_csv(ROOT_DIR + "psymukb/MasterFile_allDNMs-codingLoc_v1.5.csv",
                             usecols=["EntrezID", "PrimaryPhenotype", "PubmedID", "ExonicFunc.refGene"])

        all_genes = pd.read_csv("/home/johann/AutoCaSc_core/data/hgnc_protein_coding.tsv",
                                sep="\t",
                                usecols=["entrez_id"]).drop_duplicates(subset=["entrez_id"])

        df.columns = ["gene_id", "consequence", "pmid", "phenotype"]
        df = df.loc[df.phenotype.isin(['Uncharacterized (Mixed healthy control)', 'Sibling Control',
                                       'Congenital Heart Disease (CHD)', 'Amyotrophic Lateral Sclerosis (ALS)',
                                       'Congenital Diaphragmatic Hernia (CDH)', 'Early-onset Alzheimer Disorder (eoAD)',
                                       'Early-onset Parkinson Disorder (eoPD)', 'Obsessive-Compulsive Disorder (OCD)',
                                       'Early-onset High Myopia (eoHM)', 'Fetal non-Preterm birth (non-PTB)',
                                       'Fetal preterm birth (PTB)',
                                       'Mesial Temporal Lobe Epilepsy with Hippocampal Sclerosis (MTLE-HS)']) == False]
        df = df.loc[df.consequence.isin(["synonymous SNV", "unknown"]) == False]

        while sum(df.gene_id.str.contains(";")) != 0:
            for id in df.gene_id.unique():
                id = str(id)
                if ";" in id:
                    n_id = len(df.loc[df.gene_id == id])

                    df = pd.concat([df, df.loc[df.gene_id == id]], ignore_index=True)
                    df.iloc[- n_id:].loc[:, "gene_id"] = id.split(";")[1]
                    df.loc[df.gene_id == id, "gene_id"] = id.split(";")[0]

        denovo_count = df.groupby("gene_id").size().reset_index(name="denovo_count")
        denovo_count.gene_id = denovo_count.loc[denovo_count.gene_id != "NA"].gene_id.astype(int)
        denovo_count = all_genes.merge(denovo_count, left_on="entrez_id", right_on="gene_id", how="left").fillna(0)
        denovo_count = rank_genes(denovo_count, "denovo_count")
        denovo_count.to_csv("/home/johann/AutoCaSc_core/data/psymukb/all_genes_denovo_NDD.csv", index=False)


class GTEx:
    def __init__(self):
        gtex_data = pd.read_csv(ROOT_DIR + "gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                                sep="\t",
                                skiprows=2)

        gtex_data["ensemble_id"] = gtex_data.Name.apply(lambda x: x.split(".")[0])
        gtex_data = gtex_data.drop(columns=["Name"])

        translation_df = pd.read_csv("/home/johann/AutoCaSc_core/data/hgnc_protein_coding.tsv", sep="\t")
        gtex_data = gtex_data.merge(translation_df[["entrez_id", "ensemble_id", "hgnc_id"]], on="ensemble_id")

        columns = gtex_data.columns
        brain_columns = [x for x in columns if "Brain" in x]

        gtex_data["brain_sum"] = gtex_data[brain_columns].sum(axis=1)

        gtex_data = gtex_data[["ensemble_id", "entrez_id", "brain_sum"]]

        gtex_data = rank_genes(gtex_data, "brain_sum", 4)
        gtex_data.to_csv("/home/johann/AutoCaSc_core/data/gtex/gene_scores.csv", index=False)


class Disgenet:
    def __init__(self, n_cores=4, download=False):
        if download:
            r = requests.get("https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz",
                             stream=True)
            if r.status_code == 200:
                if os.path.isfile("/home/johann/AutoCaSc_core/data/disgenet/all_gene_disease_associations.tsv"):
                    os.rename("/home/johann/AutoCaSc_core/data/disgenet/all_gene_disease_associations.tsv",
                              "/home/johann/AutoCaSc_core/data/disgenet/all_gene_disease_associations.tsv_old")
                with open("/home/johann/AutoCaSc_core/data/disgenet/all_gene_disease_associations.tsv", 'wb') as f:
                    r.raw.decode_content = True  # just in case transport encoding was applied
                    gzip_file = gzip.GzipFile(fileobj=r.raw)
                    shutil.copyfileobj(gzip_file, f)
            else:
                print("Could not load Disgenet data! Using old data instead.")

        self.n_cores = n_cores
        self.df = pd.read_csv("/home/johann/AutoCaSc_core/data/disgenet/all_gene_disease_associations.tsv", sep="\t",
                              usecols=["geneId", "diseaseName", "diseaseClass", "score"])
        self.df = self.df.loc[self.df.diseaseClass.str.contains("C10|F03", na=False)]
        self.all_genes = list(self.df.geneId.unique())

        self.observe_term_ratio()
        self.bootstrap()
        self.score_genes()

    def observe_term_ratio(self):
        df_pos = self.df.loc[self.df.geneId.isin(sysid_primary.entrez_id.to_list())]
        observed_count_pos = df_pos.groupby("diseaseName").size().reset_index(name="count").sort_values(by="count",
                                                                                                        ascending=False)

        df_neg = self.df.loc[self.df.geneId.isin(princeton_negative)]
        observed_count_neg = df_neg.groupby("diseaseName").size().reset_index(name="count").sort_values(by="count",
                                                                                                        ascending=False)

        self.observed_df = observed_count_pos.merge(observed_count_neg, on="diseaseName", how="left",
                                                    suffixes=["_pos", "_neg"])
        self.observed_df = self.observed_df.fillna(0)

        self.observed_df["ratio"] = self.observed_df.count_pos / (self.observed_df.count_neg + 1)

    def bootstrap(self):
        self.observed_terms = self.observed_df.diseaseName.to_list()
        # sufficient for bonferroni correction:
        n_experiments = 20 * len(self.observed_terms) + 1
        # n_experiments = len(self.observed_terms) + 1

        processes_per_chunk = round(n_experiments / self.n_cores)
        process_id_lists = [list(range(n * processes_per_chunk, (n + 1) * processes_per_chunk)) for n in
                            range(self.n_cores)]

        self.df = self.df.loc[self.df.diseaseName.isin(self.observed_terms)]

        with ProcessPoolExecutor() as executor:
            self.bootstrap_df = pd.DataFrame(self.observed_terms, columns=["diseaseName"])
            results = executor.map(self.bootstrap_helper, process_id_lists)
            for chunk in results:
                self.bootstrap_df = self.bootstrap_df.merge(chunk, on="diseaseName")
            self.bootstrap_df.columns = ["diseaseName"] + [f"max_ratio_{i + 1}" for i in range(self.n_cores)]

        row_chunks = [list(range(round(n * len(self.bootstrap_df) / self.n_cores),
                                 round((n + 1) * len(self.bootstrap_df) / self.n_cores)))
                      for n in range(self.n_cores)]

        with ProcessPoolExecutor() as executor:
            self.pval_df = pd.DataFrame()
            results = executor.map(self.empirical_p_helper, row_chunks)
            for chunk in results:
                self.pval_df = pd.concat([self.pval_df, chunk])

        self.pval_df = self.pval_df.merge(self.observed_df[["diseaseName", "ratio"]], on="diseaseName", how="left")
        self.pval_df.to_csv(ROOT_DIR + "disgenet/disease_p_values.csv", index=False)

    def bootstrap_helper(self, process_ids):
        chunk = pd.DataFrame(self.observed_terms, columns=["diseaseName"])
        for n, process in enumerate(process_ids):
            print(f"{n + 1} of {len(process_ids)}")
            random.seed(process)
            pos_list = random.sample(self.all_genes, len(sysid_primary))
            neg_list = random.sample(set(self.all_genes) - set(pos_list), len(princeton_negative))

            df_pos = self.df.loc[self.df.geneId.isin(pos_list)]
            term_count_pos = df_pos.groupby("diseaseName").size().reset_index(name="count").sort_values(by="count",
                                                                                                        ascending=False)
            df_neg = self.df.loc[self.df.geneId.isin(neg_list)]
            term_count_neg = df_neg.groupby("diseaseName").size().reset_index(name="count").sort_values(by="count",
                                                                                                        ascending=False)
            term_count_neg = term_count_neg.fillna(0)

            term_count = term_count_pos.merge(term_count_neg, on="diseaseName", how="left", suffixes=["_pos", "_neg"])
            term_count[f"ratio_{process + 1}"] = term_count.count_pos / (term_count.count_neg + 1)

            chunk = chunk.merge(term_count[["diseaseName", f"ratio_{process + 1}"]], on="diseaseName", how="left")

        highest_count = pd.DataFrame(chunk.diseaseName, columns=["diseaseName", "highest_count"])
        highest_count["highest_count"] = chunk.fillna(0).max(axis=1)

        return highest_count

    def empirical_p_helper(self, rows):
        bootstraps = [f"max_ratio_{i + 1}" for i in range(len(self.bootstrap_df.columns) - 1)]
        chunk = self.bootstrap_df.iloc[rows, :].reset_index(drop=True)

        # ratio of "chunk" disease terms in observed ratio of SysID vs. Princeton
        observed_chunk = self.observed_df.loc[self.observed_df.diseaseName.isin(chunk.diseaseName)]

        for x, disease_term in enumerate(chunk.diseaseName):
            print(f"{x} of {len(chunk)}")
            chunk_twisted = chunk.loc[chunk.diseaseName == disease_term, bootstraps].T
            chunk_twisted.columns = ["records"]
            count_equal_or_greater = len(chunk_twisted.loc[chunk_twisted.records >= observed_chunk.loc[
                observed_chunk.diseaseName == disease_term, "ratio"].values[0]])
            chunk.loc[chunk.diseaseName == disease_term, "p_empirical"] = 1. * count_equal_or_greater / (
                    20 * len(self.observed_df) + 1)

        return chunk[["diseaseName", "p_empirical"]]

    def score_genes(self):
        relevant_terms = \
            self.pval_df.loc[self.pval_df.p_empirical < .05 / (len(self.pval_df) - 1)].loc[self.pval_df.ratio >= 1][
                "diseaseName"].to_list()
        df_filtered = self.df.loc[self.df.diseaseName.isin(relevant_terms)]

        gene_scores = df_filtered.groupby("geneId")["score"].apply(negative_product).reset_index(name="gene_score")
        gene_scores.columns = ["entrez_id", "gene_score"]

        gene_scores.to_csv(ROOT_DIR + "disgenet/disgenet_gene_scores.csv", index=False)

def matrix_helper(mapped_parameter):
    """This is the helper function for parallel processing of matrix chunks.
    """
    gene_list_chunk, process_id, chunk_name = mapped_parameter
    helper_chunk = pd.read_csv(ROOT_DIR + f"pubtator_central/matrix/temp/df_chunks/{chunk_name}", index_col=False)
    pivoted_df = pd.pivot_table(helper_chunk.loc[helper_chunk.gene_id.isin(gene_list_chunk)], index=["gene_id"],
                                columns=["mesh_term"], values="pmid", aggfunc=list)
    pivoted_df.to_csv(ROOT_DIR + f"pubtator_central/matrix/temp/matrix_chunks/matrix_split_chunk_{process_id}.csv")


def converter(x):
    if "[" in x:
        return [int(element) for element in x.strip("[]").replace("'", "").split(", ")]
    else:
        return []
    # return [int(element) for element in x.strip("[]").replace("'", "").split(", ")] if "[" in x else []


def process_gene_df_helper(param_tuple):
    """This is the helper function for parallel processing of gene2pubtatorcentral. It takes two parameters, the
    list of genes to process and the index of the process, packed as a tuple.
    :param param_tuple: tuple-packed list of genes and process index
    """
    gene_list, process_id = param_tuple
    # this prevents memory overflow as later processes wait longer, so earlier ones can already filter for relevant
    time.sleep(process_id * 15)
    print("initializing process")
    df_chunk = pd.read_csv(ROOT_DIR + "pubtator_central/gene2pubtatorcentral.csv",
                           index_col=False,
                           engine="c")
    gene_string_list = [str(x) for x in gene_list]
    df_chunk = df_chunk.loc[df_chunk.gene_id.isin(gene_string_list)]  # filters for genes to process
    df_to_concat = pd.DataFrame(columns=df_chunk.columns)

    # issues are lines, where the gene_id contains non-numerical characters
    df_issues = df_chunk[df_chunk["gene_id"].str.contains(";|,|\(")].reset_index(drop=True)
    df_chunk = df_chunk[df_chunk["gene_id"].str.contains(";|,|\(") == False].reset_index(drop=True)

    i = 0
    while i < len(df_issues):
        if i % 1000 == 0:
            print(f"{round(i / 1000)} of {round(len(df_issues) / 1000)}")
        try:
            target_value = str(df_issues.loc[i, "gene_id"])
            if ";" in target_value:
                target_value, new_target_value = df_issues.loc[i, "gene_id"].split(";", 1)
                df_to_concat = pd.concat([df_to_concat, df_issues.iloc[[i]]], ignore_index=True)
                df_to_concat.loc[len(df_to_concat) - 1, "gene_id"] = target_value
                df_issues.loc[i, "gene_id"] = new_target_value
            elif "," in target_value:
                target_value, new_target_value = df_issues.loc[i, "gene_id"].split(",", 1)
                df_to_concat = pd.concat([df_to_concat, df_issues.iloc[[i]]], ignore_index=True)
                df_to_concat.loc[len(df_to_concat) - 1, "gene_id"] = target_value
                df_issues.loc[i, "gene_id"] = new_target_value
            elif "(" in target_value:
                df_issues.loc[i, "gene_id"] = target_value.split("(")[0]
                i += 1
            else:
                df_issues.loc[i, "gene_id"] = target_value
                i += 1
        except TypeError:
            df_issues = df_issues.drop(i).reset_index(drop=True)
        except ValueError:
            df_issues = df_issues.drop(i).reset_index(drop=True)

    df_chunk = pd.concat([df_chunk, df_to_concat, df_issues])
    df_chunk = df_chunk[df_chunk["gene_id"] != "None"]
    df_chunk["gene_id"] = pd.to_numeric(df_chunk.gene_id, downcast="unsigned")
    print("chunk complete")

    df_chunk.to_csv(ROOT_DIR + f"pubtator_central/temp/gene_chunk_{process_id}.csv", index=False)


def process_diseases_helper(param_tuple):
    """This is the helper function for parallel processing of disease2pubtatorcentral. It takes two parameters, the
    list of mesh terms to process and the index of the process, packed as a tuple.
    :param param_tuple: tuple-packed list of genes and process index
    """
    mesh_term_chunk, process_id = param_tuple
    # this prevents memory overflow as later processes wait longer, so earlier ones can already filter for relevant
    time.sleep(process_id * 15)
    print("initializing process")
    df_chunk = pd.read_csv(ROOT_DIR + "pubtator_central/disease2pubtatorcentral.csv",
                           index_col=False,
                           engine="c")
    df_chunk = df_chunk.loc[df_chunk.mesh_term.isin(mesh_term_chunk)]
    df_chunk.reset_index(drop=True, inplace=True)

    df_to_concat = pd.DataFrame(columns=df_chunk.columns)
    df_issues = df_chunk[df_chunk["mesh_term"].str.contains(";|,|\(")].reset_index(drop=True)
    df_issues = pd.concat([df_issues, df_chunk[df_chunk["mesh_term"].str.contains("MESH|OMIM") == False]],
                          ignore_index=True).reset_index(drop=True)
    df_chunk = df_chunk[df_chunk["mesh_term"].str.contains(";|,|\(") == False]
    df_chunk = df_chunk[df_chunk["mesh_term"].str.contains("MESH|OMIM")].reset_index(drop=True)

    i = 0
    while i < len(df_issues):
        if i % 1000 == 0:
            print(f"{len(df_issues)} left")
        try:
            target_value = str(df_issues.loc[i, "mesh_term"])
            if ";" in target_value:
                target_value, new_target_value = df_issues.loc[i, "mesh_term"].split(";", 1)
                df_to_concat = pd.concat([df_to_concat, df_issues.iloc[[i]]], ignore_index=True)
                df_to_concat.loc[len(df_to_concat) - 1, "mesh_term"] = target_value
                df_issues.loc[i, "mesh_term"] = new_target_value
            elif "," in target_value:
                target_value, new_target_value = df_issues.loc[i, "mesh_term"].split(",", 1)
                df_to_concat = pd.concat([df_to_concat, df_issues.iloc[[i]]], ignore_index=True)
                df_to_concat.loc[len(df_to_concat) - 1, "mesh_term"] = target_value
                df_issues.loc[i, "mesh_term"] = new_target_value
            elif "(" in target_value:
                df_issues.loc[i, "mesh_term"] = target_value.split("(")[0]
                i += 1
            if "MESH" not in target_value and "OMIM" not in target_value:
                df_to_concat = df_to_concat.drop(len(df_to_concat) - 1).reset_index(drop=True)
            if "MESH" not in new_target_value and "OMIM" not in new_target_value:
                df_issues = df_issues.drop(i).reset_index(drop=True)
            else:
                i += 1

        except TypeError:
            df_issues = df_issues.drop(i).reset_index(drop=True)
        except ValueError:
            df_issues = df_issues.drop(i).reset_index(drop=True)

    df_chunk = pd.concat([df_chunk, df_to_concat, df_issues])

    df_chunk.to_csv(ROOT_DIR + f"pubtator_central/temp/disease_chunk_{process_id}.csv", index=False)



class PubtatorCentral:
    """This class handles PubtatorCentral data.
    """
    def __init__(self, n_cores, download=True, bootstraps=1000):
        self.bootstrap_cores = 40
        self.num_runs = bootstraps / self.bootstrap_cores
        self.n_cores = n_cores
        self.matrix = pd.DataFrame()
        # self.n_cores = 1
        self.n_matrix_column_chunks = 100
        self.cutoff = 0.001  # cutoff of max accepted empirical pvalue per term to be considered significant
        self.exponent = 0.5  # exponent for correction on total number of publications linked to a gene

        self.all_genes_df = pd.read_csv(ROOT_DIR + "hgnc_protein_coding.tsv", index_col=False,
                                        usecols=["entrez_id"], sep="\t", dtype="Int32")
        self.all_genes_df = self.all_genes_df.dropna(subset=["entrez_id"])
        self.all_genes = list(self.all_genes_df["entrez_id"].unique())
        self.ratio_pval_df = pd.DataFrame(columns=["pos_count", "neg_count", f"ratio_observed"])

        if download:
            # if os.path.isfile(ROOT_DIR + "pubtator_central/gene2pubtatorcentral"):
            #     os.rename(ROOT_DIR + "pubtator_central/gene2pubtatorcentral",
            #               ROOT_DIR + "pubtator_central/gene2pubtatorcentral_old")
            #     if os.path.isfile(ROOT_DIR + "pubtator_central/disease2pubtatorcentral"):
            #         os.rename(ROOT_DIR + "pubtator_central/disease2pubtatorcentral",
            #                   ROOT_DIR + "pubtator_central/disease2pubtatorcentral_old")
            # try:
            #     wget.download("ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/gene2pubtatorcentral.gz",
            #                   ROOT_DIR + "pubtator_central/gene2pubtatorcentral.gz")
            #     wget.download("ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/disease2pubtatorcentral.gz",
            #                   ROOT_DIR + "pubtator_central/disease2pubtatorcentral.gz")
            #     with gzip.open(ROOT_DIR + "pubtator_central/gene2pubtatorcentral.gz", "rb") as f_in:
            #         with open(ROOT_DIR + "pubtator_central/gene2pubtatorcentral", "wb") as f_out:
            #             shutil.copyfileobj(f_in, f_out)
            #     with gzip.open(ROOT_DIR + "pubtator_central/disease2pubtatorcentral.gz", "rb") as f_in:
            #         with open(ROOT_DIR + "pubtator_central/disease2pubtatorcentral", "wb") as f_out:
            #             shutil.copyfileobj(f_in, f_out)
            # except (ValueError, TypeError, AttributeError):
            #     print("Could not load Disgenet data! Using old data instead.")
            #     if os.path.isfile(ROOT_DIR + "pubtator_central/gene2pubtatorcentral_old"):
            #         os.rename(ROOT_DIR + "pubtator_central/gene2pubtatorcentral_old",
            #                   ROOT_DIR + "pubtator_central/gene2pubtatorcentral")
            #         if os.path.isfile(ROOT_DIR + "pubtator_central/disease2pubtatorcentral_old"):
            #             os.rename(ROOT_DIR + "pubtator_central/disease2pubtatorcentral_old",
            #                       ROOT_DIR + "pubtator_central/disease2pubtatorcentral")

#temrporary!
            self.gene_df = pd.read_csv(ROOT_DIR + "pubtator_central/gene_df_processed.csv")
            self.gene_disease_df = pd.read_csv(ROOT_DIR + "pubtator_central/gene_disease_df.csv",
                                               index_col=False,
                                               # nrows=10000000,
                                               dtype={"gene_id": "uint32", "pmid": "uint32", "mesh_term": "category"})

            self.num_splits = 100  # number of chunk to split the matrix into for parallel processing
            # self.preprocess_data()
            self.create_matrix()
        else:
            self.gene_df = pd.read_csv(ROOT_DIR + "pubtator_central/gene_df_processed.csv")
            self.gene_disease_df = pd.read_csv(ROOT_DIR + "pubtator_central/gene_disease_df.csv",
                                               index_col=False,
                                               # nrows=1000000,
                                               dtype={"gene_id": "uint32", "pmid": "uint32", "mesh_term": "category"})

        self.create_pmid_lists()
        # self.bootstrap()
        # self.calculate_mesh_ratios()
        # self.calculate_empirical_p()
        # self.gene_scores_exp()
        # self.lin_rank_genes()

    def preprocess_data(self):
        """This method needs two files: gene2pubtatorcentral and disease2pubtatorcentral. It first cleans the datasets
        and then merges them into the gene_disease_df.csv which is stored for later calculations.
        """
        # preprocessing_cores = 25  # the number of cores is being limited in order to prevent memory overflow
        preprocessing_cores = self.n_cores

        # cleaning "gene2pubtatorcentral" file
        print("loading gene_df...")
        gene_df = pd.read_csv(ROOT_DIR + "pubtator_central/gene2pubtatorcentral",
                              header=None,
                              engine="c",
                              sep="\t")
        gene_df.columns = ["pmid", "category", "gene_id", "gene_description", "placeholder"]

        gene_df = gene_df[["pmid", "gene_id"]]
        gene_df.to_csv(ROOT_DIR + "pubtator_central/gene2pubtatorcentral.csv", index=False)

        gene_list_chunks = [self.all_genes[round(len(self.all_genes) * i / preprocessing_cores):round(
            len(self.all_genes) * (i + 1) / preprocessing_cores)] for i in range(preprocessing_cores)]
        del gene_df
        print("chunked")

        with ProcessPoolExecutor() as executor:
            executor.map(process_gene_df_helper, zip(gene_list_chunks, list(range(preprocessing_cores))))

        gene_df_processed = pd.DataFrame()
        for chunk_name in os.listdir(ROOT_DIR + f"pubtator_central/temp/"):
            chunk = pd.read_csv(ROOT_DIR + f"pubtator_central/temp/" + chunk_name)
            gene_df_processed = pd.concat([gene_df_processed, chunk])
            os.remove(os.path.join(ROOT_DIR + f"pubtator_central/temp/", chunk_name))

        # print(gene_df_processed.head())
        print("filtering for entrez genes")
        gene_df_processed = gene_df_processed.loc[gene_df_processed.gene_id.isin(self.all_genes)].reset_index(drop=True)
        gene_df_processed.drop_duplicates(subset=["gene_id", "pmid"], inplace=True)
        gene_df_processed.to_csv(ROOT_DIR + "pubtator_central/gene_df_processed.csv",
                                 index=False)

        # cleaning "disease2pubtatorcentral" file
        print("loading disease_df...")
        disease_df = pd.read_csv(ROOT_DIR + "pubtator_central/disease2pubtatorcentral",
                                 header=None,
                                 engine="c",
                                 # nrows=100000,
                                 sep="\t")
        disease_df.columns = ["pmid", "category", "mesh_term", "disease_description", "placeholder"]

        disease_df = disease_df[["pmid", "mesh_term"]]
        disease_df.to_csv(ROOT_DIR + "pubtator_central/disease2pubtatorcentral.csv", index=False)

        # temp for debugging
        disease_df = pd.read_csv(ROOT_DIR + "pubtator_central/disease2pubtatorcentral.csv")

        all_mesh_terms = disease_df.mesh_term.unique()
        mesh_chunks = [list(all_mesh_terms)[round(len(all_mesh_terms) * i / preprocessing_cores):
                                            round(len(all_mesh_terms) * (i + 1) / preprocessing_cores)]
                       for i in range(preprocessing_cores)]
        del disease_df
        print("chunked")

        with ProcessPoolExecutor() as executor:
            executor.map(process_diseases_helper, zip(mesh_chunks, list(range(preprocessing_cores))))

        disease_df_processed = pd.DataFrame()
        for chunk_name in os.listdir(ROOT_DIR + f"pubtator_central/temp/"):
            chunk = pd.read_csv(ROOT_DIR + f"pubtator_central/temp/" + chunk_name)
            disease_df_processed = pd.concat([disease_df_processed, chunk])
            os.remove(os.path.join(ROOT_DIR + f"pubtator_central/temp/", chunk_name))

        print("concatenated")
        disease_df_processed.drop_duplicates(subset=["pmid", "mesh_term"], inplace=True)
        disease_df_processed.to_csv(ROOT_DIR + "pubtator_central/disease_df_processed.csv",
                                    index=False)

        self.gene_disease_df = disease_df_processed.merge(gene_df_processed, on="pmid")
        self.gene_disease_df = self.gene_disease_df.drop_duplicates()
        self.gene_disease_df.to_csv(ROOT_DIR + "pubtator_central/gene_disease_df.csv",
                                    index=False)

        print("preprocessing done")

    @ray.remote
    def group_df_helper(self, gene_list_chunk):
        temp = self.gene_disease_df[self.gene_disease_df.gene_id.isin(gene_list_chunk)]
        temp = temp.groupby(['gene_id', 'mesh_term']).size().reset_index(name="count")
        temp = temp.loc[temp["count"] != 0]
        # temp

    def create_pmid_lists(self):
        num_gene_splits = 250
        gene_list = list(self.gene_disease_df.gene_id.unique())
        df = self.gene_disease_df[['gene_id', 'mesh_term']]
        gene_list = list(self.gene_disease_df.gene_id.unique())
        gene_list_chunks = [gene_list[round(i * len(gene_list) / num_gene_splits):
                                      round((i + 1) * len(gene_list) / num_gene_splits)] for i in
                            range(num_gene_splits)]
        results_df = self.gene_disease_df.groupby(['gene_id', 'mesh_term']).size().reset_index(name="count")
        results_df = results_df.loc[results_df["count"] != 0]
        # [self.group_df_helper.remote(self, gene_lists_chunk) for gene_lists_chunk in gene_list_chunks]
        #
        # for df_chunk in ray.get(df_iterator):
        #     results_df = pd.concat([results_df, df_chunk])

        with open(ROOT_DIR + "pubtator_central/gene_mesh_count.pickle", "wb") as pickle_file:
            pickle.dump(results_df, pickle_file)

        print("done")


    def create_pmid_lists_delete(self):
        MATRIX_CORES = 15  # using all cores would exhaust memory

        # splitting dataframe, so the parts fit into the memory
        gene_list = list(self.gene_disease_df.gene_id.unique())
        gene_list_chunks = [gene_list[round(i * len(gene_list) / self.num_splits):
                                      round((i + 1) * len(gene_list) / self.num_splits)] for i in
                            range(self.num_splits)]

        # loads part by part and parallel processes each part on n_cores cores
        df_chunks = os.listdir(ROOT_DIR + "pubtator_central/matrix/temp/df_chunks/")

        for i, chunk_name in enumerate(df_chunks):
            print(f"chunk {i + 1} of {len(df_chunks)}")
            self.df_chunk = pd.read_csv(ROOT_DIR + f"pubtator_central/matrix/temp/df_chunks/{chunk_name}",
                                        index_col=False)
            gene_list = list(self.df_chunk.gene_id.unique())
            gene_list_chunks = [
                gene_list[round(i * len(gene_list) / MATRIX_CORES):round((i + 1) * len(gene_list) / MATRIX_CORES)] for i
                in range(MATRIX_CORES)]
            parameters_to_map = zip(gene_list_chunks, range(len(gene_list_chunks)),
                                    [chunk_name] * len(gene_list_chunks))

            print("initializing parallel processing")
            with ProcessPoolExecutor() as executor:
                executor.map(matrix_helper, parameters_to_map)

            print("concatenating chunks...")
            matrix_chunk = pd.DataFrame()
            for matrix_split_chunk_name in os.listdir(ROOT_DIR + "pubtator_central/matrix/temp/matrix_chunks/"):
                matrix_split_chunk = pd.read_csv(
                    ROOT_DIR + f"pubtator_central/matrix/temp/matrix_chunks/{matrix_split_chunk_name}")
                matrix_chunk = pd.concat([matrix_chunk, matrix_split_chunk])
            self.matrix = pd.concat([self.matrix, matrix_chunk])
            self.matrix.to_csv(ROOT_DIR + "pubtator_central/matrix/matrix.csv", index=False)
            del matrix_chunk
        print("matrix created...")
        # self.matrix.to_csv(ROOT_DIR + "pubtator_central/matrix/matrix.csv", index=False)

    # def create_matrix(self):
    #     """This method converts the raw data into a matrix, where rows are genes and columns are mesh_terms and the
    #     cells are lists of pmids containing both gene_id and mesh_term. This matrix is then used to do bootstrapping."""
    #     MATRIX_CORES = 15  # using all cores would exhaust memory
    #
    #     # splitting dataframe, so the parts fit into the memory
    #     gene_list = list(self.gene_disease_df.gene_id.unique())
    #     gene_list_chunks = [gene_list[round(i * len(gene_list) / self.num_splits):
    #                                   round((i + 1) * len(gene_list) / self.num_splits)] for i in
    #                         range(self.num_splits)]
    #
    #     print("creating temporary data splits")
    #     for i, gene_list_chunk in enumerate(gene_list_chunks):
    #         self.gene_disease_df.loc[self.gene_disease_df.gene_id.isin(gene_list_chunk)].to_csv(
    #             ROOT_DIR + f"pubtator_central/matrix/temp/df_chunks/df_chunk_{i}.csv", index=False)
    #         self.gene_disease_df = self.gene_disease_df.loc[self.gene_disease_df.gene_id.isin(gene_list_chunk) == False]
    #
    #     # loads part by part and parallel processes each part on n_cores cores
    #     df_chunks = os.listdir(ROOT_DIR + "pubtator_central/matrix/temp/df_chunks/")
    #     for i, chunk_name in enumerate(df_chunks):
    #         print(f"chunk {i + 1} of {len(df_chunks)}")
    #         self.df_chunk = pd.read_csv(ROOT_DIR + f"pubtator_central/matrix/temp/df_chunks/{chunk_name}",
    #                                     index_col=False)
    #         gene_list = list(self.df_chunk.gene_id.unique())
    #         gene_list_chunks = [
    #             gene_list[round(i * len(gene_list) / MATRIX_CORES):round((i + 1) * len(gene_list) / MATRIX_CORES)] for i
    #             in range(MATRIX_CORES)]
    #         parameters_to_map = zip(gene_list_chunks, range(len(gene_list_chunks)),
    #                                 [chunk_name] * len(gene_list_chunks))
    #
    #         print("initializing parallel processing")
    #         with ProcessPoolExecutor() as executor:
    #             executor.map(matrix_helper, parameters_to_map)
    #
    #         print("concatenating chunks...")
    #         matrix_chunk = pd.DataFrame()
    #         for matrix_split_chunk_name in os.listdir(ROOT_DIR + "pubtator_central/matrix/temp/matrix_chunks/"):
    #             matrix_split_chunk = pd.read_csv(
    #                 ROOT_DIR + f"pubtator_central/matrix/temp/matrix_chunks/{matrix_split_chunk_name}")
    #             matrix_chunk = pd.concat([matrix_chunk, matrix_split_chunk])
    #         self.matrix = pd.concat([self.matrix, matrix_chunk])
    #         self.matrix.to_csv(ROOT_DIR + "pubtator_central/matrix/matrix.csv", index=False)
    #         del matrix_chunk
    #     print("matrix created...")
    #     # self.matrix.to_csv(ROOT_DIR + "pubtator_central/matrix/matrix.csv", index=False)

    def bootstrap(self):
        self.matrix_columns = pd.read_csv(ROOT_DIR + "pubtator_central/matrix/matrix.csv", nrows=1).columns.to_list()
        len_column_chunk = len(self.matrix_columns[1:]) / self.n_matrix_column_chunks
        column_chunks = [[self.matrix_columns[0]] + self.matrix_columns[round(i * len_column_chunk) + 1:round(
            (i + 1) * len_column_chunk) + 1]
                         for i in range(self.n_column_chunks)]
        self.bootstrap_ratio_df = pd.DataFrame()

        # counter for indexing chunks
        i = 0
        for n_chunk, column_chunk in enumerate(column_chunks[0:]):
            # print(f"\n\n\nchunk {n_chunk + 1} of {len(column_chunks)}\n")
            self.column_chunk = column_chunk[1:]

            # creates empty dataframe to be filled later on
            ratio_df_chunk = pd.DataFrame.from_dict({"mesh_term": self.column_chunk})

            for run in range(self.num_runs):
                # each set of mesh terms is being bootstrapped n_cores * n_run times
                print(f"initializing process {i + 1} of {self.num_runs * len(column_chunks)}")
                i += 1
                process_numbers = list(
                    range((run * self.bootstrap_cores) + (n_chunk * self.bootstrap_cores * self.num_runs),
                          ((run + 1) * self.bootstrap_cores) + (n_chunk * self.bootstrap_cores * self.num_runs)))

                with ProcessPoolExecutor() as executor:
                    executor.map(self.bootstrap_helper, process_numbers)

            # print("concatenating chunks...")
            for bootstrap_process_number in range((n_chunk * self.bootstrap_cores * self.num_runs),
                                                  ((n_chunk + 1) * self.bootstrap_cores * self.num_runs)):
                bootstrap_df = pd.read_csv(
                    ROOT_DIR + f"pubtator_central/bootstrap/temp/bootstrap_{bootstrap_process_number}.csv")
                ratio_df_chunk = ratio_df_chunk.merge(bootstrap_df,
                                                      on="mesh_term")

            # ratio_df = pd.concat([ratio_df, ratio_df_chunk.iloc[:, 1:]], axis=1)

            ratio_df_chunk.to_csv(
                ROOT_DIR + f"pubtator_central/bootstrap/bootstrap_chunks/bootstrap_chunk_{i}.csv", index=False)
            i += 1

            filelist = os.listdir(ROOT_DIR + f"pubtator_central/bootstrap/temp/")
            for f in filelist:
                os.remove(os.path.join(ROOT_DIR + f"pubtator_central/bootstrap/temp/", f))

            # ratio_df.columns = ["mesh_term"] + [f"ratio_{n}" for n in range(self.bootstrap_cores * self.num_runs)]

        print("bootstrap successfull, concatenating straps...")
        chunk_list = os.listdir(ROOT_DIR + f"pubtator_central/bootstrap/bootstrap_chunks/")
        for chunk_name in chunk_list:
            chunk = pd.read_csv(os.path.join(ROOT_DIR + f"pubtator_central/bootstrap/bootstrap_chunks/", chunk_name))
            self.bootstrap_ratio_df = pd.concat([self.bootstrap_ratio_df, chunk])
        self.bootstrap_ratio_df.to_csv(ROOT_DIR + "pubtator_central/bootstrap/bootstrap_ratio_df.csv",
                                       index=False)
        print("bootstrap done!!")

    def bootstrap_helper(self, process_number):
        # print(f"starting process {process_number}")
        random.seed(process_number)
        matrix = pd.read_csv(ROOT_DIR + "pubtator_central/matrix/matrix.csv",
                             # nrows=100,
                             converters={column: converter for column in self.column_chunk},
                             usecols=["gene_id"] + self.column_chunk)

        pos_list = random.sample(self.all_genes, len(sysid_primary))
        neg_list = random.sample(list(set(self.all_genes) - set(pos_list)), len(princeton_negative))

        pos_pmid_total = self.gene_df[self.gene_df.gene_id.isin(pos_list)].pmid.nunique()
        neg_pmid_total = self.gene_df[self.gene_df.gene_id.isin(neg_list)].pmid.nunique()
        ratio_num = process_number % (self.bootstrap_cores * self.num_runs)

        boot_df = pd.DataFrame(columns=[f"ratio_{ratio_num}"])

        for mesh_term in self.column_chunk:
            # how many pmids contain the given mesh-term and at least one gene of our list?
            # pos_pmid = list(set(itertools.chain.from_iterable(matrix.loc[matrix.gene_id.isin(pos_list), mesh_term])))
            pos_pmid = list(itertools.chain.from_iterable(matrix.loc[matrix.gene_id.isin(pos_list), mesh_term]))
            if "" in pos_pmid:
                pos_pmid.remove("")
            pos_count = len(pos_pmid)

            # neg_pmid = list(set(itertools.chain.from_iterable(matrix.loc[matrix.gene_id.isin(neg_list), mesh_term])))
            neg_pmid = list(itertools.chain.from_iterable(matrix.loc[matrix.gene_id.isin(neg_list), mesh_term]))
            if "" in neg_pmid:
                neg_pmid.remove("")
            neg_count = len(neg_pmid)

            ratio = round((1.0 * pos_count / pos_pmid_total) / ((1.0 + neg_count) / neg_pmid_total), 4)

            # boot_df.loc[mesh_term, "pos_count"] = pos_count
            # boot_df.loc[mesh_term, "neg_count"] = neg_count
            boot_df.loc[mesh_term, f"ratio_{ratio_num}"] = ratio

        boot_df = boot_df.reset_index().rename(columns={"index": "mesh_term"})
        boot_df.to_csv(ROOT_DIR + f"pubtator_central/bootstrap/temp/bootstrap_{process_number}.csv", index=False)
        # return (boot_df, process_number)

    def calculate_mesh_ratios(self):
        print("loading data...")
        converter = lambda x: [int(element) for element in
                               x.strip("[]").replace("'", "").split(", ")] if "[" in x else []
        mesh_columns = pd.read_csv(ROOT_DIR + "pubtator_central/matrix/matrix.csv", nrows=1).columns.to_list()[1:]

        matrix = pd.read_csv(ROOT_DIR + "pubtator_central/matrix/matrix.csv",
                             converters={column: converter for column in mesh_columns})
        gene_df = pd.read_csv(ROOT_DIR + "pubtator_central/gene_df_processed.csv")

        primary_genes = sysid_primary.entrez_id.to_list()
        negative_genes = princeton_negative

        pos_pmid_total = gene_df[gene_df.gene_id.isin(primary_genes)].pmid.nunique()
        neg_pmid_total = gene_df[gene_df.gene_id.isin(negative_genes)].pmid.nunique()

        for n, mesh_term in enumerate(mesh_columns):
            print(f"column {n} of {len(mesh_columns)}")
            # how many unique pmids contain the given mesh-term and at least one gene of our list?
            pos_pmid = list(itertools.chain.from_iterable(matrix.loc[matrix.gene_id.isin(primary_genes), mesh_term]))
            if "" in pos_pmid:
                pos_pmid.remove("")
            pos_count = len(pos_pmid)

            neg_pmid = list(itertools.chain.from_iterable(matrix.loc[matrix.gene_id.isin(negative_genes), mesh_term]))
            if "" in neg_pmid:
                neg_pmid.remove("")
            neg_count = len(neg_pmid)

            ratio = round((1.0 * pos_count / pos_pmid_total) / ((1.0 + neg_count) / neg_pmid_total), 4)

            self.ratio_pval_df.loc[mesh_term, "pos_count"] = pos_count
            self.ratio_pval_df.loc[mesh_term, "neg_count"] = neg_count
            self.ratio_pval_df.loc[mesh_term, f"ratio_ratio_pval_df"] = ratio

        self.ratio_pval_df = self.ratio_pval_df.reset_index().rename(columns={"index": "mesh_term"})
        self.ratio_pval_df.to_csv(ROOT_DIR + f"pubtator_central/mesh_ratios_ratio_pval_df.csv", index=False)

    def calculate_empirical_p(self):
        self.ratio_pval_df = pd.read_csv(ROOT_DIR + "pubtator_central/mesh_ratios_observed.csv")
        self.bootstrap_ratio_df = pd.read_csv(ROOT_DIR + "pubtator_central/bootstrap/bootstrap_ratio_df.csv")
        bootstraps = [f"ratio_{i}" for i in range(len(self.bootstrap_ratio_df.columns) - 1)]
        for n, mesh_term in enumerate(self.bootstrap_ratio_df.mesh_term):
            # n += 1
            print(n)
            count_equal_or_greater = \
            self.bootstrap_ratio_df.loc[self.bootstrap_ratio_df.mesh_term == mesh_term][bootstraps].apply(
                lambda x: x >=
                          self.ratio_pval_df.loc[self.ratio_pval_df.mesh_term == mesh_term, "ratio_observed"].values[
                              0]).sum(axis=1)[n]
            self.ratio_pval_df.loc[
                self.ratio_pval_df.mesh_term == mesh_term, "p_empirical"] = 1. * count_equal_or_greater / (
                    len(self.bootstrap_ratio_df.columns) - 1)

        self.ratio_pval_df.to_csv(ROOT_DIR + f"pubtator_central/bootstrap/mesh_p_values.csv", index=False)

    def gene_scores_exp(self):
        # in order to check if ratio_pval_df exists, otherwise load from data
        try:
            _ = self.ratio_pval_df
        except AttributeError:
            self.ratio_pval_df = pd.read_csv(ROOT_DIR + f"pubtator_central/bootstrap/mesh_p_values.csv")
        cut_string = str(self.cutoff).replace(".", ",")
        exponent_string = str(self.exponent).replace(".", ",")

        self.ratio_pval_df = self.ratio_pval_df.loc[self.ratio_pval_df.p_empirical <= self.cutoff]

        gene_disease_df_filtered = self.gene_disease_df[
            self.gene_disease_df.mesh_term.isin(self.ratio_pval_df.mesh_term)]
        self.gene_scores_df = gene_disease_df_filtered.groupby("gene_id").size().reset_index(name="gene_score")

        gene_count_df = self.gene_disease_df.groupby("gene_id").size().reset_index(name="gene_count")
        self.gene_scores_df = self.gene_scores_df.merge(gene_count_df, on="gene_id")

        self.gene_scores_df["gene_score"] = self.gene_scores_df.gene_score / self.gene_scores_df.gene_count.apply(
            self.expo)

        self.gene_scores_df.to_csv(
            ROOT_DIR + f"pubtator_central/gene_scores/gene_scores_p_cutoff_{cut_string}_exp_{exponent_string}.csv",
            index=False)
        print("done calculating gene scores")

    def expo(self, x):
        return x ** self.exponent

    def lin_rank_genes(self, version="p_cutoff_0,001_exp_0,5"):
        self.gene_scores_df = pd.read_csv(ROOT_DIR + f"pubtator_central/gene_scores/gene_scores_{version}.csv",
                                          index_col=False,
                                          usecols=["gene_id", "gene_score"])
        self.gene_scores_df.columns = ["entrez_id", "pubtator_score"]
        self.gene_scores_df = add_categories(self.gene_scores_df, "entrez_id", "entrez")

        self.gene_scores_df = lin_rank(self.gene_scores_df, "pubtator_score")
        self.gene_scores_df.to_csv("/home/johann/AutoCaSc_core/data/pubtator_central/gene_scores.csv", index=False)


def update_data(n_cores=psutil.cpu_count()):
    """The function calling all the updating function. For those functions where the latest dataset download is not
    dependent on version numbers, up to date data is automatically downloaded. Manual downloads are necessary for:
        - PsyMukB
        - GTEx
        - StringDB
    :param n_cores:
    :return:
    """
    # HGNC()
    # MGI(n_cores, download=True)
    # StringDB(n_cores)
    # PsyMuKB()
    # GTEx()
    # Disgenet(n_cores, download=True)
    PubtatorCentral(n_cores=38, download=False)
    # fuse_data()


# fuse_data()

# update_gnomad_data()
update_data()
# fuse_data()
# for corr_exp in np.arange(0.3, 0.7, 0.05):
#     for rnk_exp in np.arange(1.5, 3, 0.5):
#         evaluate_parameters(MGI(10).calculate_gene_scores(corr_exp, rnk_exp), "entrez_id", "gene_score", (corr_exp, rnk_exp))
