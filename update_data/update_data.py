import pickle
import math
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


def HGNC():
    # r = requests.get(
    #     "https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_locus_group&col=gd_pub_eg_id&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit")
    # if r.ok:
    #     hgnc_table = pd.read_csv(StringIO(r.text),
    #                              sep="\t",
    #                              header=0,
    #                              names=["hgnc_id", "gene_symbol", "refseq_id", "entrez_id", "ensemble_id"])
    #     if os.path.isfile(ROOT_DIR + "hgnc_protein_coding.tsv"):
    #         os.rename(ROOT_DIR + "hgnc_protein_coding.tsv",
    #                   ROOT_DIR + "hgnc_protein_coding.tsv_old")
    hgnc_table = pd.read_csv(ROOT_DIR + "hgnc_protein_coding.tsv",
                      sep="\t")
    hgnc_table.columns = ["gene_symbol", "gene_type", "entrez_id", "ensemble_id"]
    hgnc_table = hgnc_table.loc[hgnc_table.gene_type == "protein-coding gene"][["gene_symbol",
                                                                                "entrez_id",
                                                                                "ensemble_id"]]
    hgnc_table.dropna(subset=["entrez_id"], inplace=True)
    hgnc_table.astype({"entrez_id":"Int32"}, errors="ignore")
    hgnc_table.to_csv(ROOT_DIR + "hgnc_protein_coding.tsv",
                      sep="\t",
                      index=False)
    # else:
    #     print("Could not load HGNC data! Using old data instead.")


def fuse_data(validation_run=False):
    all_genes = pd.read_csv(ROOT_DIR + "hgnc_protein_coding.tsv",
                            sep="\t",
                            usecols=["entrez_id", "ensemble_id"])

    gtex = pd.read_csv(ROOT_DIR + "gtex/gene_scores.csv", usecols=["ensemble_id", "rank_score"])
    gtex.columns = ["ensemble_id", "gtex_score"]

    denovo = pd.read_csv(ROOT_DIR + "psymukb/all_genes_denovo_NDD.csv",
                         usecols=["entrez_id", "rank_score"])
    denovo.columns = ["entrez_id", "denovo_rank_score"]

    disgenet = pd.read_csv(ROOT_DIR + "disgenet/disgenet_gene_scores.csv",
                           usecols=["entrez_id", "gene_score"])
    disgenet.columns = ["entrez_id", "disgenet_score"]

    # mgi = pd.read_csv(ROOT_DIR + "mgi/gene_annoations.csv")
    # mgi.columns = ["entrez_id", "mgi_phenotypes", "mgi_neuro_behavioral"]
    # mgi = mgi[["entrez_id", "mgi_neuro_behavioral"]]

    mgi = pd.read_csv(ROOT_DIR + "mgi/mgi_gene_scores.csv", usecols=["entrez_id", "rank_score"])
    mgi.columns = ["entrez_id", "mgi_score"]

    pubtator = pd.read_csv(ROOT_DIR + "pubtator_central/gene_scores.csv",
                           usecols=["entrez_id", "gene_score"])
    pubtator.columns = ["entrez_id", "pubtator_score"]

    string = pd.read_csv(ROOT_DIR + "string/gene_scores.csv", usecols=["gene_id", "gene_score"])
    string.columns = ["ensemble_id", "string_score"]

    gnomad = pd.read_csv(ROOT_DIR + "gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.csv", sep="\t")
    gnomad = gnomad[["gene_id", "pLI", "oe_lof",
                     "oe_lof_lower", "oe_lof_upper",
                     "oe_mis", "oe_mis_lower",
                     "oe_mis_upper", "mis_z", "obs_hom_lof"]]
    gnomad.columns = ["ensemble_id", "pLI", "oe_lof",
                      "oe_lof_lower", "oe_lof_upper",
                      "oe_mis", "oe_mis_lower",
                      "oe_mis_upper", "mis_z", "obs_hom_lof"]

    gevir = pd.read_csv(ROOT_DIR + "gevir/gene_ranking.csv",
                        usecols=["gene_id", "virlof_ar_enrichment"])
    gevir.columns = ["ensemble_id", "virlof_ar_enrichment"]


    all_data = all_genes.merge(gtex, on="ensemble_id", how="outer")
    all_data = all_data.merge(denovo, on="entrez_id", how="outer")
    all_data = all_data.merge(disgenet, on="entrez_id", how="outer")
    all_data = all_data.merge(mgi, on="entrez_id", how="outer")
    all_data = all_data.merge(pubtator, on="entrez_id", how="outer")
    all_data = all_data.merge(string, on="ensemble_id", how="outer")
    all_data = all_data.merge(gevir, on="ensemble_id", how="outer")

    all_data = all_data.fillna(0)

    # careful with filling nans, o/e = 0 would bias results
    all_data = all_data.merge(gnomad, on="ensemble_id", how="outer")



    all_data = add_categories(all_data, "entrez_id", "entrez")



    # all_data = all_data.dropna(subset=["mgi_score", "string_score", "pubtator_score", "disgenet_score", "gtex_score", "denovo_rank_score"])

    all_data[["mgi_score", "string_score", "pubtator_score", "disgenet_score", "gtex_score", "denovo_rank_score"]] = \
        all_data[
            ["mgi_score", "string_score", "pubtator_score", "disgenet_score", "gtex_score",
             "denovo_rank_score"]].fillna(0)

    # all_data[
    #     "gene_sum"] = all_data.gtex_score + all_data.denovo_rank_score + all_data.disgenet_score + all_data.mgi_score + all_data.pubtator_score + all_data.string_score
    all_data["weighted_score"] = sum(
        [spearmanr(all_data[parameter], all_data.sys_primary)[0] * all_data[parameter] for parameter in
         ["mgi_score", "string_score", "pubtator_score", "disgenet_score", "gtex_score", "denovo_rank_score"]])
    for parameter in ["mgi_score", "string_score", "pubtator_score", "disgenet_score", "gtex_score", "denovo_rank_score"]:
        print(parameter, round(6.0 * spearmanr(all_data[parameter], all_data.sys_primary)[0] / all_data.weighted_score.max(), 3))

    print(all_data.weighted_score.max())
    all_data.weighted_score = all_data.weighted_score / all_data.weighted_score.max()

    all_data = all_data.sort_values(by="weighted_score", ascending=False).drop_duplicates(
        subset=["entrez_id", "ensemble_id"], keep="first")
    if validation_run:
        all_data.to_csv(ROOT_DIR + "all_gene_data_seed_42.csv", index=False)
    else:
        all_data.to_csv(ROOT_DIR + "all_gene_data.csv", index=False)

def update_gnomad_data():
    all_genes = pd.read_csv(ROOT_DIR + "protein_gene_table.csv")
    gnomad_data = pd.read_csv(ROOT_DIR + "gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.csv",
                              sep="\t",
                              usecols=["transcript", "pLI", "oe_lof", "oe_lof_lower", "oe_lof_upper",
                                       "oe_mis", "oe_mis_lower", "oe_mis_upper", "mis_z"])

    df = all_genes.merge(gnomad_data, left_on="transcript_id", right_on="transcript", how="left").sort_values(by="pLI",
                                                                                                              ascending=False)
    df = df.drop_duplicates(subset=["entrez_id"], keep="first")
    df = df.loc[~df.entrez_id.isnull()]

    df.to_csv(ROOT_DIR + "gnomad/gnomad_data.csv", index=False)

def clean_mgi_phenotype_df(df):
    df.mgi_id = df.mgi_id.str.strip()
    i = 0
    while i < len(df):
        print(f"{i + 1} of {len(df)}")
        if "," in df.loc[i, "mgi_id"]:
            first_id, rest = df.loc[i, "mgi_id"].split(",", 1)
            df.loc[i, "mgi_id"] = first_id
            df.loc[len(df), "phenotype"] = df.loc[i, "phenotype"]
            df.loc[len(df) - 1, "mgi_id"] = rest
        i += 1
    return df

class MGI:
    """
    This class takes in two files from MGI, containing a translation from entrez to MGI identifiers and all
    phenotypes attributed to a given genotype.
    """

    def __init__(self, n_cores=4, download=False, exponent=0.5, investigate_parameters=False):
        self.exponent = exponent
        self.investigate_parameters = investigate_parameters
        if download:
            if os.path.isfile(ROOT_DIR + "mgi/HMD_HumanPhenotype.rpt"):
                os.rename(ROOT_DIR + "mgi/HMD_HumanPhenotype.rpt",
                          ROOT_DIR + "mgi/HMD_HumanPhenotype_old.rpt")
            wget.download("http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt",
                          ROOT_DIR + "mgi/HMD_HumanPhenotype.rpt")
        self.entrez_mgi = pd.read_csv(ROOT_DIR + "mgi/HMD_HumanPhenotype.rpt", sep="\t",
                                      usecols=[0, 1, 5])
        self.entrez_mgi.columns = ["gene_symbol", "entrez_id", "mgi_id"]
        self.entrez_mgi.mgi_id = self.entrez_mgi.mgi_id.str.strip()
        self.entrez_mgi.dropna(subset=["entrez_id"], inplace=True)
        self.entrez_mgi.astype({'entrez_id': 'int32'})

        self.mgi_mp = pd.read_csv(ROOT_DIR + "mgi/MGI_PhenoGenoMP.rpt.txt", sep="\t",
                                  usecols=[3, 5])
        self.mgi_mp.columns = ["phenotype", "mgi_id"]
        self.mgi_mp.mgi_id = self.mgi_mp.mgi_id.str.strip()
        # self.mgi_mp = pd.read_csv(ROOT_DIR + "mgi/MGI_PhenoGenoMP.rpt.cleaned", sep="\t")
        self.mgi_mp = clean_mgi_phenotype_df(self.mgi_mp)
        self.mgi_mp.to_csv(ROOT_DIR + "mgi/MGI_PhenoGenoMP.rpt.cleaned", sep="\t", index=False)
        # self.mgi_mp = pd.read_csv(ROOT_DIR + "mgi/MGI_PhenoGenoMP.rpt.cleaned", sep="\t")
        self.mgi_mp.dropna(inplace=True)
        self.df = self.entrez_mgi.merge(self.mgi_mp, on="mgi_id")
        self.df.drop_duplicates(inplace=True)
        self.df.loc[:, "entrez_id"] = pd.to_numeric(self.df.entrez_id, downcast="signed")
        self.n_cores = n_cores
        self.pval_df = None

    def update(self):
        self.generate_empricial_p_values()
        self.calculate_gene_scores()

    def generate_empricial_p_values(self):
        self.observed_df = self.df.loc[self.df.entrez_id.isin(sysid_primary.entrez_id.to_list())].groupby(
            "phenotype").size().reset_index(
            name="count").sort_values(by="count", ascending=False)

        # term_count_pos.to_csv(ROOT_DIR + "mgi/mp_count_observed.csv", index=False)

        """
        This part performs bootstraps on the data to identify phenotypes specific to NDD-genes.
        """

        # sufficient for bonferroni correction:
        n_experiments = 20 * len(self.observed_df) + 1
        # n_experiments = 100

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

    def calculate_gene_scores(self, rank_exp=2):
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
        gene_df["gene_score"] = gene_df.relevant_count / (gene_df.all_count ** self.exponent)
        gene_df = gene_df[["entrez_id", "gene_score"]]

        """
        Finally gene scores are calculated by ranking genes and taking their rank score to the power of two.
        """
        self.gene_scores = rank_genes(gene_df, "gene_score", rank_exp)

        if not self.investigate_parameters:
            self.gene_scores.to_csv(ROOT_DIR + "mgi/mgi_gene_scores.csv", index=False)
        else:
            self.gene_scores.to_csv(ROOT_DIR + f"mgi/mgi_gene_scores_{str(self.exponent).replace('.', ',')}.csv", index=False)

        return self.gene_scores

    def mgi_bootstrap_helper(self, process_ids):
        df_mgi_subset = mgi_df.loc[mgi_df.phenotype.isin(self.observed_df.phenotype)][["entrez_id", "phenotype"]]
        chunk = pd.DataFrame(self.observed_df, columns=["phenotype"])
        for n, process in enumerate(process_ids):
            print(f"{n + 1} of {len(process_ids)}")
            random.seed(process)
            # pos_list = random.sample(set(self.df.entrez_id.unique()) - set(
            #     sysid_candidates.entrez_id.to_list() + sysid_primary.entrez_id.to_list()),
            #                          len(sysid_primary.entrez_id))
            pos_list = random.sample(set(self.df.entrez_id.unique()) - set(
                sysid_candidates.entrez_id.to_list() + sysid_primary.entrez_id.to_list()),
                                     self.df.loc[self.df.entrez_id.isin(sysid_primary.entrez_id)].entrez_id.nunique())
            # print(pos_list)

            df_pos = df_mgi_subset.loc[df_mgi_subset.entrez_id.isin(pos_list)]
            term_count_pos = df_pos.groupby("phenotype").size().reset_index(name=f"count_{process + 1}")
            # print(term_count_pos.loc[term_count_pos.phenotype == "MP:0002083"])

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
                    20 * len(self.observed_df))

        return chunk[["phenotype", "p_empirical"]]


def evaluate_parameters(df, nomenclature, parameter, version):
    if nomenclature == "entrez_id":
        df["sys_primary"] = df.entrez_id.isin(sysid_primary.entrez_id).astype(int)
        df["sys_candidate"] = df.entrez_id.isin(sysid_candidates.entrez_id).astype(int)
        df["sys"] = df.entrez_id.isin(sysid_primary.entrez_id + sysid_candidates.entrez_id).astype(
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
        #         if os.path.isfile(ROOT_DIR + "string/9606.protein.links.detailed.v11.0.txt"):
        #             os.rename(ROOT_DIR + "string/9606.protein.links.detailed.v11.0.txt",
        #                       ROOT_DIR + "string/9606.protein.links.detailed.v11.0.txt_old")
        #         with open(ROOT_DIR + "string/9606.protein.links.detailed.v11.0.txt", 'wb') as f:
        #             r.raw.decode_content = True  # just in case transport encoding was applied
        #             gzip_file = gzip.GzipFile(fileobj=r.raw)
        #             shutil.copyfileobj(gzip_file, f)
        #     else:
        #         print("Could not load StrinDB data! Using old data instead.")
        self.string_data = pd.read_csv(ROOT_DIR + "string/9606.protein.links.detailed.v11.0.txt",
                                       sep=" ",
                                       usecols=["protein1", "protein2", "combined_score"])

        r = requests.get(
            'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_peptide_id" /></Dataset></Query>')
        if r.ok:
            self.pg_table = pd.read_csv(StringIO(r.text), sep="\t", names=["ensemble_id", "protein_id"])
        else:
            print("Could not load biomart information from Ensebml BioMart! Using old data instead.")
            self.pg_table = pd.read_csv(ROOT_DIR + "protein_gene_table.csv",
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
        gene_scores_df.to_csv(ROOT_DIR + "string/gene_scores.csv", index=False)

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

        all_genes = pd.read_csv(ROOT_DIR + "hgnc_protein_coding.tsv",
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
        denovo_count.to_csv(ROOT_DIR + "psymukb/all_genes_denovo_NDD.csv", index=False)


class GTEx:
    def __init__(self):
        gtex_data = pd.read_csv(ROOT_DIR + "gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                                sep="\t",
                                skiprows=2)

        gtex_data["ensemble_id"] = gtex_data.Name.apply(lambda x: x.split(".")[0])
        gtex_data = gtex_data.drop(columns=["Name"])

        translation_df = pd.read_csv(ROOT_DIR + "hgnc_protein_coding.tsv", sep="\t")
        gtex_data = gtex_data.merge(translation_df[["entrez_id", "ensemble_id"]], on="ensemble_id")

        columns = gtex_data.columns
        brain_columns = [x for x in columns if "Brain" in x]

        gtex_data["brain_sum"] = gtex_data[brain_columns].sum(axis=1)

        gtex_data = gtex_data[["ensemble_id", "entrez_id", "brain_sum"]]

        gtex_data = rank_genes(gtex_data, "brain_sum", 4)
        gtex_data.to_csv(ROOT_DIR + "gtex/gene_scores.csv", index=False)


class Disgenet:
    def __init__(self, n_cores=4, download=False):
        if download:
            r = requests.get(
                "https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz",
                stream=True)
            if r.status_code == 200:
                if os.path.isfile(ROOT_DIR + "disgenet/all_gene_disease_associations.tsv"):
                    os.rename(ROOT_DIR + "disgenet/all_gene_disease_associations.tsv",
                              ROOT_DIR + "disgenet/all_gene_disease_associations.tsv_old")
                with open(ROOT_DIR + "disgenet/all_gene_disease_associations.tsv", 'wb') as f:
                    r.raw.decode_content = True  # just in case transport encoding was applied
                    gzip_file = gzip.GzipFile(fileobj=r.raw)
                    shutil.copyfileobj(gzip_file, f)
            else:
                print("Could not load Disgenet data! Using old data instead.")

        self.n_cores = n_cores
        self.df = pd.read_csv(ROOT_DIR + "disgenet/all_gene_disease_associations.tsv", sep="\t",
                              usecols=["geneId", "diseaseName", "diseaseClass", "score"])
        self.df = self.df.loc[self.df.diseaseClass.str.contains("C10|F03", na=False)]
        self.all_genes = list(self.df.geneId.unique())

        self.observe_term_ratio()
        self.bootstrap()
        self.score_genes()

    def observe_term_ratio(self):
        df_pos = self.df.loc[self.df.geneId.isin(sysid_primary.entrez_id)]
        observed_count_pos = df_pos.groupby("diseaseName").size().reset_index(name="count").sort_values(by="count",
                                                                                                        ascending=False)

        df_neg = self.df.loc[self.df.geneId.isin(negative_gene_list)]
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
            neg_list = random.sample(set(self.all_genes) - set(pos_list), len(negative_gene_list))

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
                    20 * len(self.observed_df))

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


@ray.remote
def bootstrap_helper(process_id_list, df, process_id_total_counts_dict, process_id_series_dict):
    ratio_list = []
    for i, process_id in enumerate(process_id_list):
        # ToDo sollte man auch Gene mit einbeziehen, über de gar nix publiziert wurde?
        pos_series, neg_series = process_id_series_dict.get(process_id)
        pos_pmid_total, neg_pmid_total = process_id_total_counts_dict.get(process_id)

        ratio = calculate_mesh_ratios(pos_series, neg_series, df, pos_pmid_total, neg_pmid_total)

        ratio_list.append(ratio)

    return ratio_list


def calculate_mesh_ratios(pos_series, neg_series, df, pos_pmid_total, neg_pmid_total):
    pos_pmid_mesh_count = df.join(pos_series, how="inner")["count"].sum()
    neg_pmid_mesh_count = df.join(neg_series, how="inner")["count"].sum()

    ratio = round((1.0 * pos_pmid_mesh_count / pos_pmid_total) /
                  ((1.0 + neg_pmid_mesh_count) / neg_pmid_total), 4)
    return ratio

def create_pmid_lists_helper(df_chunk):
    results_chunk = df_chunk.groupby(['gene_id', 'mesh_term']).size().reset_index(name="count")
    results_chunk = results_chunk.loc[results_chunk["count"] != 0]
    return results_chunk


@ray.remote
def create_process_id_dicts(process_ids, gene_df):
    total_counts_dict_chunk = {}
    series_dict_chunk = {}
    for i, process_id in enumerate(process_ids):
        print(f"creating dict part {i + 1} of {len(process_ids)}")
        random.seed(process_id)
        pos_list = random.sample(list(gene_df.gene_id.unique()), len(sysid_primary))
        neg_list = random.sample(list(set(list(gene_df.gene_id.unique())) - set(pos_list)), len(negative_gene_list))
        pos_pmid_total = gene_df[gene_df.gene_id.isin(pos_list)].pmid.nunique()
        neg_pmid_total = gene_df[gene_df.gene_id.isin(neg_list)].pmid.nunique()
        total_counts_dict_chunk[process_id] = (pos_pmid_total, neg_pmid_total)
        series_dict_chunk[process_id] = (pd.Series(pos_list, index=pos_list, name=f"process_{process_id}_pos"),
                                         pd.Series(neg_list, index=neg_list, name=f"process_{process_id}_neg"))
    return total_counts_dict_chunk, series_dict_chunk


class PubtatorCentral:
    """This class handles PubtatorCentral data.
    """

    def __init__(self, n_cores, download=False, preprocess=False, n_bootstraps=1001, exponent=0.7, cutoff=0.0001):
        # self.bootstrap_cores = 40
        self.n_bootstraps = n_bootstraps
        # self.num_runs = bootstraps / self.bootstrap_cores
        self.n_cores = n_cores
        self.matrix = pd.DataFrame()
        self.num_splits = 100  # number of chunk to split the matrix into for parallel processing
        # self.n_cores = 1
        self.n_matrix_column_chunks = 100
        self.cutoff = cutoff  # cutoff of max accepted empirical pvalue per term to be considered significant
        self.exponent = exponent  # exponent for correction on total number of publications linked to a gene
        self.version = f'p_cutoff_{str(self.cutoff).replace(".", ",")}_exp_{str(self.exponent).replace(".", ",")}'


        self.all_genes_df = pd.read_csv(ROOT_DIR + "hgnc_protein_coding.tsv", index_col=False,
                                        usecols=["entrez_id"], sep="\t", dtype="Int32")
        self.all_genes_df = self.all_genes_df.dropna(subset=["entrez_id"])
        self.all_genes = list(self.all_genes_df["entrez_id"].unique())

        if download:
            if os.path.isfile(ROOT_DIR + "pubtator_central/gene2pubtatorcentral"):
                os.rename(ROOT_DIR + "pubtator_central/gene2pubtatorcentral",
                          ROOT_DIR + "pubtator_central/gene2pubtatorcentral_old")
                if os.path.isfile(ROOT_DIR + "pubtator_central/disease2pubtatorcentral"):
                    os.rename(ROOT_DIR + "pubtator_central/disease2pubtatorcentral",
                              ROOT_DIR + "pubtator_central/disease2pubtatorcentral_old")
            try:
                wget.download("ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/gene2pubtatorcentral.gz",
                              ROOT_DIR + "pubtator_central/gene2pubtatorcentral.gz")
                wget.download("ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/disease2pubtatorcentral.gz",
                              ROOT_DIR + "pubtator_central/disease2pubtatorcentral.gz")
                with gzip.open(ROOT_DIR + "pubtator_central/gene2pubtatorcentral.gz", "rb") as f_in:
                    with open(ROOT_DIR + "pubtator_central/gene2pubtatorcentral", "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                with gzip.open(ROOT_DIR + "pubtator_central/disease2pubtatorcentral.gz", "rb") as f_in:
                    with open(ROOT_DIR + "pubtator_central/disease2pubtatorcentral", "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
            except (ValueError, TypeError, AttributeError):
                print("Could not load Disgenet data! Using old data instead.")
                if os.path.isfile(ROOT_DIR + "pubtator_central/gene2pubtatorcentral_old"):
                    os.rename(ROOT_DIR + "pubtator_central/gene2pubtatorcentral_old",
                              ROOT_DIR + "pubtator_central/gene2pubtatorcentral")
                    if os.path.isfile(ROOT_DIR + "pubtator_central/disease2pubtatorcentral_old"):
                        os.rename(ROOT_DIR + "pubtator_central/disease2pubtatorcentral_old",
                                  ROOT_DIR + "pubtator_central/disease2pubtatorcentral")

        if preprocess:
            self.preprocess_data()
            self.create_pmid_lists()

        else:
            self.gene_df = pd.read_csv(ROOT_DIR + "pubtator_central/gene_df_processed.csv")
            self.gene_df.loc[:, "gene_id"] = pd.to_numeric(self.gene_df.loc[:, "gene_id"],
                                                           downcast="unsigned")
            self.gene_df.loc[:, "pmid"] = pd.to_numeric(self.gene_df.loc[:, "pmid"],
                                                        downcast="unsigned")
            with open(ROOT_DIR + "pubtator_central/gene_disease_df.pickle", "rb") as pickle_file:
                self.gene_disease_df = pickle.load(pickle_file)

        # self.bootstrap()
        self.gene_scores_exp()
        self.lin_rank_genes(investigate_parameters=True)

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
        with open(ROOT_DIR + "pubtator_central/gene_disease_df.pickle", "wb") as pickle_file:
            pickle.dump(self.gene_disease_df, pickle_file)

        print("preprocessing done")

    def create_pmid_lists(self):
        all_genes = list(self.gene_disease_df.gene_id.unique())
        # all_genes = all_genes[:400]
        n_splits = 6
        gene_chunks = [all_genes[round(i * len(all_genes) / n_splits):round((i + 1) * len(all_genes) / n_splits)]
                       for i in range(n_splits)]
        df_chunks = [self.gene_disease_df.loc[self.gene_disease_df.gene_id.isin(gene_chunk)].reset_index(drop=True)
                     for gene_chunk in gene_chunks]
        results_df = pd.DataFrame()
        with ProcessPoolExecutor() as executor:
            results = executor.map(create_pmid_lists_helper, df_chunks)
            for _ in results:
                results_df = pd.concat([results_df, _], ignore_index=True)

        with open(ROOT_DIR + "pubtator_central/gene_mesh_count.pickle", "wb") as pickle_file:
            pickle.dump(results_df, pickle_file)

    def bootstrap(self):
        create_new_total_count_dict = True
        ray.init(num_cpus=self.n_cores)
        with open(ROOT_DIR + "pubtator_central/gene_mesh_count.pickle", "rb") as pickle_file:
            self.gene_mesh_df = pickle.load(pickle_file)
            self.gene_mesh_df.loc[:, "gene_id"] = pd.to_numeric(self.gene_mesh_df.loc[:, "gene_id"],
                                                                downcast="unsigned")
            self.gene_mesh_df.loc[:, "count"] = pd.to_numeric(self.gene_mesh_df.loc[:, "count"],
                                                              downcast="unsigned")
        gene_df = ray.put(self.gene_df)
        len_process_chunks = self.n_bootstraps / self.n_cores
        process_id_lists = [list(range(round(i * len_process_chunks),
                                       round((i + 1) * len_process_chunks))) for i in range(self.n_cores)]

        if create_new_total_count_dict:
            create_dicts_temp = ray.get([create_process_id_dicts.remote(process_ids, gene_df)
                                              for process_ids in process_id_lists])
            process_id_total_counts_dict = {}
            process_id_series_dict = {}
            [process_id_total_counts_dict.update(total_chunk) for total_chunk, _ in create_dicts_temp]
            [process_id_series_dict.update(series_chunk) for _, series_chunk in create_dicts_temp]
            with open(ROOT_DIR + "pubtator_central/process_id_total_counts_dict.pickle", "wb") as pickle_file:
                pickle.dump(process_id_total_counts_dict, pickle_file)
            with open(ROOT_DIR + "pubtator_central/process_id_series_dict.pickle", "wb") as pickle_file:
                pickle.dump(process_id_series_dict, pickle_file)
        else:
            with open(ROOT_DIR + "pubtator_central/process_id_total_counts_dict.pickle", "rb") as pickle_file:
                process_id_total_counts_dict = pickle.load(pickle_file)
            with open(ROOT_DIR + "pubtator_central/process_id_series_dict.pickle", "rb") as pickle_file:
                process_id_series_dict = pickle.load(pickle_file)
        process_id_total_counts_dict_id = ray.put(process_id_total_counts_dict)
        process_id_series_dict_id = ray.put(process_id_series_dict)

        mesh_pval_dict = {}
        mesh_terms_to_check = list(self.gene_mesh_df.loc[
                                       self.gene_mesh_df.gene_id.isin(
                                           list(set(negative_gene_list + sysid_primary.entrez_id.to_list())))].mesh_term.unique())

        pos_pmid_total_original = self.gene_df[self.gene_df.gene_id.isin(sysid_primary.entrez_id)].pmid.nunique()
        neg_pmid_total_original = self.gene_df[self.gene_df.gene_id.isin(negative_gene_list)].pmid.nunique()

        pos_series_original = pd.Series(sysid_primary.entrez_id.to_list(), index=sysid_primary.entrez_id.to_list(),
                                        name="pos_series_original")
        neg_series_original = pd.Series(negative_gene_list, index=negative_gene_list,
                                        name="neg_series_original")

        for i, _mesh_term in enumerate(mesh_terms_to_check):
            print(f"Processing {i + 1} of {len(mesh_terms_to_check)}...")
            _mesh_df = self.gene_mesh_df.loc[self.gene_mesh_df.mesh_term == _mesh_term].set_index("gene_id")

            ratio_original = calculate_mesh_ratios(pos_series_original,
                                                   neg_series_original,
                                                   _mesh_df,
                                                   pos_pmid_total_original,
                                                   neg_pmid_total_original)

            _mesh_df_id = ray.put(_mesh_df)

            ray_returns = ray.get([bootstrap_helper.remote(process_ids,
                                                           _mesh_df_id,
                                                           process_id_total_counts_dict_id,
                                                           process_id_series_dict_id)
                                   for process_ids in process_id_lists])
            ratio_list = itertools.chain.from_iterable(ray_returns)
            ratio_list = [_ for _ in ratio_list if _ >= ratio_original]
            mesh_pval_dict[_mesh_term] = 1.0 * len(ratio_list) / self.n_bootstraps
            del ratio_list, ray_returns

        self.ratio_pval_df = pd.DataFrame.from_dict(mesh_pval_dict, orient='index')
        self.ratio_pval_df.reset_index(inplace=True)
        self.ratio_pval_df.columns = ["mesh_term", "p_empirical"]
        self.ratio_pval_df.to_csv(ROOT_DIR + f"pubtator_central/bootstrap/mesh_p_values.csv", index=False)
        print("bootstrap done!!")

    def gene_scores_exp(self):
        # in order to check if ratio_pval_df exists, otherwise load from data
        try:
            _ = self.ratio_pval_df
        except AttributeError:
            self.ratio_pval_df = pd.read_csv(ROOT_DIR + f"pubtator_central/bootstrap/mesh_p_values.csv")

        self.ratio_pval_df = self.ratio_pval_df.loc[self.ratio_pval_df.p_empirical <= self.cutoff]

        gene_disease_df_filtered = self.gene_disease_df[
            self.gene_disease_df.mesh_term.isin(self.ratio_pval_df.mesh_term)]
        self.gene_scores_df = gene_disease_df_filtered.groupby("gene_id").size().reset_index(name="gene_score")

        gene_count_df = self.gene_disease_df.groupby("gene_id").size().reset_index(name="gene_count")
        self.gene_scores_df = self.gene_scores_df.merge(gene_count_df, on="gene_id")

        self.gene_scores_df["gene_score"] = self.gene_scores_df.gene_score / self.gene_scores_df.gene_count.apply(
            self.expo)

        self.gene_scores_df.to_csv(
            ROOT_DIR + f"pubtator_central/gene_scores/gene_scores_{self.version}.csv",
            index=False)
        print("done calculating gene scores")

    def expo(self, x):
        return x ** self.exponent

    def lin_rank_genes(self, investigate_parameters=False):
        self.gene_scores_df = pd.read_csv(ROOT_DIR + f"pubtator_central/gene_scores/gene_scores_{self.version}.csv",
                                          index_col=False,
                                          usecols=["gene_id", "gene_score"])
        self.gene_scores_df.columns = ["entrez_id", "pubtator_score"]
        self.gene_scores_df = add_categories(self.gene_scores_df, "entrez_id", "entrez")

        self.gene_scores_df = lin_rank(self.gene_scores_df, "pubtator_score")
        if not investigate_parameters:
            self.gene_scores_df.to_csv(ROOT_DIR + "pubtator_central/gene_scores.csv", index=False)
        else:
            self.gene_scores_df.to_csv(ROOT_DIR + f"pubtator_central/gene_scores_{self.version}.csv", index=False)

def evaluate_pubtator_parameters():
    for exponent in [0.5, 0.6, 0.7, 0.8, 0.9]:
        for cutoff in [0.001, 0.0001, 0.00001]:
            instance = PubtatorCentral(n_cores=46, download=False, preprocess=False, exponent=exponent, cutoff=cutoff)
            df = pd.read_csv(ROOT_DIR + f"pubtator_central/gene_scores_{instance.version}.csv")
            p_val = mannwhitneyu(df.loc[df.sys_candidate == 1]['gene_score'].to_list(),
                                 df.loc[df.sys == 0]['gene_score'].to_list(), alternative='greater')[1]
            print(f"{instance.version}: {p_val}")
    """
    p_cutoff_0,001_exp_0,5: 1.0050455371702318e-120
    p_cutoff_0,0001_exp_0,5: 1.1206678948776022e-134
    p_cutoff_1e-05_exp_0,5: 2.536029138021977e-132
    p_cutoff_0,001_exp_0,6: 1.2182866524894282e-138
    p_cutoff_0,0001_exp_0,6: 8.469634478024799e-148
    p_cutoff_1e-05_exp_0,6: 3.912806508743547e-146
    p_cutoff_0,001_exp_0,7: 1.2431371571945905e-152
    p_cutoff_0,0001_exp_0,7: 3.2414372875537876e-157
    p_cutoff_1e-05_exp_0,7: 1.1597853586730318e-152
    p_cutoff_0,001_exp_0,8: 1.8797958968014228e-147
    p_cutoff_0,0001_exp_0,8: 3.1111897268660445e-140
    p_cutoff_1e-05_exp_0,8: 2.0275673897096633e-131
    p_cutoff_0,001_exp_0,9: 7.519095823436918e-97
    p_cutoff_0,0001_exp_0,9: 4.803675773940467e-95
    """

def evaluate_mgi_parameters():
    for exponent in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        MGI(n_cores=psutil.cpu_count(),
            download=False,
            exponent=exponent,
            investigate_parameters=True).calculate_gene_scores()
        df = pd.read_csv(ROOT_DIR + f"mgi/mgi_gene_scores_{str(exponent).replace('.', ',')}.csv")
        df = add_categories(df, "entrez_id", "entrez")
        p_val = mannwhitneyu(df.loc[df.sys_candidate == 1]['gene_score'].to_list(),
                             df.loc[df.sys == 0]['gene_score'].to_list(), alternative='greater')[1]
        print(f"exponent {exponent}: {p_val}")
    """
    exponent 0.2: 1.3413141828335412e-23
    exponent 0.3: 3.83851534927423e-25
    exponent 0.4: 1.8970705866357454e-26
    exponent 0.5: 1.6443960386822857e-27 --> using this one
    exponent 0.6: 9.73622834348643e-27
    exponent 0.7: 7.48339207409757e-24
    exponent 0.8: 2.6248413932730707e-19
    exponent 0.9: 1.0301802549921885e-14
    """

if __name__ == "__main__":
    ROOT_DIR = "/mnt/raid/users/johann/AutoCaSc_maintenance_data/"

    validation_run = False
    # HGNC()
    if validation_run:
        random.seed(42)

    sysid_primary = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv",
                                usecols=["Entrez id", "Ensembl id"])
    sysid_primary.columns = ["entrez_id", "ensemble_id"]
    sysid_candidates = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv",
                                   usecols=["Entrez id", "Ensembl id"])
    sysid_candidates.columns = ["entrez_id", "ensemble_id"]
    princeton_negative = pd.read_csv(ROOT_DIR + "ASD_translated_to_ensembl.csv")["entrez_id"].to_list()

    morbid_gene_symbols_list = pd.read_csv(ROOT_DIR + "MorbidGenes-Panel-v5_2020-08-26_for_varvis.csv",
                                           header=None).iloc[:, 0].to_list()
    all_genes_df = pd.read_csv(ROOT_DIR + "hgnc_protein_coding.tsv", index_col=False,
                               usecols=["entrez_id", "gene_symbol"], sep="\t",
                               dtype={"entrez_id": "Int32", "gene_symbol": str})
    morbid_genes = all_genes_df.loc[all_genes_df.gene_symbol.isin(morbid_gene_symbols_list)][
        "entrez_id"].dropna().to_list()
    panel_app_genes = pd.read_csv(ROOT_DIR + "Intellectual disability.tsv", sep="\t", usecols=["HGNC"])[
        "HGNC"].to_list()
    panel_app_genes = [int(x.strip("HGNC:")) for x in panel_app_genes if type(x) == str]
    g2p_dd_list = pd.read_csv(ROOT_DIR + "DDG2P_26_10_2020.csv",
                              usecols=["hgnc id"],
                              dtype={"hgnc id": "Int32"})
    g2p_dd_list = g2p_dd_list["hgnc id"].to_list()
    negative_gene_list = list(set(morbid_genes) - set(panel_app_genes)
                              - set(sysid_primary.entrez_id.to_list() + sysid_candidates.entrez_id.to_list())
                              - set(g2p_dd_list))

    if validation_run:
        negative_gene_list = random.sample(negative_gene_list, round(0.8 * len(negative_gene_list)))
        rows_id = random.sample(range(0, len(sysid_primary)), round(0.8 * len(sysid_primary)))
        sysid_primary = sysid_primary.loc[rows_id, :].reset_index(drop=True)
        rows_id = random.sample(range(len(sysid_candidates)), round(0.8 * len(sysid_candidates)))
        sysid_candidates = sysid_candidates.loc[rows_id, :].reset_index(drop=True)
        del rows_id
    del morbid_gene_symbols_list, validation_run, morbid_genes, panel_app_genes, g2p_dd_list, \
        princeton_negative

    n_cores = psutil.cpu_count()
    # MGI(n_cores, download=False).update()
    # StringDB(n_cores)
    # PsyMuKB()
    # GTEx()
    # Disgenet(n_cores, download=False)
    # PubtatorCentral(n_cores=46, download=False, preprocess=False, n_bootstraps=100000)
    fuse_data(validation_run=False)

    # evaluate_pubtator_parameters()
    # evaluate_mgi_parameters()