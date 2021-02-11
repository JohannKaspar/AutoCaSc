import os
import re

import numpy as np
import shlex
import subprocess

import pandas as pd
from io import StringIO


def load_omim_morbid(path="/home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/OMIM_morbidmap.tsv"):
    # loading OMIM morbid dump
    with open(path, "r") as raw_file:
        filtered_file = ""
        for line in raw_file:
            if not line[0] in ["#", "[", "{", "?"]:
                filtered_file += line

    omim_morbid = pd.read_csv(StringIO(filtered_file),
                              sep="\t",
                              header=None,
                              usecols=range(3))
    omim_morbid.columns = ["disease", "gene_symbol", "mim_gene_number"]
    omim_morbid.gene_symbol = omim_morbid.gene_symbol.apply(lambda x: x.split(",")[0] if "," in x else x)
    try:
        omim_morbid["mim_number"] = omim_morbid.disease.apply(lambda x: re.findall(", [0-9]{6}", x)[0].strip(", ") if re.findall(", [0-9]{6}", x) != [] else None)
    except IndexError:
        print("error indexing")
    omim_morbid = omim_morbid.dropna()
    return omim_morbid


def get_mim_number(gene, omim_morbid=None):
    omim_gene = omim_morbid.loc[omim_morbid.gene_symbol == gene]
    if omim_gene.empty:
        return ""
    else:
        return str(omim_gene.mim_number.to_list())


def mim_map(df, column="gene_symbol", omim_morbid=None):
    if omim_morbid is None:
        omim_morbid = load_omim_morbid()

    for i in range(len(df)):
        _gene = df.loc[i, column]
        mim_number = get_mim_number(_gene, omim_morbid)
        df.loc[i, "mim_number"] = mim_number

    return df

def add_ranks(df):
    for version in ["v1", "v2", "v3"]:
        df[f"candidate_score_{version}"] = pd.to_numeric(df[f"candidate_score_{version}"],
                                                         errors="coerce")
        df.sort_values(f"candidate_score_{version}",
                       ascending=False,
                       inplace=True,
                       ignore_index=True)
        df.loc[:, f"rank_{version}"] = df.index
        df.loc[:, f"rank_{version}"] = df.loc[:, f"rank_{version}"].apply(lambda x: int(x+1))

        temp = filter_ac_impact_mim(df)
        temp.sort_values(f"candidate_score_{version}",
                         ascending=False,
                         inplace=True,
                         ignore_index=True)
        temp.loc[:, f"rank_{version}_filtered"] = temp.index
        temp.loc[:, f"rank_{version}_filtered"] = temp.loc[:, f"rank_{version}_filtered"].apply(lambda x: int(x+1))

        df = df.merge(temp[[f"rank_{version}", f"rank_{version}_filtered"]],
                      on=f"rank_{version}",
                      how="left")

        df.loc[:, f"rank_{version}_filtered"] = pd.to_numeric(df.loc[:, f"rank_{version}_filtered"],
                                                                         downcast="unsigned",
                                                                         errors="ignore")
    return df

def filter_ac_impact_mim(df):
    temp = pd.concat([df.loc[df.sysid != ""], df.loc[df.mim_number == ""]]).drop_duplicates()
    temp = temp.loc[temp.impact.isin(["moderate", "high"])]
    temp.AC = pd.to_numeric(temp.AC, errors="coerce", downcast="unsigned").fillna(1)

    temp = pd.concat(
        [
            temp.loc[(temp.inheritance == "de_novo") & (temp.AC < 2)],
            temp.loc[(temp.inheritance == "homo") & (temp.AC < 5)],
            temp.loc[(temp.inheritance == "comphet") & (temp.AC < 3)],
            temp.loc[(temp.inheritance == "x_linked") & (temp.AC < 3)],
            temp.loc[(temp.inheritance == "ad_inherited") & (temp.AC < 3)],
        ],
        ignore_index=True)
    return temp
    #temp.reset_index(drop=True, inplace=True)


def in_sysid(df):
    sysid_primary = pd.read_csv("/home/johann/PycharmProjects/AutoCaSc_project_folder/"
                                "sonstige/data/sysid_primary_20210203.csv")["Gene symbol"].to_list()
    sysid_candidates = pd.read_csv("/home/johann/PycharmProjects/AutoCaSc_project_folder/"
                                   "sonstige/data/sysid_candidates_20210203.csv")["Gene symbol"].to_list()

    for i, row in df.iterrows():
        gene_symbol = row["gene_symbol"]
        if gene_symbol in sysid_primary:
            df.loc[i, "sysid"] = "primary"
        elif gene_symbol in sysid_candidates:
            df.loc[i, "sysid"] = "candidates"
        else:
            df.loc[i, "sysid"] = ""

    return df

def load_blacklist(path="/home/johann/PycharmProjects/AutoCaSc_project_folder/sonstige/data/gene_blacklist.txt"):
    blacklist = []
    with open(path, "r") as file:
        for line in file:
            if not line[0] == "#":
                blacklist.append(line.strip())
    return blacklist

def filter_blacklist(df, blacklist):
    for pattern in blacklist:
        df = df.loc[~df.gene_symbol.str.contains(pattern,
                                                regex=True,
                                                na=False)]
    return df.reset_index(drop=True)

def convert_to_int(x):
    if x == "-":
        return x
    try:
        return int(float(x))
    except ValueError:
        return x

blacklist = load_blacklist()
omim_morbid = load_omim_morbid()
if not os.path.exists("/home/johann/trio_scoring_results/varvis_trios/mim_mapped/"):
    os.mkdir("/home/johann/trio_scoring_results/varvis_trios/mim_mapped/")

varvis_trio_results = pd.DataFrame()
for entry in os.scandir("/home/johann/trio_scoring_results/varvis_trios/"):
    if entry.is_file() and not "lock" in entry.name:
        df = pd.read_csv(entry.path,
                         decimal=",",
                         sep="\t")

        try:
            df_whitelisted = filter_blacklist(df, blacklist)
            mim_mapped = mim_map(df_whitelisted, omim_morbid=omim_morbid)
            mim_mapped_sysid = in_sysid(mim_mapped)
            mim_mapped_sysid_ranks = add_ranks(mim_mapped_sysid)

            mim_mapped_sysid_ranks.sort_values("rank_v1_filtered",
                                               ascending=True,
                                               inplace=True)
        except AttributeError:
            print(f"some Attribute Error occured for {entry.name}")
        family_id = entry.name.split(".")[0]
        mim_mapped_sysid_ranks.loc[:, "family_id"] = family_id

        if not "other_variant" in mim_mapped_sysid_ranks.columns:
            mim_mapped_sysid_ranks.loc[:, "other_variant"] = ""

        for x in ["DP_moth", "AD_moth"]:
            if x in mim_mapped_sysid_ranks.columns:
                mim_mapped_sysid_ranks = mim_mapped_sysid_ranks.rename(columns={x: f"{x}er"})
        mim_mapped_sysid_ranks = mim_mapped_sysid_ranks[["family_id", 'variant', 'gene_symbol', 'hgvsc', 'hgvsp', 'impact', 'inheritance',
                                          'candidate_score_v1', 'candidate_score_v2', 'candidate_score_v3',
                                          'rank_v1_filtered', 'rank_v2_filtered',
                                          'rank_v3_filtered', 'literature_score', 'CADD_phred', 'sysid',
                                          'rank_v1', 'rank_v2', 'rank_v3', 'transcript', 'status_code', 'AC',
                                                         #'AF',
                                          'QUAL', 'GQ_index', 'AD_index', 'AD_father', 'AD_mother', 'DP_index',
                                          'DP_father', 'DP_mother', 'factors', 'other_variant', 'mim_number']]
        mim_mapped_sysid_ranks.to_csv(f"/home/johann/trio_scoring_results/varvis_trios/mim_mapped/{entry.name}.mim",
                          sep="\t",
                          index=False,
                          decimal=",")
        varvis_trio_results = pd.concat([varvis_trio_results, mim_mapped_sysid_ranks],
                                        ignore_index=True)


for column in ["QUAL", "AC", "status_code", "GQ_index", "DP_index", "DP_father", "DP_mother"]:
    varvis_trio_results.loc[:, column] = varvis_trio_results.loc[:, column].apply(convert_to_int)

for column in varvis_trio_results.columns:
    if "score" in column or column == "CADD_phred":
        varvis_trio_results.loc[:, column] = varvis_trio_results.loc[:, column].astype(float)

varvis_trio_results = varvis_trio_results.convert_dtypes()
varvis_trio_results.to_csv(f"/home/johann/trio_scoring_results/varvis_trios/mim_mapped/varvis_trio_results_2021-02-10.csv",
                          sep="\t",
                          index=False,
                          decimal=",")

#with open("/home/johann/trio_scoring_results/varvis_trios/mim_mapped/varvis_trio_results_20210205.csv", "w") as outfile:
#    subprocess.run(shlex.split("awk '(NR == 1) || (FNR > 1)' /home/johann/trio_scoring_results/varvis_trios/mim_mapped/L*.csv.mim"),
#                   universal_newlines=True,
#                   stdout=outfile)
# for file in *.csv.mim; do sed -i "s/$/\t$file/" $file; done
# awk '(NR == 1) || (FNR > 1)' L*.csv.mim > concat.csv

