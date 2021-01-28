import os
from sys import path

import pandas as pd
from io import StringIO


def load_omim_morbid(path="/home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/OMIM_morbidmap.tsv"):
    # loading OMIM morbid dump
    with open(path, "r") as raw_file:
        filtered_file = ""
        for line in raw_file:
            if not line[0] in ["#", "[", "{"]:
                filtered_file += line

    omim_morbid = pd.read_csv(StringIO(filtered_file),
                              sep="\t",
                              header=None,
                              usecols=range(3))
    omim_morbid.columns = ["disease", "gene_symbol", "mim_number"]
    omim_morbid.gene_symbol = omim_morbid.gene_symbol.apply(lambda x: x.split(",")[0] if "," in x else x)
    return omim_morbid


def get_mim_number(gene, omim_morbid=None):
    omim_gene = omim_morbid.loc[omim_morbid.gene_symbol == gene]
    if omim_gene.empty:
        return ""
    else:
        return str(omim_gene.mim_number.to_list())


def process_df(df, column="gene_symbol", omim_morbid=None):
    if omim_morbid is None:
        omim_morbid = load_omim_morbid()

    for i in range(len(df)):
        _gene = df.loc[i, column]
        mim_number = get_mim_number(_gene, omim_morbid)
        df.loc[i, "mim_number"] = mim_number

    return df


omim_morbid = load_omim_morbid()
if not os.path.exists("/home/johann/trio_scoring_results/varvis_trios/mim_mapped/"):
    os.mkdir("/home/johann/trio_scoring_results/varvis_trios/mim_mapped/")

for entry in os.scandir("/home/johann/trio_scoring_results/varvis_trios/"):
    if entry.is_file():
        df = pd.read_csv(entry.path)
        mim_mapped = process_df(df, omim_morbid=omim_morbid)
        mim_mapped.to_csv(f"/home/johann/trio_scoring_results/varvis_trios/mim_mapped/{entry.name}.csv",
                          #sep="\t",
                          index=False)
