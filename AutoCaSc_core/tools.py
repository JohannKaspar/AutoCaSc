import random

import pandas as pd

# ROOT_DIR = '/home/johann/AutoCaSc_core/data/'
ROOT_DIR = "/Users/johannkaspar/OneDrive/Promotion/AutoCaSc_project_folder/AutoCaSc_maintenance/data/"



def filterTheDict(dictObj, filter_value, key_to_pop=None):
  newDict = dict()
  for (key, value) in dictObj.items():
    if value == filter_value:
      newDict[key] = value
  if len(newDict) > 2:
    print("\n\n\nATTENTION! More than 2 entries in the comphet dict!!!\n\n\n")
  if key_to_pop:
    newDict.pop(key_to_pop)
  return newDict

def safe_get(l, idx, default=None):
  try:
    return l[idx]
  except IndexError:
    return default
  except AttributeError:
    return default

def rank_genes(df, poi="gene_score", exponent=1):
  gene_scores_df = df.sort_values(by=poi, ignore_index=True, ascending=False)
  value_counts = gene_scores_df[poi].value_counts().reset_index()
  value_counts.columns = [poi, "count"]
  value_counts = value_counts.sort_values(by=poi, ascending=True, ignore_index=True)
  # rank_scores = np.linspace(0., 1., len(value_counts))
  # total = value_counts.count.sum()

  for n in range(len(value_counts)):
    value_counts.loc[n,"rank_score"] = value_counts.iloc[:n].loc[:,"count"].sum() / value_counts["count"].sum()

  gene_scores_df = gene_scores_df.merge(value_counts[[poi, "rank_score"]], how="left", left_on=poi, right_on=poi)

  gene_scores_df.rank_score = gene_scores_df.rank_score ** exponent

  return gene_scores_df

def lin_rank(df, poi, mode="median"):
  if mode == "median":
    df["gene_score"] = df[poi] - df.loc[df.sys == 0][poi].median()
    df.loc[df.gene_score < 0, "gene_score"] = 0
    df.gene_score = df.gene_score / df.loc[df.sys_primary == 1].gene_score.median()
    df.loc[df.gene_score > 1, "gene_score"] = 1
  elif mode == "mean":
    df["gene_score"] = df[poi] - df.loc[df.sys == 0][poi].mean()
    df.loc[df.gene_score < 0, "gene_score"] = 0
    df.gene_score = df.gene_score / df.loc[df.sys_primary == 1].gene_score.mean()
    df.loc[df.gene_score > 1, "gene_score"] = 1
  return df


def add_categories(original_df, column, identifiers, negative_list="princeton"):
  ROOT_DIR = "/home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_maintenance/data/"
  # ROOT_DIR = "/Users/johannkaspar/OneDrive/Promotion/AutoCaSc_project_folder/AutoCaSc_maintenance/data/"

  sysid_primary = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv",
                              usecols=["Entrez id", "Ensembl id"])
  sysid_primary.columns = ["entrez_id", "ensemble_id"]
  sysid_candidates = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv",
                                 usecols=["Entrez id", "Ensembl id"])
  sysid_candidates.columns = ["entrez_id", "ensemble_id"]
  princeton_negative = pd.read_csv(ROOT_DIR + "ASD_translated_to_ensembl.csv", usecols=["entrez_id", "ensemble_id"])

  sysid_primary_gene_symbols = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv", usecols=["Gene symbol"])[
    "Gene symbol"].to_list()
  sysid_candidates_gene_symbols = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv", usecols=["Gene symbol"])[
    "Gene symbol"].to_list()
  morbid_gene_symbols_list = list(set(pd.read_csv("/home/johann/AutoCaSc/data/pubtator_central/MorbidGenes-Panel"
                                                  "-v5_2020-08-26_for_varvis.csv", header=None).iloc[:, 0].to_list()) -
                                  set(sysid_primary_gene_symbols + sysid_candidates_gene_symbols))
  all_genes_df = pd.read_csv(ROOT_DIR + "hgnc_protein_coding.tsv", index_col=False,
                             usecols=["entrez_id", "gene_symbol", "ensemble_id"], sep="\t",
                             dtype={"entrez_id": "Int32", "gene_symbol": str})
  morbid_genes = all_genes_df.loc[all_genes_df.gene_symbol.isin(morbid_gene_symbols_list)][["entrez_id", "ensemble_id"]]

  del all_genes_df, morbid_gene_symbols_list, sysid_candidates_gene_symbols, sysid_primary_gene_symbols

  dr_ranked = original_df.copy()
  if identifiers == "ensemble":
      sysid_primary = sysid_primary["ensemble_id"].to_list()
      sysid_candidates = sysid_candidates["ensemble_id"].to_list()
      if negative_list == "princeton":
        negative_genes = princeton_negative["ensemble_id"].to_list()
      else:
        negative_genes = morbid_genes["ensemble_id"].to_list()

  else:
      sysid_primary = sysid_primary["entrez_id"].to_list()
      sysid_candidates = sysid_candidates["entrez_id"].to_list()
      if negative_list == "princeton":
        negative_genes = princeton_negative["entrez_id"].to_list()
      else:
        negative_genes = morbid_genes["entrez_id"].to_list()

  dr_ranked["sys_primary"] = dr_ranked[column].isin(sysid_primary).astype(int)
  dr_ranked["sys_candidate"] = dr_ranked[column].isin(sysid_candidates).astype(int)
  dr_ranked["negative"] = dr_ranked[column].isin(negative_genes).astype(int)
  dr_ranked["sys"] = dr_ranked[column].isin(sysid_primary + sysid_candidates).astype(int)


  dr_ranked["sys_category"] = "unknown"
  dr_ranked.loc[dr_ranked.sys_candidate == 1, "sys_category"] = "candidate"
  dr_ranked.loc[dr_ranked.sys_primary == 1, "sys_category"] = "known NDD"
  dr_ranked.loc[dr_ranked.negative == 1, "sys_category"] = "negative control"

  return dr_ranked


def negative_product(score_list):
  temp = 1
  for score in score_list:
    temp *= (1 - score)
  return 1 - temp