from concurrent.futures.process import ProcessPoolExecutor
from statistics import mean

import pandas as pd
from AutoCaSc_core.AutoCaSc import AutoCaSc, AUTOCASC_VERSION
import re
from concurrent.futures import ThreadPoolExecutor
import xlrd
import time

def calculate_comphet(variant, variant_other, inheritance, family_history):
    instance = AutoCaSc(variant, inheritance, family_history, other_impact="unknkown")
    other_instance = AutoCaSc(variant_other, inheritance, family_history, other_impact=instance.impact)
    instance.other_impact = other_instance.impact
    instance.get_scores()
    if instance.status_code == 200 and other_instance.status_code == 200:
        instance.candidate_score = mean([instance.candidate_score, other_instance.candidate_score])
    else:
        instance.status_code == 123
    return instance

def dot2comma(x):
    return str(x).replace(".", ",")


def export_results(instance, df_instance, i, strand_shift):
    if instance.status_code == 200:
        try:
            instance.calculate_candidate_score()
            df_instance.loc[i, "AutoCaSc_gene_symbol"] = dot2comma(instance.gene_symbol)
            df_instance.loc[i, "AutoCaSc_candidate_score"] = dot2comma(instance.candidate_score)
            df_instance.loc[i, "explanation"] = dot2comma("|".join([str(x) for x in instance.factors]))
            df_instance.loc[i, "AutoCaSc_literature_score"] = dot2comma(instance.literature_score)
            df_instance.loc[i, "AutoCaSc_mgi_score"] = dot2comma(instance.mgi_score)
            df_instance.loc[i, "AutoCaSc_pubtator_score"] = dot2comma(instance.pubtator_score)
            df_instance.loc[i, "AutoCaSc_gtex_score"] = dot2comma(instance.gtex_score)
            df_instance.loc[i, "AutoCaSc_denovo_rank_score"] = dot2comma(instance.denovo_rank_score)
            df_instance.loc[i, "AutoCaSc_disgenet_score"] = dot2comma(instance.disgenet_score)
            df_instance.loc[i, "AutoCaSc_string_score"] = dot2comma(instance.string_score)
            df_instance.loc[i, "strand_shift"] = dot2comma(strand_shift)
            df_instance.loc[i, "impact"] = instance.impact
            df_instance.loc[i, "AutoCaSc_version"] = dot2comma(AUTOCASC_VERSION)
        except AttributeError:
            for key in ["AutoCaSc_gene_symbol",
                        "AutoCaSc_candidate_score",
                        "AutoCaSc_literature_score",
                        "AutoCaSc_mgi_score",
                        "AutoCaSc_pubtator_score",
                        "AutoCaSc_gtex_score",
                        "AutoCaSc_denovo_rank_score",
                        "AutoCaSc_disgenet_score",
                        "strand_shift",
                        "AutoCaSc_string_score"]:
                df_instance.loc[i, key] = "-"
            df_instance.loc[i, "AutoCaSc_status_code"] = dot2comma(instance.status_code)
            df_instance.loc[i, "AutoCaSc_version"] = dot2comma(AUTOCASC_VERSION)

    else:
        for key in ["AutoCaSc_gene_symbol",
                    "AutoCaSc_candidate_score",
                    "AutoCaSc_literature_score",
                    "AutoCaSc_mgi_score",
                    "AutoCaSc_pubtator_score",
                    "AutoCaSc_gtex_score",
                    "AutoCaSc_denovo_rank_score",
                    "AutoCaSc_disgenet_score",
                    "strand_shift",
                    "AutoCaSc_string_score"]:
            df_instance.loc[i, key] = "-"
        df_instance.loc[i, "AutoCaSc_status_code"] = dot2comma(instance.status_code)
        df_instance.loc[i, "AutoCaSc_version"] = dot2comma(AUTOCASC_VERSION)
    return df_instance

def score_list(param_tuple):
    index_list, chunk_num = param_tuple

    other_variant = None
    strand_shift = False

    chunk_df = df.loc[index_list, :].reset_index(drop=True)
    for i in range(len(chunk_df)):
        print(f"{i} of {len(chunk_df)}")
        inheritance = chunk_df.loc[i, "inheritance"].lower()
        has_sibling = chunk_df.loc[i, "has_sibling"]
        cosegregating = chunk_df.loc[i, "cosegregating"]

        if ":" in str(chunk_df.loc[i, "hgvsc"]):
            print("found : in hgvsc!!! --> message from line 222")
            # chunk_df.loc[i, "hgvsc"] = chunk_df.loc[i, "hgvsc"].split(":")[-1]

        variant = chunk_df.loc[i, "transcript"].strip() + ":" + chunk_df.loc[i, "hgvsc"].strip()

        if inheritance == "comphet":
            other_variant = str(chunk_df.loc[i, "transcript"]).strip() + ":" + str(chunk_df.loc[i, "hgvsc_other"]).strip()

        instance = AutoCaSc(variant=variant,
                            inheritance=inheritance,
                            other_variant=other_variant,
                            has_sibling=has_sibling,
                            cosegregating=cosegregating)

        if instance.status_code == 200:
            chunk_df = export_results(instance=instance,
                                      df_instance=chunk_df,
                                      i=i,
                                      strand_shift=strand_shift)
        else:
            print("transcript not found")
            variant = str(chunk_df.loc[i, "gene_symbol_varvis"]).strip() + ":" + str(chunk_df.loc[i, "hgvsc"]).strip()
            if inheritance == "comphet":
                other_variant = str(chunk_df.loc[i, "gene_symbol_varvis"]).strip() + ":" + str(
                    chunk_df.loc[i, "hgvsc_other"]).strip()

            instance = AutoCaSc(variant=variant,
                                inheritance=inheritance,
                                other_variant=other_variant,
                                has_sibling=has_sibling,
                                cosegregating=cosegregating)

            if instance.status_code == 200:
                print("found with gene_symbol!")
                chunk_df = export_results(instance=instance,
                                          df_instance=chunk_df,
                                          i=i,
                                          strand_shift=strand_shift)
            else:
                print("gene not found")
                refseq, altseq = None, None
                try:
                    if re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()) is not None:
                        altseq = re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()).group()
                        if re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc"].strip()) is not None:
                            refseq = re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc"].strip()).group()
                    # if (refseq is not None) and (altseq is not None):
                            variant = ":".join([str(chunk_df.loc[i, "chr"]).strip(),
                                                str(chunk_df.loc[i, "pos_start"]).strip().replace(",", ""),
                                                refseq,
                                                altseq])
                except AttributeError:
                    print(f"\n\nproblems with {variant}\n\n")

                if inheritance == "comphet":
                    try:
                        refseq_other, altseq_other = None, None
                        if re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()) is not None:
                            altseq_other = re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc_other"].strip()).group()
                            if re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc_other"].strip()) is not None:
                                refseq_other = re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc_other"].strip()).group()
                            # if (refseq_other is not None) and (altseq_other is not None):
                                other_variant = ":".join([str(chunk_df.loc[i, "chr"]).strip(),
                                                          str(chunk_df.loc[i, "pos_start_other"]).strip().replace(",", ""),
                                                          refseq_other,
                                                          altseq_other])
                                instance = AutoCaSc(variant=variant,
                                                    inheritance=inheritance,
                                                    other_variant=other_variant,
                                                    has_sibling=has_sibling,
                                                    cosegregating=cosegregating)
                    except AttributeError:
                        print(f"\n\nproblems with {variant}\n\n")

                if instance.status_code == 201:
                    corresponding_base = {"C": "G",
                                          "G": "C",
                                          "T": "A",
                                          "A": "T",
                                          "": ""}

                    if (refseq is not None) and (altseq is not None):
                        variant = ":".join([chunk_df.loc[i, "chr"].strip(),
                                            chunk_df.loc[i, "pos_start"].strip(),
                                            corresponding_base[refseq],
                                            corresponding_base[altseq]])
                    if inheritance == "comphet":
                        if (refseq_other is not None) and (altseq_other is not None):
                            other_variant = ":".join([chunk_df.loc[i, "chr"].strip(),
                                                      chunk_df.loc[i, "pos_start_other"].strip(),
                                                      corresponding_base[refseq_other],
                                                      corresponding_base[altseq_other]])
                    instance = AutoCaSc(variant=variant,
                                        inheritance=inheritance,
                                        other_variant=other_variant,
                                        has_sibling=has_sibling,
                                        cosegregating=cosegregating)

                if instance.status_code == 200:
                    print("VCF found!")
                chunk_df = export_results(instance=instance,
                                          df_instance=chunk_df,
                                          i=i,
                                          strand_shift=strand_shift)
    return chunk_num, chunk_df

def inheritance_function(row):
    segregation = row["segregation"].lower()
    zygosity = row["zygosity"].lower()
    if "de novo" in segregation:
        inheritance = "de_novo"
    elif "het" in zygosity and ("maternal" in segregation or "paternal" in segregation) and "comphet" not in zygosity:
        inheritance = "ad_inherited"
    elif zygosity == "homo":
       inheritance = "homo"
    elif "comphet" in zygosity:
       inheritance = "comphet"
    elif "hemi" in zygosity:
       inheritance = "x_linked"
    else:
       inheritance = "unknown"
    return inheritance

def family_history_function(row):
    has_sibling = False
    cosegregating = False
    family_history = str(row["family_history"]).lower()
    if any([a in family_history for a in ["sister", "brother", "sibling"]]):
        has_sibling = True
        cosegregating = True
    return has_sibling, cosegregating

def process_excel_df(df):
    fill_unknown = ['zygosity', "segregation", "hgvsc", "chr", "pos_start", "pos_end", "transcript"]
    df.loc[:, fill_unknown] = df.loc[:, fill_unknown].fillna("unknown")
    df.fillna("", inplace=True)

    df.loc[:, "hgvsc"] = df.loc[:, "hgvsc"].apply(lambda x: x.split(":")[-1] if ":" in x else x)
    df.loc[:, "hgvsc_other"] = df.loc[:, "hgvsc_other"].apply(
        lambda x: x.split(":")[-1] if isinstance(x, str) and ":" in x else x)

    df["manual_CaSc"] = df["manual_CaSc"].astype(str).str.replace(".",",")
    df["inheritance"] = [inheritance_function(row) for _, row in df.iterrows()]
    df[["has_sibling", "cosegregating"]] = [family_history_function(row) for _, row in
                                                      df.iterrows()]
    return df


# column_dict_old = {
#     # "variant_id":"variant_id",
#     "HGNC symbol of candidate gene\n\n(alternative gene name)": 'gene_symbol',
#     'Zygosity\n\nhet /\nhomo /\nhemi /\nmosaic': 'zygosity',
#     'Origin\n\nde novo / paternal & maternal': 'segregation',
#     'Family history': 'family_history',
#     'hints to be involved in neuronal functions \nno= 0\nyes= 0,5\nsignalling/development= 1\n': 'neuronal_functions',
#     'animal models with neuronal-phenotype \nno=0, \nyes= 0,5\nyes and fits=1': 'animal_model',
#     'The gene has been reported in another database\nHGMD (handle restrictively), DDD, personal communication, GeneMatcher) as a candidate for a related phenotype (autism, epilepsy, neuronal issues, etc,\n0,33 per hit, max value is 2': 'reported_elsewhere',
#     'related genes with neuronal phenotype, or interactions with genes that has neuronal influence?\nyes: 1': 'interactions',
#     'expression\nnot in cns= 0\nlow in cns, and more in other tissues= 0,4\nin several tissues, including cns= 0,7\nmost in cns= 1': 'expression',
#     'Impact': 'impact',
#     'Chr': 'chr',
#     'Pos\n\nStart\n\n(formerly it was just one specification for the position)': 'pos_start',
#     'Pos\n\nEnd': 'pos_end',
#     'Transcript': 'transcript',
#     'cDNA': 'hgvsc',
#     'AAChange': 'hgvsp',
#     'Gene': 'gene',
#     'Chr.1': 'chr_other',
#     'Pos\n\nStart': 'pos_start_other',
#     'Pos\n\nEnd.1': 'pos_end_other',
#     'Transcript.1': 'transcript_other',
#     'cDNA.1': 'hgvsc_other',
#     'Impact.1': 'impact_other'
# }
column_dict = {
    "HGNC symbol": 'gene_symbol',
    "variant 1": "variant_1",
    "variant 2": "variant_2",
    "max, score: 15": "manual_CaSc",
    'Zygosity\n\nhet /\nhomo /\nhemi /\nmosaic': 'zygosity',
    'Origin\n\nde novo / paternal & maternal / paternal / maternal / unknown\n\nBitte sonst keine anderen Angaben. Weitere Überlegungen kommen unter Kommentare oder Thoughts on Research': 'segregation',
    'Family history': 'family_history',
    'Variant 1\n\nChr': 'chr',
    'Pos\n\nStart': 'pos_start',
    'Pos\n\nEnd': 'pos_end',
    'Transcript': 'transcript',
    'cDNA': 'hgvsc',
    'Gene': 'gene_symbol_varvis',
    'Pos\n\nStart.1': 'pos_start_other',
    'Pos\n\nEnd.1': 'pos_end_other',
    'cDNA.1': 'hgvsc_other'
}

df = pd.read_excel("/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/sonstige/data/candidate-scores_20201124_essential.xlsx",
                   dtype=str,
                   usecols=list(column_dict.keys()),
                   # nrows=100
                   )
df.rename(columns=column_dict, inplace=True)

# df = df.iloc[:50].reset_index(drop=True)

df = process_excel_df(df)

new_scores = pd.DataFrame()
# problem_df = pd.DataFrame(columns=df.columns)

n_chunks = 5
# n_chunks = 1
index_list = [list(range(len(df)))[round(len(df) / n_chunks * i):round(len(df) / n_chunks * (i + 1))] for i in range(n_chunks)]
chunk_num = list(range(n_chunks))
with ThreadPoolExecutor() as executor:
    results = executor.map(score_list, zip(index_list, chunk_num))
    result_df = pd.DataFrame()
print("processing done!")
result_dict = {chunk_num: df_chunk for chunk_num, df_chunk in results}

for _num_part in range(n_chunks):
    print(f"part {_num_part}")
    result_df = pd.concat([result_df, result_dict.get(_num_part)])
    # problem_df = pd.concat([problem_df, problem_chunk])

result_df = result_df.drop_duplicates()
result_df.to_csv("/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/sonstige/data/candidate-scores_reevaluated.csv", index=False)
# problem_df.to_csv("/media/johann/ROTHSTICK/candidate_table/candidate-scores_problems.xlsx", index=False)


print("done")
