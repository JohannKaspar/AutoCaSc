from concurrent.futures.process import ProcessPoolExecutor
from statistics import mean

import pandas as pd
from AutoCaSc_core.AutoCaSc import AutoCaSc, AUTOCASC_VERSION
import re
from concurrent.futures import ThreadPoolExecutor
import xlrd
import time

n_cores = 1
column_dict = {
    # "variant_id":"variant_id",
    "HGNC symbol of candidate gene\n\n(alternative gene name)": 'gene_symbol',
    'Zygosity\n\nhet /\nhomo /\nhemi /\nmosaic': 'zygosity',
    'Origin\n\nde novo / paternal & maternal': 'segregation',
    'Family history': 'family_history',
    'hints to be involved in neuronal functions \nno= 0\nyes= 0,5\nsignalling/development= 1\n': 'neuronal_functions',
    'animal models with neuronal-phenotype \nno=0, \nyes= 0,5\nyes and fits=1': 'animal_model',
    'The gene has been reported in another database\nHGMD (handle restrictively), DDD, personal communication, GeneMatcher) as a candidate for a related phenotype (autism, epilepsy, neuronal issues, etc,\n0,33 per hit, max value is 2': 'reported_elsewhere',
    'related genes with neuronal phenotype, or interactions with genes that has neuronal influence?\nyes: 1': 'interactions',
    'expression\nnot in cns= 0\nlow in cns, and more in other tissues= 0,4\nin several tissues, including cns= 0,7\nmost in cns= 1': 'expression',
    'Impact': 'impact',
    'Chr': 'chr',
    'Pos\n\nStart\n\n(formerly it was just one specification for the position)': 'pos_start',
    'Pos\n\nEnd': 'pos_end',
    'Transcript': 'transcript',
    'cDNA': 'hgvsc',
    'AAChange': 'hgvsp',
    'Gene': 'gene',
    'Chr.1': 'chr_other',
    'Pos\n\nStart': 'pos_start_other',
    'Pos\n\nEnd.1': 'pos_end_other',
    'Transcript.1': 'transcript_other',
    'cDNA.1': 'hgvsc_other',
    'Impact.1': 'impact_other'
}
df = pd.read_excel("/Volumes/ROTHSTICK/candidate_table/candidate-scores_clean.xlsx",
                   sheet_name="CaSc", dtype=str,
                   # usecols=list(column_dict.keys()),
                   # nrows=25
                   )
# _cols = list(df.columns)
# df = df.iloc[145:,:].reset_index(drop=True)
# df.columns = ["variant_id"] + _cols

# df = pd.read_csv("/home/johann/AutoCaSc/data/candidate_table/candidate-scores-20200810-relevant_2.csv",
#                  sep="\t",
#                  dtype=str,
#                  # encoding="UT",
#                  # skiprows=105,
#                  # nrows=10
#                  )
# df.columns = ['gene_symbol', 'percentile', 'candidate_score', 'zygosity',
#        'segregation', 'family_history', 'neuronal_functions', 'animal_model',
#        'reported_elsewhere', 'interactions', 'expression', 'literature_score',
#        'impact', 'chr', 'pos_start', 'pos_end', 'transcript', 'hgvsc', 'hgvsp',
#        'gene', 'chr_other', 'pos_start_othe', 'pos_end_other',
#        'transcript_other', 'hgvsc_other', 'impact_other']
# df = df.rename(columns=column_dict)

df.loc[:, ['zygosity', "segregation", "hgvsc", "chr", "pos_start", "pos_end", "transcript"]] = df.loc[:,
                                                                                               ['zygosity',
                                                                                                "segregation", "hgvsc",
                                                                                                "chr", "pos_start",
                                                                                                "pos_end",
                                                                                                "transcript"]].fillna(
    "unknown")
df.loc[:, "hgvsc"] = df.loc[:, "hgvsc"].apply(lambda x: x.split(":")[-1] if ":" in x else x)
df.loc[:, "hgvsc_other"] = df.loc[:, "hgvsc_other"].apply(lambda x: x.split(":")[-1] if isinstance(x, str) and ":" in x else x)

new_scores = pd.DataFrame()
# problem_df = pd.DataFrame(columns=df.columns)
print("stop")

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


def export_results(instance, df_instance, i):
    if instance.status_code == 200:
        df_instance.loc[i, "AutoCaSc_gene_symbol"] = dot2comma(instance.gene_symbol)
        df_instance.loc[i, "AutoCaSc_candidate_score"] = dot2comma(instance.candidate_score)
        df_instance.loc[i, "AutoCaSc_inheritance_score"] = dot2comma(instance.inheritance_score)
        df_instance.loc[i, "AutoCaSc_gene_attribute_score"] = dot2comma(instance.gene_attribute_score)
        df_instance.loc[i, "AutoCaSc_variant_attribute_score"] = dot2comma(instance.variant_score)
        df_instance.loc[i, "AutoCaSc_literature_score"] = dot2comma(instance.literature_score)
        df_instance.loc[i, "AutoCaSc_mgi_score"] = dot2comma(instance.mgi_score)
        df_instance.loc[i, "AutoCaSc_pubtator_score"] = dot2comma(instance.pubtator_score)
        df_instance.loc[i, "AutoCaSc_gtex_score"] = dot2comma(instance.gtex_score)
        df_instance.loc[i, "AutoCaSc_denovo_rank_score"] = dot2comma(instance.denovo_rank_score)
        df_instance.loc[i, "AutoCaSc_disgenet_score"] = dot2comma(instance.disgenet_score)
        df_instance.loc[i, "AutoCaSc_string_score"] = dot2comma(instance.string_score)
        df_instance.loc[i, "AutoCaSc_version"] = dot2comma(AUTOCASC_VERSION)

    else:
        for key in ["AutoCaSc_variant_attribute_score", "AutoCaSc_gene_attribute_score", "AutoCaSc_literature_score",
                    "AutoCaSc_candidate_score", "AutoCaSc_mgi_score", "AutoCaSc_pubtator_score", "AutoCaSc_gtex_score",
                    "AutoCaSc_denovo_rank_score", "AutoCaSc_disgenet_score", "AutoCaSc_string_score"]:
            df_instance.loc[i, key] = "-"
        df_instance.loc[i, "AutoCaSc_status_code"] = dot2comma(instance.status_code)
        df_instance.loc[i, "AutoCaSc_version"] = dot2comma(AUTOCASC_VERSION)
    return df_instance


# def score_variant(variant, inheritance, family_history, impact_other, chunk_df, i):
#     instance = AutoCaSc(variant=variant, inheritance=inheritance, family_history=family_history,
#                         other_impact=impact_other)
#
#     if instance.status_code in [200, 201]:
#         chunk_df = export_results(instance, chunk_df, i)
#     else:
#         print("transcript not found")
#         variant = str(chunk_df.loc[i, "gene_symbol"]).strip() + ":" + str(chunk_df.loc[i, "hgvsc"]).strip()
#         instance = AutoCaSc(variant, inheritance, family_history, impact_other)
#         if instance.status_code == 200:
#             print("found with gene_symbol!")
#             chunk_df = export_results(instance, chunk_df, i)
#         else:
#             print("gene not found")
#             refseq, altseq = "", ""
#             if re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc"].strip()) is not None:
#                 refseq = re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc"].strip()).group()
#             if re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()) is not None:
#                 altseq = re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()).group()
#
#             variant = ":".join([str(chunk_df.loc[i, "chr"]).strip(),
#                                 str(chunk_df.loc[i, "pos_start"]).strip().replace(",", ""),
#                                 refseq,
#                                 altseq])
#             instance = AutoCaSc(variant, inheritance, family_history, impact_other)
#             if instance.status_code == 201:
#                 corresponding_base = {"C": "G",
#                                       "G": "C",
#                                       "T": "A",
#                                       "A": "T",
#                                       "": ""}
#
#                 variant = ":".join([chunk_df.loc[i, "chr"].strip(),
#                                     chunk_df.loc[i, "pos_start"].strip(),
#                                     corresponding_base[refseq],
#                                     corresponding_base[altseq]])
#                 instance = AutoCaSc(variant, inheritance, family_history, impact_other)
#             if instance.status_code == 200:
#                 print("VCF found!")
#                 chunk_df = export_results(instance, chunk_df, i)
#             else:
#                 print("persisting problem")
#                 chunk_df = export_results(instance, chunk_df, i)
                # problem_df = pd.concat([problem_df, pd.DataFrame(chunk_df.iloc[i, :]).T])


# def helper_function(chunk_df, i):
#     segregation = chunk_df.loc[i, "segregation"].lower()
#     zygosity = chunk_df.loc[i, "zygosity"].lower()
#     if "de novo" in segregation:
#         inheritance = "de_novo"
#     elif "het" in zygosity and segregation == "unknown":
#         inheritance = "other"
#     elif "het" in zygosity and ("maternal" in segregation or "paternal" in segregation):
#         inheritance = "ad_inherited"
#     elif zygosity == "homo":
#         inheritance = "homo"
#     elif zygosity == "comphet":
#         inheritance = "comphet"
#     elif "hemi" in zygosity:
#         inheritance = "x_linked"
#     else:
#         inheritance = "unknown"
#
#     if ":" in str(chunk_df.loc[i, "hgvsc"]):
#         chunk_df.loc[i, "hgvsc"] = chunk_df.loc[i, "hgvsc"].split(":")[-1]
#
#     variant = chunk_df.loc[i, "transcript"].strip() + ":" + chunk_df.loc[i, "hgvsc"].strip()
#     family_history = chunk_df.loc[i, "family_history"] == "yes"
#     impact_other = "unknown"
#     strand_shift = False
#
#     instance = AutoCaSc(variant, inheritance, family_history, impact_other)
#
#     if inheritance == "comphet":
#         try:
#             variant_other = chunk_df.loc[i, "transcript"].strip() + ":" + chunk_df.loc[i, "hgvsc_other"].strip()
#             instance = calculate_comphet(variant, variant_other, inheritance, family_history)
#         except (AttributeError, IndexError, TypeError):
#             return
#
#     if instance.status_code == 200:
#         chunk_df = export_results(instance, chunk_df, i)
#     else:
#         print("transcript not found")
#         variant = str(chunk_df.loc[i, "gene_symbol"]).strip() + ":" + str(chunk_df.loc[i, "hgvsc"]).strip()
#         instance = AutoCaSc(variant, inheritance, family_history, impact_other)
#         if instance.status_code == 200:
#             print("found with gene_symbol!")
#             try:
#                 variant_other = chunk_df.loc[i, "gene_symbol"].strip() + ":" + chunk_df.loc[i, "hgvsc_other"].strip()
#                 instance = calculate_comphet(variant, variant_other, inheritance, family_history)
#             except (AttributeError, IndexError, TypeError):
#                 continue
#             chunk_df = export_results(instance, chunk_df, i)
#         else:
#             print("gene not found")
#             refseq, altseq = "", ""
#             if re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc"].strip()) is not None:
#                 refseq = re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc"].strip()).group()
#             if re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()) is not None:
#                 altseq = re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()).group()
#
#             variant = ":".join([str(chunk_df.loc[i, "chr"]).strip(),
#                                 str(chunk_df.loc[i, "pos_start"]).strip().replace(",", ""),
#                                 refseq,
#                                 altseq])
#             instance = AutoCaSc(variant, inheritance, family_history, impact_other)
#             if instance.status_code == 201:
#                 corresponding_base = {"C": "G",
#                                       "G": "C",
#                                       "T": "A",
#                                       "A": "T",
#                                       "": ""}
#
#                 variant = ":".join([chunk_df.loc[i, "chr"].strip(),
#                                     chunk_df.loc[i, "pos_start"].strip(),
#                                     corresponding_base[refseq],
#                                     corresponding_base[altseq]])
#                 instance = AutoCaSc(variant, inheritance, family_history, impact_other)
#                 strand_shift = True
#             if instance.status_code == 200:
#                 print("VCF found!")
#                 try:
#                     if re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc_other"].strip()) is not None:
#                         refseq_other = re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc_other"].strip()).group()
#                     if re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()) is not None:
#                         altseq_other = re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc_other"].strip()).group()
#
#                     variant_other = ":".join([str(chunk_df.loc[i, "chr"]).strip(),
#                                               str(chunk_df.loc[i, "pos_start_other"]).strip().replace(",", ""),
#                                               refseq_other,
#                                               altseq_other])
#                     if strand_shift:
#                         variant_other = ":".join([chunk_df.loc[i, "chr"].strip(),
#                                                   chunk_df.loc[i, "pos_start_other"].strip(),
#                                                   corresponding_base[refseq_other],
#                                                   corresponding_base[altseq_other]])
#                     instance = calculate_comphet(variant, variant_other, inheritance, family_history)
#                 except (AttributeError, IndexError, TypeError):
#                     continue
#                 chunk_df = export_results(instance, chunk_df, i)
#             else:
#                 print("persisting problem, starting second try in 30 seconds...")
#                 time.sleep(30)

def score_list(param_tuple):
    index_list, chunk_num = param_tuple
    chunk_df = df.loc[index_list, :].reset_index(drop=True)
    for i in range(len(chunk_df)):
        print(f"{i} of {len(chunk_df)}")
        segregation = chunk_df.loc[i, "segregation"].lower()
        zygosity = chunk_df.loc[i, "zygosity"].lower()
        if "de novo" in segregation:
            inheritance = "de_novo"
        elif "het" in zygosity and segregation == "unknown":
            inheritance = "other"
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

        if ":" in str(chunk_df.loc[i, "hgvsc"]):
            chunk_df.loc[i, "hgvsc"] = chunk_df.loc[i, "hgvsc"].split(":")[-1]

        variant = chunk_df.loc[i, "transcript"].strip() + ":" + chunk_df.loc[i, "hgvsc"].strip()
        family_history = chunk_df.loc[i, "family_history"] == "yes"
        impact_other = "unknown"
        strand_shift = False

        instance = AutoCaSc(variant, inheritance, family_history, impact_other)

        if inheritance == "comphet":
            try:
                variant_other = str(chunk_df.loc[i, "transcript"]).strip() + ":" + str(chunk_df.loc[i, "hgvsc_other"]).strip()
                instance = calculate_comphet(variant, variant_other, inheritance, family_history)
            except (AttributeError, IndexError, TypeError):
                continue

        if instance.status_code == 200:
            chunk_df = export_results(instance, chunk_df, i)
        else:
            print("transcript not found")
            variant = str(chunk_df.loc[i, "gene_symbol"]).strip() + ":" + str(chunk_df.loc[i, "hgvsc"]).strip()
            instance = AutoCaSc(variant, inheritance, family_history, impact_other)
            if instance.status_code == 200:
                print("found with gene_symbol!")
                if inheritance == "comphet":
                    try:
                        variant_other = chunk_df.loc[i, "gene_symbol"].strip() + ":" + chunk_df.loc[i, "hgvsc_other"].strip()
                        instance = calculate_comphet(variant, variant_other, inheritance, family_history)
                    except (AttributeError, IndexError, TypeError):
                        continue
                chunk_df = export_results(instance, chunk_df, i)
            else:
                print("gene not found")
                refseq, altseq = "", ""
                if re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc"].strip()) is not None:
                    refseq = re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc"].strip()).group()
                if re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()) is not None:
                    altseq = re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()).group()

                variant = ":".join([str(chunk_df.loc[i, "chr"]).strip(),
                                    str(chunk_df.loc[i, "pos_start"]).strip().replace(",", ""),
                                    refseq,
                                    altseq])
                instance = AutoCaSc(variant, inheritance, family_history, impact_other)
                if instance.status_code == 201:
                    corresponding_base = {"C": "G",
                                          "G": "C",
                                          "T": "A",
                                          "A": "T",
                                          "": ""}

                    variant = ":".join([chunk_df.loc[i, "chr"].strip(),
                                        chunk_df.loc[i, "pos_start"].strip(),
                                        corresponding_base[refseq],
                                        corresponding_base[altseq]])
                    instance = AutoCaSc(variant, inheritance, family_history, impact_other)
                    strand_shift = True
                if instance.status_code == 200:
                    print("VCF found!")
                    if inheritance == "comphet":
                        try:
                            if re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc_other"].strip()) is not None:
                                refseq_other = re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc_other"].strip()).group()
                            if re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()) is not None:
                                altseq_other = re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc_other"].strip()).group()

                            variant_other = ":".join([str(chunk_df.loc[i, "chr"]).strip(),
                                        str(chunk_df.loc[i, "pos_start_other"]).strip().replace(",", ""),
                                        refseq_other,
                                        altseq_other])
                            if strand_shift:
                                variant_other = ":".join([chunk_df.loc[i, "chr"].strip(),
                                          chunk_df.loc[i, "pos_start_other"].strip(),
                                          corresponding_base[refseq_other],
                                          corresponding_base[altseq_other]])
                            instance = calculate_comphet(variant, variant_other, inheritance, family_history)
                        except (AttributeError, IndexError, TypeError):
                            continue
                        # chunk_df = export_results(instance, chunk_df, i)
                chunk_df = export_results(instance, chunk_df, i)
                # else:
                #     print("persisting problem, starting second try in 30 seconds...")
                #     time.sleep(5)
                #     instance = AutoCaSc(variant, inheritance, family_history, impact_other)
                #
                #     if inheritance == "comphet":
                #         try:
                #             variant_other = chunk_df.loc[i, "transcript"].strip() + ":" + chunk_df.loc[
                #                 i, "hgvsc_other"].strip()
                #             instance = calculate_comphet(variant, variant_other, inheritance, family_history)
                #         except (AttributeError, IndexError, TypeError):
                #             continue
                #
                #     if instance.status_code == 200:
                #         chunk_df = export_results(instance, chunk_df, i)
                #     else:
                #         print("transcript not found")
                #         variant = str(chunk_df.loc[i, "gene_symbol"]).strip() + ":" + str(
                #             chunk_df.loc[i, "hgvsc"]).strip()
                #         instance = AutoCaSc(variant, inheritance, family_history, impact_other)
                #         if instance.status_code == 200:
                #             print("found with gene_symbol!")
                #             try:
                #                 variant_other = chunk_df.loc[i, "gene_symbol"].strip() + ":" + chunk_df.loc[
                #                     i, "hgvsc_other"].strip()
                #                 instance = calculate_comphet(variant, variant_other, inheritance, family_history)
                #             except (AttributeError, IndexError, TypeError):
                #                 continue
                #             chunk_df = export_results(instance, chunk_df, i)
                #         else:
                #             print("gene not found")
                #             refseq, altseq = "", ""
                #             if re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc"].strip()) is not None:
                #                 refseq = re.search(r"(?<=\d)[CTGA]+(?=>)", chunk_df.loc[i, "hgvsc"].strip()).group()
                #             if re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()) is not None:
                #                 altseq = re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()).group()
                #
                #             variant = ":".join([str(chunk_df.loc[i, "chr"]).strip(),
                #                                 str(chunk_df.loc[i, "pos_start"]).strip().replace(",", ""),
                #                                 refseq,
                #                                 altseq])
                #             instance = AutoCaSc(variant, inheritance, family_history, impact_other)
                #             if instance.status_code == 201:
                #                 corresponding_base = {"C": "G",
                #                                       "G": "C",
                #                                       "T": "A",
                #                                       "A": "T",
                #                                       "": ""}
                #
                #                 variant = ":".join([chunk_df.loc[i, "chr"].strip(),
                #                                     chunk_df.loc[i, "pos_start"].strip(),
                #                                     corresponding_base[refseq],
                #                                     corresponding_base[altseq]])
                #                 instance = AutoCaSc(variant, inheritance, family_history, impact_other)
                #                 strand_shift = True
                #             if instance.status_code == 200:
                #                 print("VCF found!")
                #                 try:
                #                     if re.search(r"(?<=\d)[CTGA]+(?=>)",
                #                                  chunk_df.loc[i, "hgvsc_other"].strip()) is not None:
                #                         refseq_other = re.search(r"(?<=\d)[CTGA]+(?=>)",
                #                                                  chunk_df.loc[i, "hgvsc_other"].strip()).group()
                #                     if re.search(r"(?<=>)[CTGA]+", chunk_df.loc[i, "hgvsc"].strip()) is not None:
                #                         altseq_other = re.search(r"(?<=>)[CTGA]+",
                #                                                  chunk_df.loc[i, "hgvsc_other"].strip()).group()
                #
                #                     variant_other = ":".join([str(chunk_df.loc[i, "chr"]).strip(),
                #                                               str(chunk_df.loc[i, "pos_start_other"]).strip().replace(
                #                                                   ",", ""),
                #                                               refseq_other,
                #                                               altseq_other])
                #                     if strand_shift:
                #                         variant_other = ":".join([chunk_df.loc[i, "chr"].strip(),
                #                                                   chunk_df.loc[i, "pos_start_other"].strip(),
                #                                                   corresponding_base[refseq_other],
                #                                                   corresponding_base[altseq_other]])
                #                     instance = calculate_comphet(variant, variant_other, inheritance, family_history)
                #                 except (AttributeError, IndexError, TypeError):
                #                     continue
                #                 chunk_df = export_results(instance, chunk_df, i)
                #             else:
                #                 print("persisting problem")
                #                 chunk_df = export_results(instance, chunk_df, i)
    return chunk_num, chunk_df

n_chunks = 10
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
result_df.to_csv("/Volumes/ROTHSTICK/candidate_table/candidate-scores_reevaluated.csv", index=False)
# problem_df.to_csv("/media/johann/ROTHSTICK/candidate_table/candidate-scores_problems.xlsx", index=False)


print("done")
