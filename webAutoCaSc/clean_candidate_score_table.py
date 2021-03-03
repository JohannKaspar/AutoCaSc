import pandas as pd
from numpy import nan
import xlrd
import os
import re
from AutoCaSc_core.AutoCaSc import AutoCaSc
from AutoCaSc_core.AutoCaSc_vcf import update_request_cache

def load_manual_candidate_table(path_to_excel="/Users/johannkaspar/Documents/Promotion/AutoCaSc_analytics/candidate-scores-for-varvis-analysis.xlsx"):
    candidate_table = pd.read_excel(path_to_excel)
    candidate_table["family_id"] = candidate_table["PLIGU number"].apply(lambda x: x.rsplit("-", 1)[0] if str(x)[0] == "L" else None)
    candidate_table.dropna(subset=["family_id"], inplace=True)
    candidate_table.drop(columns=["PLIGU number"], inplace=True)
    candidate_table.rename(columns={"""Zygosity

het /
homo /
hemi /
mosaic""": "zygosity",
                                    """Zygosity
de novo = 2
homo or comphet with multiple affected children = 3
homo with one affected child = 2
comphet with one affected child = 1
X-linked and a boy = 1
X-linked and an affected maternal male relative = 2
other inheritance = 0""": "zygosity_points"}, inplace=True)
    candidate_table.zygosity_points = pd.to_numeric(candidate_table.zygosity_points,
                                                    downcast="unsigned",
                                                    errors="coerce")
    return candidate_table

def filter_family_ids(candidate_table, family_ids):
    candidate_table = candidate_table.loc[candidate_table.family_id.isin(family_ids)]
    candidate_table.reset_index(inplace=True)
    return candidate_table

def clean_variants(candidate_table):
    for i, row in candidate_table.iterrows():
        try:
            candidate_table.loc[i, "variant_1"] = re.findall("^NM.+(?=p\.)", row["variant 1"])[0]
            candidate_table.loc[i, "variant_1"] = re.sub("\s", "", candidate_table.loc[i, "variant_1"])
            if row["variant 2"] != None and not pd.isnull(row["variant 2"]):
                if re.findall("^NM.+(?=p\.)", row["variant 2"]) != []:
                    candidate_table.loc[i, "variant_2"] = re.findall("^NM.+(?=p\.)", row["variant 2"])[0]
        except IndexError:
            pass
    candidate_table.drop(columns=["variant 1", "variant 2"], inplace=True)
    return candidate_table


def interpret_inheritance(candidate_table):
    for i, row in candidate_table.iterrows():
        zygosity = row["zygosity"]
        zygosity_points = row["zygosity_points"]
        if "hemi" in zygosity:
            inheritance = "x_linked"
        elif zygosity[:3] == "het":
            if zygosity_points == 2:
                inheritance = "de_novo"
            elif zygosity_points == 0:
                inheritance = "unknown"
            elif zygosity_points == 1:
                inheritance = "ad_inherited"
        elif "comphet" in zygosity.lower():
            inheritance = "comphet"
        elif zygosity == "homo":
            inheritance = "homo"
        else:
            inheritance = None
        candidate_table.loc[i, "inheritance"] = inheritance
    return candidate_table


def get_vcf_strings(candidate_table):
    for i, row in candidate_table.iterrows():
        instance = make_vep_vcf_string_requests(row["variant_1"], row["Gene"], row['Pos\n\nStart'])
        if instance == None:
            candidate_table.loc[i, "variant_1_vcf"] = None
            candidate_table.loc[i, "variant_1_vcf"] = None
        else:
            candidate_table.loc[i, "variant_1_vcf"] = instance.vcf_string
        if not pd.isnull(row["variant_2"]):
            instance = make_vep_vcf_string_requests(row["variant_2"], row["Gene"], row['Pos\n\nStart.1'])
            if instance == None:
                candidate_table.loc[i, "variant_2_vcf"] = None
            else:
                candidate_table.loc[i, "variant_2_vcf"] = instance.vcf_string
    return candidate_table

def get_candidate_scores(candidate_table):
    for i, row in candidate_table.iterrows():
        print(i)
        instance = make_vep_requests(hgvs_string=row["variant_1"],
                                     gene_symbol=row["Gene"],
                                     start_pos=row['Pos\n\nStart'],
                                     inheritance=row["inheritance"],
                                     hgvs_string_2=row["variant_2"],
                                     start_pos_2=row['Pos\n\nStart'])
        if instance == None:
            candidate_table.loc[i, "variant_1_vcf"] = None
            candidate_table.loc[i, "variant_2_vcf"] = None
            candidate_table.loc[i, "candidate_score"] = None
        else:
            instance.calculate_candidate_score()
            candidate_table.loc[i, "variant_1_vcf"] = instance.vcf_string
            candidate_table.loc[i, "candidate_score"] = instance.candidate_score
            if instance.other_autocasc_obj is not None:
                candidate_table.loc[i, "variant_2_vcf"] = instance.other_autocasc_obj.vcf_string
    return candidate_table

def make_vep_vcf_string_requests(hgvs_string, gene_symbol, start_pos):
    instance = AutoCaSc(hgvs_string, mode="web")
    instance.retrieve_data(gnomad=False)
    if instance != None:
        if instance.status_code == 200:
            return instance
    instance = AutoCaSc(gene_symbol + ":" + hgvs_string.split(":")[1], mode="web")
    instance.retrieve_data(gnomad=False)
    if instance.status_code == 200:
        if str(instance.vcf_string.split(":")[1]) == str(start_pos):
            return instance
    return None

def make_vep_requests(hgvs_string,
                      gene_symbol,
                      start_pos,
                      inheritance,
                      hgvs_string_2=None,
                      start_pos_2=None):
    request_cache = "/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/sonstige/data/"
    if not isinstance(hgvs_string, str):
        return None
    if inheritance != "comphet":
        instance = AutoCaSc(variant=hgvs_string,
                            inheritance=inheritance,
                            path_to_request_cache_dir=request_cache)
        if instance != None:
            if instance.status_code == 200:
                return instance
        instance = AutoCaSc(gene_symbol + ":" + hgvs_string.split(":")[1],
                            inheritance=inheritance,
                            path_to_request_cache_dir=request_cache)
        if instance.status_code == 200:
            if str(instance.vcf_string.split(":")[1]) == str(start_pos):
                return instance
        return None
    else:
        if not isinstance(hgvs_string_2, str):
            return None
        instance_2 = AutoCaSc(variant=hgvs_string_2,
                              inheritance=inheritance,
                              path_to_request_cache_dir=request_cache)
        if instance_2 != None:
            if instance_2.status_code == 200:
                instance_1 = AutoCaSc(variant=hgvs_string,
                                      inheritance=inheritance,
                                      other_autocasc_obj=instance_2,
                                      path_to_request_cache_dir=request_cache)
                if instance_1 != None:
                    if instance_1.status_code == 200:
                        return instance_1
                instance_1 = AutoCaSc(gene_symbol + ":" + hgvs_string.split(":")[1],
                                      inheritance=inheritance,
                                      other_autocasc_obj=instance_2,
                                      path_to_request_cache_dir=request_cache)
                if instance_1.status_code == 200:
                    if str(instance_1.vcf_string.split(":")[1]) == str(start_pos):
                        return instance_1
                return None

        instance_2 = AutoCaSc(gene_symbol + ":" + hgvs_string_2.split(":")[1])
        if instance_2.status_code == 200:
            if str(instance_2.vcf_string.split(":")[1]) == str(start_pos_2):
                instance_1 = AutoCaSc(variant=hgvs_string,
                                      inheritance=inheritance,
                                      other_autocasc_obj=instance_2,
                                      path_to_request_cache_dir=request_cache)
                if instance_1 != None:
                    if instance_1.status_code == 200:
                        return instance_1
                instance_1 = AutoCaSc(gene_symbol + ":" + hgvs_string.split(":")[1],
                                      inheritance=inheritance,
                                      other_autocasc_obj=instance_2,
                                      path_to_request_cache_dir=request_cache)
                if instance_1.status_code == 200:
                    if str(instance_1.vcf_string.split(":")[1]) == str(start_pos):
                        return instance_1
                return None
        return None




# extracted from TriosReal-List excel
family_ids = ['L20-1369', 'L20-0085', 'L20-0325', 'L20-0435', 'L20-1963', 'L19-2165', 'L20-1044', 'L20-0660', 'L20-0553', 'L20-0348', 'L20-0359', 'L20-0505', 'L19-0289', 'L20-0538', 'L20-0231', 'L20-0240', 'L20-0263', 'L20-0309', 'L20-0548', 'L19-0699', 'L16-0467', 'L20-0700', 'L20-0132', 'L20-0213', 'L20-0225', 'L19-1273', 'L20-1275', 'L19-2249', 'L19-2250', 'L19-2201', 'L20-0163', 'L19-2203', 'L19-2271', 'L20-0160', 'L19-2340', 'L20-0068', 'L20-0488', 'L19-1222', 'L19-2020', 'L19-2529', 'L19-2281', 'L19-2453', 'L19-2235', 'L19-2257', 'L19-2372', 'L17-0697', 'L19-2534', 'L19-2097', 'L20-0486', 'L19-2373', 'L20-0373', 'L19-2170', 'L19-2145', 'L20-0705', 'L19-2112', 'L20-0563', 'L17-1132', 'L19-2108', 'L20-0337', 'L19-1053', 'L19-1785', 'L20-0711', 'L20-0677', 'L19-2306', 'L18-1122', 'L19-1732', 'L18-1803', 'L18-1669', 'L20-0238', 'L20-0227', 'L19-1081', 'L19-1784', 'L19-0961', 'L19-2313', 'L19-1800', 'L19-1255', 'L19-0016', 'L19-2521', 'L19-1274', 'L19-0798', 'L18-0207', 'L19-1216', 'L19-1301', 'L19-1166', 'L19-2320', 'L19-2304', 'L19-0936', 'L19-0763', 'L19-1830', 'L19-1179', 'L19-0963', 'L19-0580', 'L19-0011', 'L19-0461']

candidate_table = load_manual_candidate_table()
#candidate_table = filter_family_ids(candidate_table, family_ids)
candidate_table = interpret_inheritance(candidate_table)
candidate_table = clean_variants(candidate_table)
# candidate_table = get_vcf_strings(candidate_table)
update_request_cache("/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/sonstige/data/")
candidate_table = get_candidate_scores(candidate_table)
update_request_cache("/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/sonstige/data/")
#candidate_table = candidate_table.drop(columns=["index"])
candidate_table.columns = [x.replace("\n", "") for x in candidate_table.columns]


candidate_table.to_csv("/Users/johannkaspar/Documents/Promotion/AutoCaSc_analytics/data/candidate_scores_cases_filtered.csv",
                       index=False)
