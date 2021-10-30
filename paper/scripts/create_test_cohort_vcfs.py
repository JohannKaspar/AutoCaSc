import os
import shutil
import time

import pandas as pd
import shlex
import subprocess

def create_vcfs_with_inserted_variants():
    variant_excel = pd.read_excel(
        "/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/published-variants-for-simulation_2021-01-22.xlsx",
        sheet_name="variants")
    variant_df = variant_excel.loc[variant_excel.Use == "yes"][
        ["Case_Number", "Variant_hg19_VCF", "Sex_Index", "GT_Index", "GT_Father", "GT_Mother", "Segregation_Siblings"]]
    for case in variant_df.Case_Number.unique():
        if case != "sim61":
            continue
        case_df = variant_df.loc[variant_df.Case_Number == case].reset_index(drop=True)
        for trio, trio_path in [("CEU", "/mnt/raid/users/johann/VCFs/CeuTrio.hc-joint.MergeVcf.recalibrated.split.vcf.gz"),
                                ("ASH", "/mnt/raid/users/johann/VCFs/AshTrio.hc-joint.MergeVcf.recalibrated.split.vcf.gz")
                                ]:
            command = f"python add_variant_to_vcf.py " \
                      f"-v {case_df.Variant_hg19_VCF.values[0]} " \
                      f"--gt_index {case_df.GT_Index.values[0]} " \
                      f"--gt_father {case_df.GT_Father.values[0]} " \
                      f"--gt_mother {case_df.GT_Mother.values[0]} " \
                      f"-i {trio_path} " \
                      f"-o /mnt/raid/users/johann/VCFs/modified_VCFs/_new_try_sim40/{trio}_{case}.vcf.gz " \
                      f"-p /home/johann/PEDs/{trio}_h.ped "
            if len(case_df) > 1:
                command += f"-ov {case_df.Variant_hg19_VCF.values[1]} " \
                           f"--gt_index_other_variant {case_df.GT_Index.values[1]} " \
                           f"--gt_father_other_variant {case_df.GT_Father.values[1]} " \
                           f"--gt_mother_other_variant {case_df.GT_Mother.values[1]} "
            subprocess.run(shlex.split(command))

def score_modified_vcfs():
    variant_excel = pd.read_excel(
        "/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/published-variants-for-simulation_2021-10-28-mod.xlsx",
        sheet_name="variants")
    variant_df = variant_excel.loc[variant_excel.Use == "yes"][
        ["Case_Number", "Publication_Name", "Variant_hg19_VCF", "Sex_Index", "GT_Index", "GT_Father", "GT_Mother", "Segregation_Siblings"]]

    inheritance_excel = pd.read_excel(
        "/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/published-variants-for-simulation_2021-01-22.xlsx",
        sheet_name="publications")
    inheritance_excel = inheritance_excel.loc[:, ["Publication_Name", "InheritancePattern_reported"]]

    for case in variant_df.Case_Number.unique():

        case_df = variant_df.loc[variant_df.Case_Number == case].reset_index(drop=True)

        publication_name = case_df["Publication_Name"].values[0]
        reported_inheritance = inheritance_excel.loc[inheritance_excel.Publication_Name == publication_name,
                                                     "InheritancePattern_reported"].values[0]
        index_sex = case_df["Sex_Index"].values[0]
        if index_sex not in ["male", "female"]:
            index_sex = "male"

        if "dominant" in reported_inheritance:  # there are no cases where two parents are affected
            if "1" in case_df["GT_Mother"].values[0]:
                ped_suffix = "_mother_affected"
            elif "1" in case_df["GT_Father"].values[0]:
                ped_suffix = "_father_affected"
            else:
                # has to be denovo --> parents not affected
                ped_suffix = ""
        elif ("recessive" in reported_inheritance) or ("X-linked" in reported_inheritance):
            if "1/1" in case_df["GT_Mother"].values[0]:
                ped_suffix = "_mother_affected"
            elif "1/1" in case_df["GT_Father"].values[0]:
                ped_suffix = "_father_affected"
            else:
                ped_suffix = ""
            if len(case_df) == 2:  # then it's comphet
                if (case_df["GT_Mother"].values[0][2] == "1") and (case_df["GT_Mother"].values[1][2] == "1"):
                    ped_suffix = "_mother_affected"
                elif (case_df["GT_Father"].values[0][2] == "1") and (case_df["GT_Father"].values[1][2] == "1"):
                    ped_suffix = "_father_affected"
                else:
                    ped_suffix = ""

        for trio in ["ASH", "CEU"]:
            print(f"working on {trio}{case}")
            os.makedirs(f'/home/johann/trio_scoring_results/synthetic_trios/{date}/cache/{trio}_{case}/',
                        exist_ok=True)
            subprocess.run(shlex.split("python /home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_core/vcfAutoCaSc.py "
                                       f"score_vcf "
                                       f"-v /mnt/raid/users/johann/VCFs/modified_VCFs/annotated/concat/{trio}_{case}.vcf.gz "
                                       # f"-vcf_non_ch '/home/johann/trio_scoring_results/synthetic_trios/{date}/cache/{trio}_{case}/{trio}_a_{index_sex}{ped_suffix}_non_comphets' "
                                       # f"-vcf_ch '/home/johann/trio_scoring_results/synthetic_trios/{date}/cache/{trio}_{case}/{trio}_a_{index_sex}{ped_suffix}_comphets' "
                                       f"-p /home/johann/PEDs/{trio}_a_{index_sex}{ped_suffix}.ped "
                                       f"-g /home/johann/tools/slivar/gnotate/gnomad.hg37.zip "
                                       f"-j /home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_core/data/slivar-functions.js "
                                       f"-o /home/johann/trio_scoring_results/synthetic_trios/{date}/{trio}_{case}.csv "
                                       f"-a GRCh37 "
                                       f"-s /home/johann/tools/slivar/slivar "
                                       f"-blp '/home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_core/data/gene_blacklist.txt' "
                                       f"-omim '/home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_core/data/OMIM_morbidmap.tsv' "
                                       f"-sys_prim '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/sysid_old/sysid_primary_2020_09_24.csv' "
                                       f"-sys_cand '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/sysid_old/sysid_candidates_2020_09_24.csv' "
                                       # f"-ssli "
                                       # f"-dbed "
                                       f"-req_cache '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/' "
                                       f"--cache '/home/johann/trio_scoring_results/synthetic_trios/{date}/cache/{trio}_{case}' "
                                       f"-dp 20 "
                                       f"-ab 0.3 "
                                  ))

def score_original_trios():
    for entry in os.scandir("/mnt/raid/users/johann/PEDs/varvis"):
        if entry.is_file():
            print(f"working on family {entry.name}")

            # if (entry.name + ".csv") in os.listdir("/home/johann/trio_scoring_results/varvis_trios/2021-10-28/"):
            if entry.name != "L19-2257.ped":
                continue

            os.makedirs(f"/home/johann/trio_scoring_results/varvis_trios/{date}/cache/{entry.name.strip('.ped')}/",
                      exist_ok=True)
            subprocess.run(shlex.split(
                "python /home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_core/vcfAutoCaSc.py "
                f"score_vcf "
                f"-v /mnt/raid/users/johann/VCFs/AutoCaScValidationCohort.ann.vcf.gz.bed_filtered.AC_filtered.impact_filtered.vcf.gz "
                # f"-vcf_non_ch '/home/johann/trio_scoring_results/varvis_trios/{date}/cache/{entry.name.strip('.ped')}/{entry.name.strip('.ped')}_non_comphets' "
                # f"-vcf_ch '/home/johann/trio_scoring_results/varvis_trios/{date}/cache/{entry.name.strip('.ped')}/{entry.name.strip('.ped')}_comphets' "         
                f"-p {entry.path} "
                f"-g /home/johann/tools/slivar/gnotate/gnomad.hg37.zip "
                f"-j /home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_core/data/slivar-functions.js "
                f"-o /home/johann/trio_scoring_results/varvis_trios/{date}/{entry.name}.csv "
                f"-a GRCh37 "
                f"-s /home/johann/tools/slivar/slivar "
                f"-blp '/home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_core/data/gene_blacklist.txt' "
                f"-omim '/home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_core/data/OMIM_morbidmap.tsv' "
                f"-sys_prim '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/sysid_primary_20210203.csv' "
                f"-sys_cand '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/sysid_candidates_20210203.csv' "
                f"-dbed "
                # f"-ssli "
                f"-req_cache '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/' "
                f"--cache '/home/johann/trio_scoring_results/varvis_trios/{date}/cache/{entry.name.strip('.ped')}/' "
                f"-dp 20 "
                f"-ab 0.3 "
            ))

def concat_results(path):
    concat_df = pd.DataFrame()
    for entry in os.scandir(path):
        if ".csv" in entry.name and not "lock" in entry.name:
            _df = pd.read_csv(entry.path, sep="\t")
            _df["family_id"] = entry.name.split(".")[0]
            concat_df = pd.concat([concat_df, _df],
                                  ignore_index=True)
    concat_df.loc[concat_df.autocasc_filter != "PASS", "autocasc_filter"] = "FAIL"
    concat_df.to_csv(path + f"varvis_trio_results_{date}.csv",
                     index=False)


# not used but useful
def get_seq_difference(ref, alt):
    n_same_start = 0
    while ref[:n_same_start+1] == alt[:n_same_start+1]:
        n_same_start += 1
    ref_new = ref[n_same_start:]
    alt_new = alt[n_same_start:]

    n_same_end = 0
    while ref_new[n_same_end-1:] == alt_new[n_same_end-1:]:
        n_same_end -= 1

    if n_same_end == 0:
        ref_new = ref[n_same_start:]
        alt_new = alt[n_same_start:]
    else:
        ref_new = ref_new[:n_same_end]
        alt_new = alt_new[:n_same_end]

    if alt_new == "":
        alt_new = "-"
    if ref_new == "":
        ref_new = "-"

    return f"{ref_new}>{alt_new}"

def score_clinvar():
    subprocess.run(shlex.split(
        "python /home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_core/vcfAutoCaSc.py "
        f"score_vcf -vcf_non_ch '/home/johann/VCFs/clinvar/clinvar.vcf' "
        f"--cache '/home/johann/VCFs/clinvar/tmp/' "
        f"-o /home/johann/trio_scoring_results/clinvar/{date}/clinvar.csv "
        f"-a GRCh37 "
        f"-dbed "
        f"-ssli "
        f"-req_cache '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/' "
        f"--trio_name clinvar"
    ))


date = "2021-10-28"

# create_vcfs_with_inserted_variants()
score_modified_vcfs()
# score_original_trios()
# score_modified_vcfs()
# score_clinvar()

# score_modified_vcfs()
# concat_results(f"/home/johann/trio_scoring_results/varvis_trios/{date}/")