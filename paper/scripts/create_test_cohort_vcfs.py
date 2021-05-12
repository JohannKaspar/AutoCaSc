import os
import shutil

import pandas as pd
import shlex
import subprocess
from io import StringIO

def create_vcfs_with_inserted_variants():
    variant_excel = pd.read_excel(
        "/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/published-variants-for-simulation_2021-01-22.xlsx",
        sheet_name="variants")
    variant_df = variant_excel.loc[variant_excel.Use == "yes"][
        ["Case_Number", "Variant_hg19_VCF", "Sex_Index", "GT_Index", "GT_Father", "GT_Mother", "Segregation_Siblings"]]
    for case in variant_df.Case_Number.unique():
        case_df = variant_df.loc[variant_df.Case_Number == case].reset_index(drop=True)
        for trio, trio_path in [("CEU", "/home/johann/VCFs/CeuTrio.hc-joint.MergeVcf.recalibrated.split.vcf.gz"),
                                ("ASH", "/home/johann/VCFs/AshTrio.hc-joint.MergeVcf.recalibrated.split.vcf.gz")
                                ]:
            command = f"python add_variant_to_vcf.py " \
                      f"-v {case_df.Variant_hg19_VCF.values[0]} " \
                      f"--gt_index {case_df.GT_Index.values[0]} " \
                      f"--gt_father {case_df.GT_Father.values[0]} " \
                      f"--gt_mother {case_df.GT_Mother.values[0]} " \
                      f"-i {trio_path} " \
                      f"-o /home/johann/VCFs/modified_VCFs/{trio}_{case}.vcf.gz " \
                      f"-p /home/johann/PEDs/{trio}_h.ped "
            if len(case_df) > 1:
                command += f"-ov {case_df.Variant_hg19_VCF.values[1]} " \
                           f"--gt_index_other_variant {case_df.GT_Index.values[1]} " \
                           f"--gt_father_other_variant {case_df.GT_Father.values[1]} " \
                           f"--gt_mother_other_variant {case_df.GT_Mother.values[1]} "
            subprocess.run(shlex.split(command))

def score_modified_vcfs():
    variant_excel = pd.read_excel(
        "/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/published-variants-for-simulation_2021-01-22.xlsx",
        sheet_name="variants")
    variant_df = variant_excel.loc[variant_excel.Use == "yes"][
        ["Case_Number", "Variant_hg19_VCF", "Sex_Index", "GT_Index", "GT_Father", "GT_Mother", "Segregation_Siblings"]]
    for case in variant_df.Case_Number.unique():
        for requests_file in ["gnomad_requests_GRCh37",
                              "gnomad_requests_GRCh38",
                              "vep_requests_GRCh37",
                              "vep_requests_GRCh38"]:
            try:
                shutil.copy("/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/" + requests_file,
                            "/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/" + requests_file + ".bak")
            except FileNotFoundError:
                pass

        case_df = variant_df.loc[variant_df.Case_Number == case].reset_index(drop=True)
        if len(case_df) == 2:
            # then it's comphet
            ped_suffix = ""
        elif "1" in case_df["GT_Mother"].values[0]:
            ped_suffix = "_mother_affected"
        elif "1" in case_df["GT_Father"].values[0]:
            ped_suffix = "_father_affected"
        else:
            # has to be denovo --> parents not affected
            ped_suffix = ""

        for trio in ["ASH"]:
            os.makedirs(f'/home/johann/trio_scoring_results/synthetic_trios/{date}/cache/{trio}_{case}/',
                        exist_ok=True)
            subprocess.run(shlex.split("python /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/vcfAutoCaSc.py "
                                       f"score_vcf "
                                       # f"-v /mnt/raid/users/johann/VCFs/modified_VCFs/annotated/{trio}_{case}.vcf.gz "
                                       f"-vcf_non_ch '/home/johann/trio_scoring_results/synthetic_trios/2021-04-05/cache/{trio}_{case}/{trio}_a{ped_suffix}_non_comphets' "
                                       f"-vcf_ch '/home/johann/trio_scoring_results/synthetic_trios/2021-04-05/cache/{trio}_{case}/{trio}_a{ped_suffix}_comphets' "
                                       f"-p /home/johann/PEDs/{trio}_a{ped_suffix}.ped "
                                       f"-g /home/johann/tools/slivar/gnotate/gnomad.hg37.zip "
                                       f"-j /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/slivar-functions.js "
                                       f"-o /home/johann/trio_scoring_results/synthetic_trios/{date}/{trio}_{case}.csv "
                                       f"-a GRCh37 "
                                       f"-s /home/johann/tools/slivar/slivar "
                                       f"-blp '/home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/gene_blacklist.txt' "
                                       f"-omim '/home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/OMIM_morbidmap.tsv' "
                                       f"-sys_prim '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/sysid_primary_20210203.csv' "
                                       f"-sys_cand '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/sysid_candidates_20210203.csv' "
                                       #f"-q 500 "
                                       f"-ssli "
                                       f"-dbed "
                                       f"-req_cache '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/' "
                                       f"--cache '/home/johann/trio_scoring_results/synthetic_trios/{date}/cache/{trio}_{case}/' "
                                       #f"-pass "
                                       f"-dp 20 "
                                       f"-ab 0.3 "
                                  ))

def score_original_trios():
    for entry in os.scandir("/mnt/raid/users/johann/PEDs/varvis"):
        if entry.is_file():
            print(f"working on family {entry.name}")
            os.makedirs(f"/home/johann/trio_scoring_results/varvis_trios/{date}/cache/{entry.name.strip('.ped')}/",
                      exist_ok=True)
            subprocess.run(shlex.split(
                "python /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/vcfAutoCaSc.py "
                f"score_vcf -v /mnt/raid/users/johann/VCFs/AutoCaScValidationCohort.ann.vcf.gz.bed_filtered.AC_filtered.impact_filtered.vcf.gz "
                f"-p {entry.path} "
                f"-g /home/johann/tools/slivar/gnotate/gnomad.hg37.zip "
                f"-j /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/slivar-functions.js "
                f"-o /home/johann/trio_scoring_results/varvis_trios/{date}/{entry.name}.csv "
                f"-a GRCh37 "
                f"-s /home/johann/tools/slivar/slivar "
                f"-blp '/home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/gene_blacklist.txt' "
                f"-omim '/home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/OMIM_morbidmap.tsv' "
                f"-sys_prim '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/sysid_primary_20210203.csv' "
                f"-sys_cand '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/sysid_candidates_20210203.csv' "
                f"-dbed "
                #f"-ssli "
                f"-req_cache '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/' "
                f"--cache '/home/johann/trio_scoring_results/varvis_trios/{date}/cache/{entry.name.strip('.ped')}/' "
                #f"-pass "
                f"-dp 20 "
                f"-ab 0.3 "
            ))

def concat_results(path="/home/johann/trio_scoring_results/varvis_trios/new_script_test/"):
    concat_df = pd.DataFrame()
    for entry in os.scandir(path):
        if ".csv" in entry.name and not "lock" in entry.name:
            _df = pd.read_csv(entry.path, sep="\t")
            _df["family_id"] = entry.name.split(".")[0]
            concat_df = pd.concat([concat_df, _df],
                                  ignore_index=True)
    concat_df.loc[concat_df.autocasc_filter != "PASS", "autocasc_filter"] = "FAIL"
    concat_df.to_csv(path + "varvis_trio_results_2021-03-04.csv",
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
        "python /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/vcfAutoCaSc.py "
        f"score_vcf -vcf_non_ch '/home/johann/VCFs/clinvar/clinvar.vcf' "
        f"--cache '/home/johann/VCFs/clinvar/tmp/' "
        f"-o /home/johann/trio_scoring_results/clinvar/{date}/clinvar.csv "
        f"-a GRCh37 "
        f"-dbed "
        f"-ssli "
        f"-req_cache '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/' "
        f"--trio_name clinvar"
    ))


date = "2021-04-09"

score_modified_vcfs()
# score_original_trios()
# score_clinvar()

# score_modified_vcfs()
# concat_results()