import os
import shutil

import pandas as pd
import shlex
import subprocess

variant_excel = pd.read_excel(
    "/home/johann/PycharmProjects/AutoCaSc_project_folder/sonstige/data/published-variants-for-simulation_2021-01-22.xlsx",
    sheet_name="variants")
variant_df = variant_excel.loc[variant_excel.Use == "yes"][
    ["Case_Number", "Variant_hg19_VCF", "Sex_Index", "GT_Index", "GT_Father", "GT_Mother", "Segregation_Siblings"]]

def create_vcfs_with_inserted_variants():
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
    for case in variant_df.Case_Number.unique():
        for requests_file in ["gnomad_requests_GRCh37",
                              "gnomad_requests_GRCh38",
                              "vep_requests_GRCh37",
                              "vep_requests_GRCh38"]:
            try:
                shutil.copy("/home/johann/PycharmProjects/AutoCaSc_project_folder/sonstige/data/" + requests_file,
                            "/home/johann/PycharmProjects/AutoCaSc_project_folder/sonstige/data/" + requests_file + ".bak")
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

        for trio in ["CEU", "ASH"]:
            subprocess.run(shlex.split("python /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/AutoCaSc_vcf.py "
                                       f"score_vcf -v /home/johann/VCFs/modified_VCFs/annotated/{trio}_{case}.vcf.gz "
                                       f"-p /home/johann/PEDs/{trio}_a{ped_suffix}.ped "
                                       f"-g /home/johann/tools/slivar/gnotate/gnomad.hg37.zip "
                                       f"-j /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/slivar-functions.js "
                                       f"-o /home/johann/trio_scoring_results/{trio}_{case}.csv "
                                       f"-a GRCh37 "
                                       f"-s /home/johann/tools/slivar/slivar "
                                  ))

def score_original_trios():
    for entry in os.scandir("/home/johann/PycharmProjects/AutoCaSc_project_folder/sonstige/data/ped_files"):
        if entry.is_file():
            print(f"working on family {entry.name}")
            subprocess.run(shlex.split(
                "python /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/AutoCaSc_vcf.py "
                f"score_vcf -v /home/johann/VCFs/vcf_filtered_temp_split "
                f"-p {entry.path} "
                f"-g /home/johann/tools/slivar/gnotate/gnomad.hg37.zip "
                f"-j /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/slivar-functions.js "
                f"-o /home/johann/trio_scoring_results/varvis_trios/{entry.name}.csv "
                f"-a GRCh37 "
                f"-s /home/johann/tools/slivar/slivar "
                f"-dbed "
                #f"-ssli"
                ))

score_original_trios()
