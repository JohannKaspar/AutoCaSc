import os
import shutil

import pandas as pd
import shlex
import subprocess
from io import StringIO

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
                                       f"-o /home/johann/trio_scoring_results/modified_trios/{trio}_{case}.csv "
                                       f"-a GRCh37 "
                                       f"-s /home/johann/tools/slivar/slivar "
                                       f"-q 200"
                                  ))

def score_original_trios():
    for entry in os.scandir("/home/johann/PycharmProjects/AutoCaSc_project_folder/sonstige/data/ped_files"):
        if entry.is_file():
            print(f"working on family {entry.name}")
            subprocess.run(shlex.split(
                "python /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/AutoCaSc_vcf.py "
                f"score_vcf -v /home/johann/VCFs/AutoCaScValidationCohort.ann.vcf.gz.bed_filtered.AC_filtered.impact_filtered.vcf.gz "
                f"-p {entry.path} "
                f"-g /home/johann/tools/slivar/gnotate/gnomad.hg37.zip "
                f"-j /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/slivar-functions.js "
                f"-o /home/johann/trio_scoring_results/varvis_trios/{entry.name}.csv "
                f"-a GRCh37 "
                f"-s /home/johann/tools/slivar/slivar "
                f"-dbed "
                f"-ssli"
            ))

    """for entry in os.scandir("/home/johann/trio_scoring_results/varvis_trios/rescore_2/"):
        if entry.is_file():
            family_id = entry.name.split(".")[0]
            print(f"working on family {entry.name}")
            subprocess.run(shlex.split(
                "python /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/AutoCaSc_vcf.py "
                f"score_vcf -v /home/johann/VCFs/AutoCaScValidationCohort.ann.vcf.gz.bed_filtered "
                f"-p /home/johann/PycharmProjects/AutoCaSc_project_folder/sonstige/data/ped_files/{family_id}.ped "
                f"-g /home/johann/tools/slivar/gnotate/gnomad.hg37.zip "
                f"-j /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/slivar-functions.js "
                f"-o /home/johann/trio_scoring_results/varvis_trios/{entry.name} "
                f"-a GRCh37 "
                f"-s /home/johann/tools/slivar/slivar "
                f"-dbed "
                #f"-ssli"
               ))"""


def get_ids(family_id, ped_directory="/home/johann/PycharmProjects/AutoCaSc_project_folder/sonstige/data/ped_files/"):
    ids = pd.read_csv(f"{ped_directory}{family_id}.ped",
                header=None,
                usecols=[1],
                sep="\t")[1].to_list()
    return ids

def get_seq_difference(ref, alt):
    """if val == "+":
        if ref in alt:
            try:
                x_1, x_2 = alt.split(ref)
            except:
                ValueError
            if x_1 == "":
                return val + x_2
            else:
                return val + x_1

    if val == "-":
        if alt in ref:
            if alt == "-":
                return val + ref
            try:
                x_1,  x_2 = ref.split(alt)
            except ValueError:
                return()
            if x_1 == "":
                return val + x_2
            else:
                return val + x_1"""

    # if substitution of different length
    #if val != "/":
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

    #n_same_end += 1




def concat_vcf(directory, family_id):
    files = [f"{directory}{family_id}_comphets", f"{directory}{family_id}_non_comphets"]
    for file in files:
        subprocess.run(shlex.split(f"tabix -p vcf {file}"))
    merge_proc = subprocess.run(shlex.split(f'bcftools concat -a -o {directory}.vcf -O v {files[0]} {files[1]}'),
                                universal_newlines=True)

    subprocess.run(shlex.split(f'bcftools sort -o {directory}.vcf -O z {directory}.vcf'))

def lookup_parameters():
    for entry in os.scandir("/home/johann/VCFs/slivar_filtered/"):
        if not entry.is_file():
            family_id = entry.name
            print(f"working on {family_id}")
            ids = get_ids(family_id)

            mim_mapped = pd.read_csv(f"/home/johann/trio_scoring_results/varvis_trios/mim_mapped/{entry.name}")
            variants = pd.DataFrame()
            mim_mapped = mim_mapped.loc[mim_mapped.variant.str[0].isin(["X", "Y"] + [str(x) for x in range(10)])]
            variants["CHROM"] = mim_mapped.variant.apply(lambda x: x.split(":")[0])
            variants["POS"] = mim_mapped.variant.apply(lambda x: x.split(":")[1])
            variants.to_csv("positions", sep="\t", index=False)

            variants["REF"] = mim_mapped.variant.apply(lambda x: x.split(":")[2])
            variants["ALT"] = mim_mapped.variant.apply(lambda x: x.split(":")[3])

            concat_vcf(f"/home/johann/VCFs/slivar_filtered/{family_id}")
            basic_command = f"vcftools --vcf /home/johann/VCFs/AutoCaScValidationCohort.ann.vcf.gz.bed_filtered " \
                            f"--indv {ids[0]} " \
                            f"--indv {ids[1]} " \
                            f"--indv {ids[2]} " \
                            f"--positions positions "

            ac_process = subprocess.run(basic_command + f"--get-INFO AC --stdout",
                                        shell=True,
                                        universal_newlines=True,
                                        stdout=subprocess.PIPE)
            ac_df = pd.read_csv(StringIO(ac_process.stdout), sep="\t")
            ac_df = ac_df.astype(str).reset_index()
            ac_df["ALT"] = ac_df.ALT.str.replace("*", "-")

            for i, row in variants.iterrows():
                variants.loc[i, "seq_diff"] = get_seq_difference(row["REF"], row["ALT"])
            for i, row in ac_df.iterrows():
                ac_df.loc[i, "seq_diff"] = get_seq_difference(row["REF"], row["ALT"])

            variants = variants.merge(ac_df[["CHROM", "POS", "seq_diff", "index"]], on=["CHROM", "POS", "seq_diff"],
                                      how="left")

            index_to_use = variants["index"].to_list()

            gq_process = subprocess.run(basic_command + f"--extract-FORMAT-info GQ --stdout",
                shell=True,
                universal_newlines=True,
                stdout=subprocess.PIPE)
            gq_df = pd.read_csv(StringIO(gq_process.stdout), sep="\t")
            ad_process = subprocess.run(basic_command + f"--extract-FORMAT-info AD --stdout",
                                        shell=True,
                                        universal_newlines=True,
                                        stdout=subprocess.PIPE)
            ad_df = pd.read_csv(StringIO(ad_process.stdout), sep="\t")
            dp_process = subprocess.run(basic_command + f"--extract-FORMAT-info DP --stdout",
                                        shell=True,
                                        universal_newlines=True,
                                        stdout=subprocess.PIPE)
            dp_df = pd.read_csv(StringIO(dp_process.stdout), sep="\t")
            """
            qual_process = subprocess.run(basic_command + f"--site-quality --stdout",
                                          shell=True,
                                          universal_newlines=True,
                                          stdout=subprocess.PIPE)
            qual_df = pd.read_csv(StringIO(qual_process.stdout), sep="\t")"""

            for df in [gq_df, ad_df, dp_df]:
                df = df.astype(str)
                df = df.loc[df.index.isin(index_to_use)]
                variants = variants.merge(df, on=["CHROM", "POS"],
                                          how="left")



            param_mapped = mim_mapped.concat(variants, axis=0)  #todo check axis, drop CHROM POS ALT, check if fits
            print("stop")

score_original_trios()
# score_modified_vcfs()
