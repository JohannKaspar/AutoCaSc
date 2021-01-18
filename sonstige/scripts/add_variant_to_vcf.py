import pandas as pd
import click
import shutil
import re
import subprocess, shlex

def get_genotype_alleledepth(inheritance,
                             father_affected=False,
                             mother_affected=False,
                             var_num=0):
    index_gt_ad = ("0/0", "99,0")
    father_gt_ad = ("0/0", "99,0")
    mother_gt_ad = ("0/0", "99,0")
    if inheritance == "de_novo":
        index_gt_ad = ("0/1", "50,49")
    elif inheritance == "homo":
        index_gt_ad = ("1/1", "0,99")
        father_gt_ad = ("0/1", "50,49")
        mother_gt_ad = ("1/0", "50,49")
    elif inheritance == "ad_inherited":
        index_gt_ad = ("0/1", "50,49")
        if father_affected:
            father_gt_ad = ("0/1", "50,49")
        if mother_affected:
            mother_gt_ad = ("0/1", "50,49")
        else:
            raise IOError("Either father or mother has to be affected if variant has been inherited!")
    elif inheritance == "comphet":
        if var_num == 0:
            index_gt_ad = ("0/1", "50,49")
            father_gt_ad = ("0/1", "50,49")
        else:
            index_gt_ad = ("1/0", "50,49")
            mother_gt_ad = ("1/0", "50,49")
    elif inheritance == "x_linked":
        index_gt_ad = ("1/1", "0,99")
        mother_gt_ad = ("1/0", "50, 49")
        father_gt_ad = ("0/0", "99,0")

    return index_gt_ad, father_gt_ad, mother_gt_ad


def convert_variant(input_variant, inheritance, family_order, var_num=0):
    input_variant = re.sub(r"^[\W]", "", input_variant)
    input_variant = input_variant.strip("Chr").strip("chr")
    input_variant = re.sub(r"[^A-Za-z0-9](?!$)", ":", input_variant)

    CHROM, POS, REF, ALT = input_variant.split(":")
    ID = "."
    QUAL = "999"
    FILTER = "PASS"
    INFO = "."  # todo gene symbol == XYZ oder sowas für comphet
    FORMAT = ":".join(["GT:AD:DP:GQ:PL"])

    GQ = "99"
    DP = "99"
    PL = ".,.,."
    index_gt_ad, father_gt_ad, mother_gt_ad = get_genotype_alleledepth(inheritance,
                                                                       var_num=var_num)
    GT_index, AD_index = index_gt_ad
    GT_father, AD_father = father_gt_ad
    GT_mother, AD_mother = mother_gt_ad
    SAMPLE_INFO_index = ":".join([GT_index, AD_index, DP, GQ, PL])
    SAMPLE_INFO_father = ":".join([GT_father, AD_father, DP, GQ, PL])
    SAMPLE_INFO_mother = ":".join([GT_mother, AD_mother, DP, GQ, PL])

    sample_info_dict = {}
    sample_info_dict["index"] = SAMPLE_INFO_index
    sample_info_dict["father"] = SAMPLE_INFO_father
    sample_info_dict["mother"] = SAMPLE_INFO_mother

    sample_info_ordered = [sample_info_dict.get(_id) for _id in family_order]

    record = "\t".join([CHROM, POS, ID, REF, ALT,
                        QUAL, FILTER,
                        INFO, FORMAT] + sample_info_ordered)

    return record


def is_parent(ped_df, sample_id):
    if ped_df.loc[ped_df.individual_id == sample_id, "paternal_id"] not in ped_df.individual_id.to_list() \
            and ped_df.loc[ped_df.individual_id == sample_id, "maternal_id"] not in ped_df.individual_id.to_list():
        return True
    else:
        return False


def get_family_order(ped_file, header):
    pedigree = pd.read_csv(ped_file, sep="\t", header=None)
    pedigree.columns = ["family_id", "individual_id", "paternal_id", "maternal_id", "sex", "affected_status"]

    sample_order = list(header[:-1].split("\t")[:-3])

    who_is_who = {}

    for i, row in pedigree.iterrows():
        if not is_parent(pedigree, row["sample_id"]):
            who_is_who[row["individual_id"]] = "index"
            who_is_who[row["paternal_id"]] = "father"
            who_is_who[row["maternal_id"]] = "mother"

    family_order = [who_is_who.get(_id) for _id in sample_order]

    return family_order


def create_single_var_vcf(header, record, output_path, other_record=None):
    with open(output_path, "wb") as vcf_file:
        for line in header:
            vcf_file.write(line)

            print(line)

        vcf_file.write(record)
        if other_record:
            vcf_file.write(other_record)


def merge_vcfs(input_path, output_path):
    with open(output_path, "wb") as output_file:
        merge_proc = subprocess.run(shlex.split(f'vcftools merge {input_path} {output_path}'),
                       universal_newlines=True,
                       stdout=output_file)

    if merge_proc.returncode == 0:
        subprocess.run(shlex.split(f'vcftools sort {output_path}'))

    else:
        raise IOError


def get_header(input_path):
    proc = subprocess.run(shlex.split(f"zcat {input_path} | head -n 1000 | grep '^#'"),
                   universal_newlines=True,
                   stdout=subprocess.PIPE)
    if proc.returncode == 0:
        return proc.stdout.read()
    else:
        raise IOError


@click.option('--variant', '-v',
              required=True,
              help="Variant(s) to insert.")
@click.option('--other_variant', '-ov',
              required=False,
              help="Second variant in case of comphet.")
@click.option('--inheritance', '-ih',
              required=True,
              help="The way in which the variant has been inherited.")
@click.option('--input_vcf', '-i',
              required=True,
              type=click.Path(exists=True),
              help="Path to original VCF file.")
@click.option('--input_ped', '-p',
              required=True,
              type=click.Path(exists=True),
              help="Path to ped file.")
@click.option('--output_vcf', '-o',
              required=True,
              help="Path to output file.")
def main(variant, other_variant, inheritance, input_vcf, input_ped, output_vcf):
    header = get_header(input_vcf)

    family_order = get_family_order(input_ped, header)
    record = convert_variant(variant, inheritance, family_order)

    if other_variant:
        other_record = convert_variant(other_variant, inheritance, family_order)
        create_single_var_vcf(header, record, output_vcf, other_record=other_record)
    else:
        create_single_var_vcf(header, record, output_vcf)

    merge_vcfs(input_vcf, output_vcf)

    print(f"Variant has been inserted and VCF stored as {output_vcf}")

if __name__ == "__main__":
    main(obj={})