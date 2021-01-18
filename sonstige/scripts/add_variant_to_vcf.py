import pandas as pd
import click
import shutil
import re


def copy_vcf(input_path, output_path):
    shutil.copy(input_path, output_path)


def get_genotype_code(inheritance, father_affected=False, mother_affected=False):
    father_gt = "0/0"
    mother_gt = "0/0"
    if inheritance == "de_novo":
        index_gt = "0/1"
    elif inheritance == "homo":
        index_gt = "1/1"
        father_gt = "0/1"
        mother_gt = "1/0"
    elif inheritance in ["ad_inherited", "comphet"]:
        index_gt = "0/1"
        if father_affected:
            father_gt = "0/1"
        if mother_affected:
            mother_gt = "0/1"
        else:
            raise IOError("Either father or mother has to be affected if variant has been inherited!")
    elif inheritance == "x_linked":
        index_gt = "1"
        mother_gt = "1/0"
        father_gt = "0"

    return index_gt, father_gt, mother_gt

def get_allele_depth(inheritance):
    if inheritance == "homo":
        return "99,0"
    else:
        return "50,49"


def convert_variant(input_variant, inheritance):
    input_variant = re.sub(r"^[\W]", "", input_variant)
    input_variant = input_variant.strip("Chr").strip("chr")
    input_variant = re.sub(r"[^A-Za-z0-9](?!$)", ":", input_variant)

    chromosome, position, reference, alternative = input_variant.split(":")
    quality = "999"
    gq = "99"
    dp = "99"
    ad = get_allele_depth(inheritance)
    pl = ".,.,."

    gt_index, gt_father, gt_mother = get_genotype_code(inheritance)

    mod_line

def write_line(mod_line, output_path):
    with open(output_path, "a") as vcf_file:
        vcf_file.write(mod_line)

# todo VCF mit HEader entsprechend original VCF erstellen und nur der einzufügenden Variante, mit VCF tools mergen und sortieren