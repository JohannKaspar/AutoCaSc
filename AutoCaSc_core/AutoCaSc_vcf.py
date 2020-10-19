import os
import shlex
import subprocess
from statistics import mean
import click
import time
from AutoCaSc_core import AutoCaSc
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import re
import shutil

def create_bed_file(assembly="GRCh37", ensembl_version="101"):
    print("create bed currently not working")
    """This function expands the exon positions by 50 to each side. Then one has to use the UNIX commands:
    sort -k 1,1 -k2,2n Homo_sapiens.GRCh38.bed | bedtools merge -i - | awk '{print "chr"$0}' - > GRCh38.bed"""
    df = pd.read_csv(f"/home/johann/AutoCaSc_core/data/BED/Homo_sapiens.{assembly}.{ensembl_version}.gtf", sep="\t", skiprows=5)
    df.columns = ["chr", "source", "type", "start", "stop", "unknown_1", "unknown_2", "unknown_3", "info"]
    exons = df.loc[df.type == "exon"]
    exons = exons.loc[exons["info"].str.contains("protein_coding")]
    exons = exons[["chr", "start", "stop"]]
    exons.start = exons.start - 50
    exons.stop = exons.stop + 50
    "chr" + exons.chr.astype(str)
    exons.to_csv(f"/home/johann/AutoCaSc_core/data/BED/Homo_sapiens.{assembly}_unsorted.bed", sep="\t", index=False, header=False)

    # cache = "/home/johann/AutoCaSc_core/data/BED/temp/"
    # if not os.path.exists(cache):
    #     os.mkdir(cache)
    # exons.to_csv(f"{cache}{assembly}_unsorted.bed", sep="\t", index=False, header=False)
    # subprocess.run(shlex.split(f"echo hi"),
    #                stdout=subprocess.PIPE,
    #                universal_newlines=True,
    #                shell=True)
    # subprocess.run(shlex.split(f"sort -k 1,1 -k2,2n {cache}{assembly}_unsorted.bed"),
    #                stdout=subprocess.PIPE,
    #                universal_newlines=True)
    # subprocess.Popen(shlex.split(f"bedtools merge -i - > data/BED/{assembly}.bed"),
    #                stdout=subprocess.PIPE,
    #                universal_newlines=True,
    #                shell=True)
    # shutil.rmtree(cache)

def thread_function_AutoCaSc_classic(param_tuple):
    row, assembly = param_tuple
    chrom = row["#CHROM"]
    pos = row["POS"]
    ref = row["REF"]
    alt = row["ALT"]
    alt = alt.replace("*", "-")
    variant_vcf = ":".join(map(str, [chrom, pos, ref, alt]))
    if "de_novo" in row["INFO"]:
        inheritance = "de_novo"
    elif "homo" in row["INFO"]:
        inheritance = "homo"
    elif "x_linked" in row["INFO"]:
        inheritance = "x_linked"
    else:
        inheritance = "other"
    AutoCaSc_instance = AutoCaSc(variant_vcf, inheritance=inheritance, assembly=assembly)
    if AutoCaSc_instance.status_code in [503, 497, 496, 201]:
        print("There has been an issue with a variant. Retrying...")
        for i in range(10):
            if AutoCaSc_instance.status_code in [503, 497, 496, 201]:
                time.sleep(3)
                AutoCaSc_instance = AutoCaSc(variant_vcf, inheritance=inheritance, assembly=assembly)
            else:
                break
    return variant_vcf, AutoCaSc_instance

def thread_function_comphet_classic(param_tuple):
    """function for bases line annotating cmpound heterozygous variants
    """
    row, assembly = param_tuple
    chrom = row["#CHROM"]
    pos = row["POS"]
    ref = row["REF"]
    alt = row["ALT"]
    alt = alt.replace("*", "-")
    variant_vcf = ":".join(map(str, [chrom, pos, ref, alt]))
    inheritance = "comphet"
    AutoCaSc_instance = AutoCaSc(variant_vcf, inheritance=inheritance, assembly="GRCh38")
    if AutoCaSc_instance.status_code in [503, 497, 496, 201]:
        print("There has been an issue with a variant. Retrying...")
        for i in range(10):
            if AutoCaSc_instance.status_code in [503, 497, 496, 201]:
                time.sleep(3)
                AutoCaSc_instance = AutoCaSc(variant_vcf, inheritance=inheritance, assembly="GRCh38")
            else:
                break
    return variant_vcf, AutoCaSc_instance

def score_non_comphets(filtered_vcf, cache, trio_name, assembly):
    # this loads the vcf containing all variants but compound heterozygous ones and converts it to a DataFrame
    with open(filtered_vcf, "r") as inp, open(
            f"{cache}/temp_{trio_name}.tsv",
            "w") as out:
        for row in inp:
            if "##" not in row:
                out.write(row)

    vcf_annotated = pd.read_csv(f"{cache}/temp_{trio_name}.tsv", sep="\t")
    # vcf_annotated = vcf_annotated.loc[vcf_annotated["#CHROM"] == "chr1"]
    vcf_chunks = [vcf_annotated.loc[i, :] for i in range(len(vcf_annotated))]

    # this starts annotation of all variants that are comphet
    with ThreadPoolExecutor() as executor:
        instances_iterator = executor.map(thread_function_AutoCaSc_classic,
                                          zip(vcf_chunks, [assembly] * len(vcf_chunks)))
        # instances_iterator = executor.map(thread_function_AutoCaSc_classic, vcf_chunks)
    result_df = pd.DataFrame()
    i = 0
    for _variant_vcf, _AutoCaSc_instance in instances_iterator:
        result_df.loc[i, "variant"] = _variant_vcf
        try:
            result_df.loc[i, "gene"] = _AutoCaSc_instance.gene_symbol
        except AttributeError:
            result_df.loc[i, "gene"] = "-"
        result_df.loc[i, "hgvsc"] = _AutoCaSc_instance.hgvsc_change
        result_df.loc[i, "hgvsp"] = _AutoCaSc_instance.hgvsp_change
        result_df.loc[i, "candidate_score"] = _AutoCaSc_instance.candidate_score
        result_df.loc[i, "literature_score"] = _AutoCaSc_instance.literature_score
        result_df.loc[i, "CADD_phred"] = _AutoCaSc_instance.__dict__.get("cadd_phred") or 0
        result_df.loc[i, "impact"] = _AutoCaSc_instance.impact
        result_df.loc[i, "inheritance"] = _AutoCaSc_instance.inheritance
        result_df.loc[i, "status_code"] = _AutoCaSc_instance.status_code
        i += 1

    return result_df

def score_comphets(comphets_vcf, cache, trio_name, output_path, assembly):
    #this loads the vcf containing all variants but compound heterozygous ones and converts it to a DataFrame
    with open(comphets_vcf, "r") as inp, open(f"{cache}/temp_{trio_name}.tsv", "w") as out:
        for row in inp:
            if "##" not in row:
                out.write(row)

    comphets_df = pd.read_csv(f"{cache}/temp_{trio_name}.tsv", sep="\t")
    # comphets_df = comphets_df.loc[comphets_df["#CHROM"] == "chr1"]
    comphet_chunks = [comphets_df.loc[i, :] for i in range(len(comphets_df))]
    comphet_cross_df = pd.DataFrame()

    for i in range(len(comphets_df)):
        chrom_1 = comphets_df.loc[i,"#CHROM"]
        pos_1 = comphets_df.loc[i,"POS"]
        ref_1 = comphets_df.loc[i,"REF"]
        alt_1 = comphets_df.loc[i,"ALT"]
        vcf_1 = ":".join(map(str, [chrom_1, pos_1, ref_1, alt_1]))
        info = comphets_df.loc[i,"INFO"]
        slivar_substrings = re.findall(r'(?:[^\/\,]+\/){6}[CTGA\-*]+', info.split("slivar_comphet=")[1])
        for substring in slivar_substrings:
            cross_df_index = len(comphet_cross_df)
            comphet_cross_df.loc[cross_df_index, "var_1"] = vcf_1
            chrom_2, pos_2, ref_2, alt_2 = substring.split("/")[-4:]
            vcf_2 = ":".join([chrom_2, pos_2, ref_2, alt_2])
            comphet_cross_df.loc[cross_df_index, "var_2"] = vcf_2

    with ThreadPoolExecutor() as executor:
        instances_iterator = executor.map(thread_function_comphet_classic,
                                          zip(comphet_chunks, [assembly] * len(comphet_chunks)))
        instances_dict = {vcf: instance for vcf, instance in instances_iterator}

    # loop through the variants and re_rate the impact based on the other corresponding variant
    for i in range(len(comphet_cross_df)):
        var_1 = comphet_cross_df.loc[i, "var_1"]
        var_2 = comphet_cross_df.loc[i, "var_2"]
        instance_1 = instances_dict.get(var_1)
        instance_2 = instances_dict.get(var_2)

        if instance_1.status_code == 200:
            try:
                instance_1.other_impact = instance_2.impact
                instance_1.rate_impact()
            except AttributeError:
                pass

    # loop through the pairs and calculate the mean as variant scores
    for i in range(len(comphet_cross_df)):
        var_1 = comphet_cross_df.loc[i, "var_1"]
        var_2 = comphet_cross_df.loc[i, "var_2"]
        instance_1 = instances_dict.get(var_1)
        instance_2 = instances_dict.get(var_2)
        combined_variant_score = mean([instance_1.variant_score,
                                       instance_2.variant_score])
        instance_1.variant_score = combined_variant_score
        instance_2.variant_score = combined_variant_score

        try:
            comphet_cross_df.loc[i, "gene"] = instance_1.gene_symbol
        except AttributeError:
            comphet_cross_df.loc[i, "gene"] = "-"
        comphet_cross_df.loc[i, "hgvsc"] = instance_1.hgvsc_change
        comphet_cross_df.loc[i, "hgvsp"] = instance_1.hgvsp_change
        comphet_cross_df.loc[i, "candidate_score"] = instance_1.candidate_score
        comphet_cross_df.loc[i, "literature_score"] = instance_1.literature_score
        comphet_cross_df.loc[i, "CADD_phred"] = instance_1.__dict__.get("cadd_phred") or 0
        comphet_cross_df.loc[i, "impact"] = instance_1.impact
        comphet_cross_df.loc[i, "other_impact"] = instance_1.other_impact
        comphet_cross_df.loc[i, "inheritance"] = instance_1.inheritance
        comphet_cross_df.loc[i, "status_code"] = instance_1.status_code

    comphet_cross_df = comphet_cross_df.sort_values(by="candidate_score", ascending=False).reset_index(drop=True)
    comphet_cross_df = comphet_cross_df.drop_duplicates(subset=["var_1"], keep="first")
    comphet_cross_df.to_csv(f"{output_path}_scored_comphets.csv",
                            index=False)

    return comphet_cross_df

@click.group(invoke_without_command=True)  # Allow users to call our app without a command
@click.pass_context
@click.option('--verbose', '-v', is_flag=True, help="Increase output verbosity level")
def main(ctx, verbose):
    group_commands = ['candidates']
    """
            AutoCaSc_core is a command line tool that helps scoring
            variants of NDD patients for their plausible pathogenicity.\n
            example: python AutoCaSc_core.py batch -i input/file/path
    """

    if ctx.invoked_subcommand is None:
        # No command supplied
        # Inform user on the available commands when running the app
        click.echo("Specify one of the commands below")
        print(*group_commands, sep='\n')

    ctx.obj['VERBOSE'] = verbose

@main.command("score_vcf")
@click.option("--vcf_file", "-v",
              required=True,
              type=click.Path(exists=True),
              help="Path to called vcf file.")
@click.option("--ped_file", "-p",
              required=True,
              type=click.Path(exists=True),
              help="Path to pedigree file.")
@click.option("--gnotate_file", "-g",
              required=True,
              type=click.Path(exists=True),
              help="Path to slivar gnotate file.")
@click.option("--javascript_file", "-j",
              required=True,
              type=click.Path(exists=True),
              help="Path to slivar javascript file.")
@click.option("--output_path", "-o",
              help="Output path for table with scored variants.")
@click.option("--denovo_af_max", "-daf",
              default="0.00000001",
              help="Popmax AF in gnomad for denovo variants.")
@click.option("--x_recessive_af_max", "-raf",
              default="0.0001",
              help="Popmax AF in gnomad for x-linked recessive variants.")
@click.option("--compf_af_max", "-raf",
              default="0.0001",
              help="Popmax AF in gnomad for compound heterozygous variants.")
@click.option("--nhomalt", "-n",
              default="0",
              help="Max number of alternative sequence homozygotes in gnomad.")
@click.option("--cache", "-c",
              type=click.Path(exists=True),
              help="Path to cache for temporary data.")
@click.option("--assembly", "-a",
              default="GRCh37",
              help="Reference genome to use, either GRCh37 or GRCh38.")
@click.option("--non_coding", "-n",
              is_flag=True,
              help="Add to deactivate filter for protein coding variants. (-50 and +50bp of exons)")
def score_vcf(vcf_file, ped_file, gnotate_file, javascript_file, output_path,
              denovo_af_max, x_recessive_af_max, compf_af_max, nhomalt, cache, assembly, non_coding):
    trio_name = vcf_file.split("/")[-1]
    #ToDo noch mal genau mit regex überlegen, damit es nicht "family_XYZ.vcf.gz_filtered" heißt
    if not cache:
        cache = vcf_file.rstrip(trio_name) + "temp/"
        if not os.path.exists(cache):
            os.mkdir(cache)
        click.echo(f"Using cache directory: {cache}")
    if not non_coding:
        bed_filter_command = f'bedtools intersect -wa -a {vcf_file} -b /home/johann/AutoCaSc_core/data/BED/Homo_sapiens.{assembly}.bed > {cache}_temp'
        bed_filter_process = subprocess.run(shlex.split(bed_filter_command),
                                        stdout=subprocess.PIPE,
                                        universal_newlines=True)
        if bed_filter_process.returncode == 0:
            vcf_file = f"{cache}vcf_protein_coding"
        else:
            click.echo("There has been a problem filtering for protein coding variants only! Continuing with all variants.")

    if float(denovo_af_max) == 0.0:
        denovo_af_max = "0.00000001"

    slivar_noncomp_command = f'slivar_folder/slivar expr -v {vcf_file} -j {javascript_file} -p {ped_file} ' \
                             f'--pass-only -g {gnotate_file} -o {cache}{trio_name}_filtered ' \
                             f'--trio "de_novo:INFO.gnomad_popmax_af <= {denovo_af_max} ' \
                             f'&& (trio_denovo(kid, mom, dad) || (variant.CHROM==\'chrX\' ' \
                             f'&& trio_x_linked_recessive(kid, dad, mom)))" ' \
                             f'--trio "homo:INFO.gnomad_nhomalt <= {nhomalt} ' \
                             f'&& trio_autosomal_recessive(kid, mom, dad)" ' \
                             f'--trio "x_linked_recessive:variant.CHROM==\'chrX\' ' \
                             f'&& INFO.gnomad_popmax_af <= {x_recessive_af_max} ' \
                             f'&& trio_x_linked_recessive(kid, dad, mom)"'
    slivar_noncomp_process = subprocess.run(shlex.split(slivar_noncomp_command),
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True)

    slivar_precomp_command = f'slivar_folder/slivar expr -v {vcf_file} -j {javascript_file} -p {ped_file} --pass-only -g {gnotate_file} -o {cache}{trio_name}_comp_prefiltered --info "INFO.gnomad_popmax_af <= {compf_af_max}"'
    slivar_precomp_process = subprocess.run(shlex.split(slivar_precomp_command),
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True)
    slivar_comp_command = f'slivar_folder/slivar compound-hets -v {cache}/{trio_name}_comp_prefiltered -p {ped_file} -o {cache}{trio_name}_comphets'
    slivar_comp_process = subprocess.run(shlex.split(slivar_comp_command),
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True)

    # slivar_noncomp_process = subprocess.run(shlex.split("echo test"),
    #                                 stdout=subprocess.PIPE,
    #                                 universal_newlines=True)
    # slivar_precomp_process = subprocess.run(shlex.split("echo test"),
    #                                 stdout=subprocess.PIPE,
    #                                 universal_newlines=True)
    # slivar_comp_process = subprocess.run(shlex.split("echo test"),
    #                                 stdout=subprocess.PIPE,
    #                                 universal_newlines=True)

    if not all(c == 0 for c in [slivar_noncomp_process.returncode, slivar_precomp_process.returncode, slivar_comp_process.returncode]):
        click.echo("There has some error with the slivar subprocess! Discontinuing further actions!")
    else:
        print("Slivar subprocesses successfull, starting scoring!")
        non_comphet_results = score_non_comphets(f"{cache}/{trio_name}_filtered", cache, trio_name, assembly)
        print("de_novos, x_linked & recessive scored!")
        comphet_results = score_comphets(f"{cache}/{trio_name}_comphets", cache, trio_name, output_path, assembly)
        print("comphets scored!")

        comphet_columns = list(comphet_results.columns)
        comphet_results.columns = ["variant"] + comphet_columns[1:]
        result_df = pd.concat([non_comphet_results, comphet_results.drop(columns="var_2")])
        result_df = result_df.sort_values(by="candidate_score", ascending=False).reset_index(drop=True)
        result_df.to_csv(f"{output_path}", index=False)
        print(result_df.head())
        shutil.rmtree(cache)

if __name__ == "__main__":
    score_vcf(shlex.split("-v /home/johann/AutoCaSc_core/slivar_folder/Bernt_VCF/conNDDcohort.annoated.vcf.gz "
                          "-p /home/johann/AutoCaSc_core/slivar_folder/Bernt_VCF/pedigree_MR144.ped"
                          "-g /home/johann/AutoCaSc_core/slivar_folder/gnomad.hg38.v2.zip"
                          "-j /home/johann/AutoCaSc_core/slivar_folder/slivar-functions.js"
                          "-o /home/johann/AutoCaSc_core/slivar_folder/scored/MR144_scored.csv"))
    # main(obj={})