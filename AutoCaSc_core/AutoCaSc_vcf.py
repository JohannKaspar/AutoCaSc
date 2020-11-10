import os
import random
import shlex
import subprocess
from statistics import mean
import click
import time
from AutoCaSc import AutoCaSc
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import re
import shutil
from tenacity import retry, stop_after_attempt, wait_exponential

@retry(reraise=True, stop=(stop_after_attempt(10)|wait_exponential(multiplier=1, min=1, max=5)))
def iteration_func(variant_vcf, inheritance, assembly, transcript_num):
     instance = AutoCaSc(variant_vcf,
                               inheritance=inheritance,
                               assembly=assembly,
                               transcript_num=transcript_num)
     if instance.status_code == 200:
         return instance
     else:
         raise IOError("There has been an issue with a variant.")

def annotate_variant(variant_vcf, inheritance, assembly, transcript_num=None):
 try:
     instance = iteration_func(variant_vcf,
                               inheritance=inheritance,
                               assembly=assembly,
                               transcript_num=transcript_num)
 except IOError:
     instance = AutoCaSc(variant_vcf,
                         inheritance=inheritance,
                         assembly=assembly,
                         transcript_num=transcript_num)
 return instance


def create_bed_file(assembly="GRCh37", ensembl_version="101"):
    print("create bed currently not working")
    """This function expands the exon positions by 50 to each side. Then one has to use the UNIX commands:
    sort -k 1,1 -k2,2n Homo_sapiens.GRCh38.bed | bedtools merge -i - | awk '{print "chr"$0}' - > GRCh38.bed"""
    df = pd.read_csv(f"/home/johann/AutoCaSc/data/BED/Homo_sapiens.{assembly}.{ensembl_version}.gtf", sep="\t", skiprows=5)
    df.columns = ["chr", "source", "type", "start", "stop", "unknown_1", "unknown_2", "unknown_3", "info"]
    exons = df.loc[df.type == "exon"]
    exons = exons.loc[exons["info"].str.contains("protein_coding")]
    exons = exons[["chr", "start", "stop"]]
    exons.start = exons.start - 50
    exons.stop = exons.stop + 50
    "chr" + exons.chr.astype(str)
    exons.to_csv(f"/home/johann/AutoCaSc/data/BED/Homo_sapiens.{assembly}_unsorted.bed", sep="\t", index=False, header=False)

    # cache = "/home/johann/AutoCaSc/data/BED/temp/"
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
    vcf_chunk, assembly = param_tuple
    vcf_chunk.reset_index(drop=True, inplace=True)
    return_dict = {}
    for row_id in range(len(vcf_chunk)):
        try:
            chrom = vcf_chunk.loc[row_id, "#CHROM"]
        except pd.core.indexing.IndexingError:
            print("error")
        pos = vcf_chunk.loc[row_id, "POS"]
        ref = vcf_chunk.loc[row_id, "REF"]
        alt = vcf_chunk.loc[row_id, "ALT"]
        alt = alt.replace("*", "-")
        variant_vcf = ":".join(map(str, [chrom, pos, ref, alt]))
        if "de_novo" in str(vcf_chunk.loc[row_id, "INFO"]):
            inheritance = "de_novo"
        elif "homo" in str(vcf_chunk.loc[row_id, "INFO"]):
            inheritance = "homo"
        elif "x_linked" in str(vcf_chunk.loc[row_id, "INFO"]):
            inheritance = "x_linked"
        elif "comphet" in str(vcf_chunk.loc[row_id, "INFO"]):
            inheritance = "comphet"
        else:
            inheritance = "other"
        autocasc_instance = annotate_variant(variant_vcf,
                                             inheritance=inheritance,
                                             assembly=assembly)
        return_dict[variant_vcf] = autocasc_instance

        # autocasc_instance = AutoCaSc(variant_vcf, inheritance=inheritance, assembly=assembly)
        # return_dict[variant_vcf] = autocasc_instance
        # if autocasc_instance.status_code in [503, 497, 496, 201]:
        #     print("There has been an issue with a variant. Retrying...")
        #     for i in range(10):
        #         if autocasc_instance.status_code in [503, 497, 496, 201]:
        #             time.sleep(3)
        #             autocasc_instance = AutoCaSc(variant_vcf, inheritance=inheritance, assembly=assembly)
        #             return_dict[variant_vcf] = autocasc_instance
        #         else:
        #             break

        # if there are multiple transcripts of same significance, calculate CaSc for all of them
        if autocasc_instance.multiple_transcripts:
            for i in range(autocasc_instance.multiple_transcripts - 1):
                transcript_num = i + 1
                alternative_instance = annotate_variant(variant_vcf,
                                                        inheritance=inheritance,
                                                        assembly=assembly,
                                                        transcript_num=transcript_num)
                return_dict[variant_vcf + f"({transcript_num})"] = alternative_instance

    return return_dict


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
    num_threads = 10
    vcf_chunks = [vcf_annotated.loc[i * round(len(vcf_annotated) / num_threads):(i+1) * round(len(vcf_annotated) / num_threads),:] for i in range(num_threads)]

    # this starts annotation of all variants that are comphet
    with ThreadPoolExecutor() as executor:
        dicts_iterator = executor.map(thread_function_AutoCaSc_classic,
                                          zip(vcf_chunks, [assembly] * len(vcf_chunks)))
        instances_iterator = {}
        for _dict in dicts_iterator:
            instances_iterator.update(_dict)
        # instances_iterator = executor.map(thread_function_AutoCaSc_classic, vcf_chunks)
    result_df = pd.DataFrame()
    i = 0
    for _variant_vcf, _AutoCaSc_instance in instances_iterator.items():
        result_df.loc[i, "variant"] = _variant_vcf
        try:
            result_df.loc[i, "gene"] = _AutoCaSc_instance.gene_symbol
        except AttributeError:
            result_df.loc[i, "gene"] = "-"
        result_df.loc[i, "hgvsc"] = _AutoCaSc_instance.hgvsc_change
        result_df.loc[i, "transcript"] = _AutoCaSc_instance.transcript
        result_df.loc[i, "hgvsp"] = _AutoCaSc_instance.hgvsp_change
        result_df.loc[i, "candidate_score"] = _AutoCaSc_instance.candidate_score
        result_df.loc[i, "literature_score"] = _AutoCaSc_instance.literature_score
        result_df.loc[i, "CADD_phred"] = _AutoCaSc_instance.__dict__.get("cadd_phred") or 0
        result_df.loc[i, "impact"] = _AutoCaSc_instance.impact
        result_df.loc[i, "inheritance"] = _AutoCaSc_instance.inheritance
        result_df.loc[i, "status_code"] = _AutoCaSc_instance.status_code
        i += 1

    return result_df

def score_comphets(comphets_vcf, cache, trio_name, assembly, num_threads=10):
    #this loads the vcf containing all variants but compound heterozygous ones and converts it to a DataFrame
    with open(comphets_vcf, "r") as inp, open(f"{cache}/temp_{trio_name}.tsv", "w") as out:
        for row in inp:
            if "##" not in row:
                out.write(row)

    comphets_df = pd.read_csv(f"{cache}/temp_{trio_name}.tsv", sep="\t")
    # comphets_df = comphets_df.loc[comphets_df["#CHROM"] == "chr1"]
    # comphet_chunks = [comphets_df.loc[i, :] for i in range(len(comphets_df))]
    comphet_chunks = [comphets_df[round(i * len(comphets_df) / num_threads):round((i+1) * len(comphets_df) / num_threads)] for i in range(num_threads)]
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

    # score all variants
    with ThreadPoolExecutor() as executor:
        dicts_iterator = executor.map(thread_function_AutoCaSc_classic,
                                          zip(comphet_chunks, [assembly] * len(comphet_chunks)))
        instances_dict = {}
        for _dict in dicts_iterator:
            instances_dict.update(_dict)

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
        combined_candidate_score = mean([instance_1.candidate_score,
                                         instance_2.candidate_score])

        try:
            comphet_cross_df.loc[i, "gene"] = instance_1.gene_symbol
        except AttributeError:
            comphet_cross_df.loc[i, "gene"] = "-"

        comphet_cross_df.loc[i, "hgvsc"] = instance_1.hgvsc_change
        comphet_cross_df.loc[i, "transcript"] = instance_1.transcript
        comphet_cross_df.loc[i, "hgvsp"] = instance_1.hgvsp_change
        comphet_cross_df.loc[i, "candidate_score"] = combined_candidate_score
        comphet_cross_df.loc[i, "literature_score"] = instance_1.literature_score
        comphet_cross_df.loc[i, "CADD_phred"] = instance_1.__dict__.get("cadd_phred") or 0
        comphet_cross_df.loc[i, "impact"] = instance_1.impact
        comphet_cross_df.loc[i, "other_impact"] = instance_1.other_impact
        comphet_cross_df.loc[i, "other_transcript"] = instance_2.transcript
        comphet_cross_df.loc[i, "inheritance"] = instance_1.inheritance
        comphet_cross_df.loc[i, "status_code"] = instance_1.status_code

        # ToDo: delete this line
        comphet_cross_df.loc[i, "impact_score"] = instance_1.impact_score


    comphet_cross_df = comphet_cross_df.rename(columns={"var_1":"variant", "var_2":"other_variant"})
    # comphet_cross_df = comphet_cross_df.sort_values(by="candidate_score", ascending=False).reset_index(drop=True)
    # comphet_cross_df = comphet_cross_df.drop_duplicates(subset=["var_1"], keep="first")
    # comphet_cross_df.to_csv(f"{output_path}_scored_comphets.csv",
    #                         index=False)
    return comphet_cross_df

def edit_java_script_functions(js_file, dp_filter=20, gq_filter=10):
    with open(js_file, "r") as original_file, open(js_file + ".temp", "w") as new_file:
        for line in original_file:
            new_file.write(line)
        new_file.write(
            "function hq(sample) {\n"
                "if(sample.alts == -1) { return false; }\n"
                "if(sample.DP < 10) { return false; }\n"
                f"if(sample.GQ < {gq_filter}) {{ return false; }}\n"
                "if(sample.alts == 0) {\n"
                    f"if(sample.DP > {dp_filter} && sample.AB > 0.02) {{ return false; }}\n"
                    f"if(sample.DP <= {dp_filter} && sample.AD[1] > 1) {{ return false; }}\n"
                    "return true\n"
                "}\n"
                "if(sample.alts == 1) {\n"
                    "if(sample.AB < 0.25 || sample.AB > 0.75) { return false; }\n"
                    "return true\n"
                "}\n"
                "if(sample.alts == 2) {\n"
                    f"if(sample.DP > {dp_filter} && sample.AB < 0.98) {{ return false; }}\n"
                    f"if(sample.DP <= {dp_filter} && sample.AD[0] > 1) {{ return false; }}\n"
                    "return true\n"
                "}\n"
            "}\n")
    js_file += ".temp"
    return js_file

def clean_up_duplicates(df):
    df = df.sort_values(by="variant", ascending=True)
    df = df.drop_duplicates(subset=list(df.columns)[1:], keep="first")
    df = df.sort_values(by="candidate_score", ascending=False).reset_index(drop=True)
    return df

@click.group(invoke_without_command=True)  # Allow users to call our app without a command
@click.pass_context
@click.option('--verbose', '-v', is_flag=True, help="Increase output verbosity level")
def main(ctx, verbose):
    group_commands = ['candidates']
    """
            AutoCaSc is a command line tool that helps scoring
            variants of NDD patients for their plausible pathogenicity.\n
            example: python AutoCaSc.py batch -i input/file/path
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
              default="0.0000001",
              help="Popmax AF in gnomad for denovo variants.")
@click.option("--x_recessive_af_max", "-raf",
              default="0.0001",
              help="Popmax AF in gnomad for x-linked recessive variants.")
@click.option("--compf_af_max", "-raf",
              default="0.0001",
              help="Popmax AF in gnomad for compound heterozygous variants.")
@click.option("--nhomalt", "-nh",
              default="0",
              help="Max number of alternative sequence homozygotes in gnomad.")
@click.option("--quality", "-q",
              default="100",
              help="Minimum quality of variants.")
@click.option("--gq_filter", "-gq",
              default="10",
              help="Minimum GQ (quality) for variants.")
@click.option("--dp_filter", "-dp",
              default="20",
              help="Minimum DP (reads) for variants.")
@click.option("--cache", "-c",
              type=click.Path(exists=True),
              help="Path to cache for temporary data.")
@click.option("--slivar_dir", "-s",
              type=click.Path(exists=True),
              help="Path to slivar executable.")
@click.option("--assembly", "-a",
              default="GRCh37",
              help="Reference genome to use, either GRCh37 or GRCh38.")
@click.option("--non_coding", "-nc",
              is_flag=True,
              help="Add to deactivate filter for protein coding variants. (-50 and +50bp of exons)")
def score_vcf(vcf_file, ped_file, gnotate_file, javascript_file, output_path,
              denovo_af_max, x_recessive_af_max, compf_af_max, nhomalt, quality,
              gq_filter, dp_filter, cache, slivar_dir, assembly, non_coding):
    trio_name = vcf_file.split("/")[-1]
    #ToDo noch mal genau mit regex überlegen, damit es nicht "family_XYZ.vcf.gz_filtered" heißt
    if not cache:
        cache = vcf_file.rstrip(trio_name) + "tmp"
        if not os.path.exists(cache):
            os.mkdir(cache)
    click.echo(f"Using cache directory: {cache}")
    if not non_coding:
        with open(f"{cache}/vcf_filtered_temp", "wb") as vcf_filtered:
            bed_filter_command = f'bedtools intersect -wa -header -a {vcf_file} -b /home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_core/data/BED/{assembly}.bed'
            bed_filter_process = subprocess.run(shlex.split(bed_filter_command),
                                            stdout=vcf_filtered,
                                            universal_newlines=True)
        if bed_filter_process.returncode == 0:
            vcf_file = f"{cache}/vcf_filtered_temp"
            click.echo("BED-filtering successfull!")
        else:
            click.echo("There has been a problem filtering for protein coding variants only! Continuing with all variants.")

    javascript_file = edit_java_script_functions(javascript_file, gq_filter=gq_filter, dp_filter=dp_filter)
    if float(x_recessive_af_max) == 0.0:
        x_recessive_af_max = "0.000000001"
    if float(denovo_af_max) == 0.0:
        denovo_af_max = "0.000000001"
    if float(compf_af_max) == 0.0:
        compf_af_max = "0.000000001"

    slivar_noncomp_command = f'{slivar_dir} expr -v {vcf_file} -j {javascript_file} -p {ped_file} ' \
                             f'--pass-only -g {gnotate_file} -o {cache}/{trio_name}_non_comphets ' \
                             f'--info "variant.QUAL >= {quality}" ' \
                             f'--trio "homo:INFO.gnomad_nhomalt <= {nhomalt} ' \
                             f'&& trio_autosomal_recessive(kid, mom, dad)" ' \
                             '--trio "x_linked_recessive:variant.CHROM==\'chrX\' ' \
                             f'&& INFO.gnomad_popmax_af <= {x_recessive_af_max} ' \
                             f'&& trio_x_linked_recessive(kid, dad, mom)" ' \
                             f'--trio "de_novo:INFO.gnomad_popmax_af <= {denovo_af_max} ' \
                             '&& (trio_denovo(kid, mom, dad) || (variant.CHROM==\'chrX\' ' \
                             f'&& trio_x_linked_recessive(kid, dad, mom)))" '
    slivar_noncomp_process = subprocess.run(shlex.split(slivar_noncomp_command),
                            stdout=subprocess.PIPE,
                            universal_newlines=True)

    slivar_precomp_command = f'{slivar_dir} expr -v {vcf_file} -j {javascript_file} ' \
                             f'-p {ped_file} --pass-only ' \
                             f'-g {gnotate_file} -o {cache}/{trio_name}_comp_prefiltered ' \
                             f'--info "INFO.gnomad_popmax_af <= {compf_af_max} && variant.QUAL >= {quality}"'
    slivar_precomp_process = subprocess.run(shlex.split(slivar_precomp_command),
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True)
    slivar_comp_command = f'{slivar_dir} compound-hets -v {cache}/{trio_name}_comp_prefiltered ' \
                          f'-p {ped_file} -o {cache}/{trio_name}_comphets'
    slivar_comp_process = subprocess.run(shlex.split(slivar_comp_command),
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True)

    # slivar_noncomp_process = subprocess.run(shlex.split("echo test"))
    # slivar_precomp_process = subprocess.run(shlex.split("echo test"))
    # slivar_comp_process = subprocess.run(shlex.split("echo test"))

    if not all(c == 0 for c in [slivar_noncomp_process.returncode, slivar_precomp_process.returncode, slivar_comp_process.returncode]):
        click.echo("There has some error with the slivar subprocess! Discontinuing further actions!")
    else:
        print("Slivar subprocesses successfull, starting scoring!")
        non_comphet_results = score_non_comphets(f"{cache}/{trio_name}_non_comphets", cache, trio_name, assembly)
        print("de_novos, x_linked & recessive scored!")
        comphet_results = score_comphets(f"{cache}/{trio_name}_comphets", cache, trio_name, assembly)
        print("comphets scored!")

        comphet_columns = list(comphet_results.columns)
        comphet_results.columns = ["variant"] + comphet_columns[1:]
        # result_df = pd.concat([non_comphet_results, comphet_results.drop(columns="var_2")])
        result_df = pd.concat([non_comphet_results, comphet_results])
        result_df = clean_up_duplicates(result_df)
        result_df.to_csv(f"{output_path}", index=False)
        # print(result_df.head())
        shutil.rmtree(cache)
        os.remove(javascript_file)


if __name__ == "__main__":
    # score_vcf(shlex.split("-v /home/johann/Bernt_VCF/bernt_vcf_coding "
    #                       "-p /home/johann/Bernt_VCF/pedigrees/pedigree_MR144.ped "
    #                       "-g /home/johann/tools/slivar/gnotate/gnomad.hg38.v2.zip "
    #                       "-j /home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_core/data/slivar-functions.js "
    #                       "-o /home/johann/Bernt_VCF/scored/MR144_scored.csv "
    #                       "-a GRCh38 "
    #                       "-s /home/johann/tools/slivar/slivar "
    #                       "-nc"))
    main(obj={})

