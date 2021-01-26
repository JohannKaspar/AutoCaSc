import os
import sys
import random
import shlex
import subprocess
from io import StringIO
from statistics import mean
import click
import time
from AutoCaSc import AutoCaSc
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import re
import shutil
from tenacity import retry, stop_after_attempt, wait_exponential
import pickle

@retry(reraise=True, stop=stop_after_attempt(10), wait=wait_exponential(multiplier=1, min=1, max=5))
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


def load_omim_morbid(path="/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/OMIM_morbidmap.tsv"):
    # loading OMIM morbid dump
    with open(path, "r") as raw_file:
        filtered_file = ""
        for line in raw_file:
            if not line[0] in ["#", "[", "{"]:
                filtered_file += line

    omim_morbid = pd.read_csv(StringIO(filtered_file),
                              sep="\t",
                              header=None,
                              usecols=range(3))
    omim_morbid.columns = ["disease", "gene_symbol", "mim_number"]
    omim_morbid.gene_symbol = omim_morbid.gene_symbol.apply(lambda x: x.split(",")[0] if "," in x else x)
    return omim_morbid


def get_mim_number(gene, omim_morbid=None):
    if omim_morbid == None:
        omim_morbid = load_omim_morbid()
    omim_gene = omim_morbid.loc[omim_morbid.gene_symbol == gene]
    if omim_gene.empty:
        return ""
    else:
        return omim_gene.mim_number.to_list()

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

def thread_function_AutoCaSc_non_comp(param_tuple):
    vcf_chunk, assembly = param_tuple
    vcf_chunk.reset_index(drop=True, inplace=True)
    return_df = pd.DataFrame(columns=["variant", "instance", "transcript_num"])

    # iterrating through the rows of the chunk, variant by variant
    for row_id in range(len(vcf_chunk)):
        print(row_id)
        try:
            chrom = vcf_chunk.loc[row_id, "#CHROM"]
        except pd.core.indexing.IndexingError:
            print("error")
        pos = vcf_chunk.loc[row_id, "POS"]
        ref = vcf_chunk.loc[row_id, "REF"]
        alt = vcf_chunk.loc[row_id, "ALT"]
        alt = alt.replace("*", "-")
        variant_vcf = ":".join(map(str, [chrom, pos, ref, alt]))


        # todo change this back
        #if "," in variant_vcf:
        #    continue


        if "de_novo" in str(vcf_chunk.loc[row_id, "INFO"]):
            inheritance = "de_novo"
        elif "homo" in str(vcf_chunk.loc[row_id, "INFO"]):
            inheritance = "homo"
        elif "x_linked" in str(vcf_chunk.loc[row_id, "INFO"]):
            inheritance = "x_linked"
        elif "autosomal_dominant" in str(vcf_chunk.loc[row_id, "INFO"]):
            inheritance = "ad_inherited"
        else:
            inheritance = "unknown"

        autocasc_instance = AutoCaSc(variant=variant_vcf,
                                     inheritance=inheritance,
                                     assembly=assembly)

        if autocasc_instance.status_code == 200:
            autocasc_instance.calculate_candidate_score()

        i = len(return_df)
        return_df.loc[i, "variant"] = variant_vcf
        return_df.loc[i, "instance"] = autocasc_instance
        return_df.loc[i, "transcript_num"] = 0


        # if there are multiple transcripts of same significance, calculate CaSc for all of them
        if autocasc_instance.multiple_transcripts:
            for j in range(autocasc_instance.multiple_transcripts - 1):
                transcript_num = j + 1
                alternative_instance = AutoCaSc(variant=variant_vcf,
                                                inheritance=inheritance,
                                                assembly=assembly,
                                                transcript_num=transcript_num)
                if autocasc_instance.status_code == 200:
                    autocasc_instance.calculate_candidate_score()

                ix = len(return_df)
                return_df.loc[ix, "variant"] = variant_vcf
                return_df.loc[ix, "instance"] = alternative_instance
                return_df.loc[ix, "transcript_num"] = transcript_num

    return return_df


def thread_function_AutoCaSc_comp(param_tuple):
    variants, assembly = param_tuple
    return_dict = {}

    # iterrating through the variants
    for variant_vcf in variants:

        print(f"working {variant_vcf}")
        autocasc_instance = AutoCaSc(variant=variant_vcf,
                                     inheritance="comphet",
                                     assembly=assembly)
        return_dict[variant_vcf] = autocasc_instance
    # multiple transcripts will be ignored for compound heterozygous variants
    return return_dict


def score_non_comphets(filtered_vcf, cache, trio_name, assembly, num_threads=1):
    # this loads the vcf containing all variants but compound heterozygous ones and converts it to a DataFrame
    with open(filtered_vcf, "r") as inp, open(
            f"{cache}/temp_{trio_name}.tsv",
            "w") as out:
        for row in inp:
            if "##" not in row:
                out.write(row)

    vcf_annotated = pd.read_csv(f"{cache}/temp_{trio_name}.tsv", sep="\t")
    # vcf_annotated = vcf_annotated.loc[vcf_annotated["#CHROM"] == "chr1"]


    vcf_chunks = [vcf_annotated.loc[i * round(len(vcf_annotated) / num_threads):(i+1) * round(len(vcf_annotated) / num_threads),:] for i in range(num_threads)]

    # this starts annotation of all variants that are comphet
    with ThreadPoolExecutor() as executor:
        df_iterator = executor.map(thread_function_AutoCaSc_non_comp,
                                   zip(vcf_chunks, [assembly] * len(vcf_chunks)))
        variant_df = pd.DataFrame()
        for _df in df_iterator:
        # instances_iterator = executor.map(thread_function_AutoCaSc_classic, vcf_chunks)
            variant_df = pd.concat([variant_df, _df], ignore_index=True)

    return variant_df

def score_comphets(comphets_vcf, cache, trio_name, assembly, num_threads=1):
    #this loads the vcf containing all compound heterozygous variants and converts it to a DataFrame
    with open(comphets_vcf, "r") as inp, open(f"{cache}/temp_{trio_name}.tsv", "w") as out:
        for row in inp:
            if "##" not in row:
                out.write(row)
    comphets_df = pd.read_csv(f"{cache}/temp_{trio_name}.tsv", sep="\t")

    # comphet_chunks = [comphets_df[round(i * len(comphets_df) / num_threads):round((i+1) * len(comphets_df) / num_threads)] for i in range(num_threads)]

    # make comphet_cross_df
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
            j = len(comphet_cross_df)
            chrom_2, pos_2, ref_2, alt_2 = substring.split("/")[-4:]
            vcf_2 = ":".join([chrom_2, pos_2, ref_2, alt_2])

            # todo change this back
            #if "," in vcf_1 + vcf_2:
            #    continue

            comphet_cross_df.loc[j, "var_1"] = vcf_1
            comphet_cross_df.loc[j, "var_2"] = vcf_2

    all_variants = list(list(set(comphet_cross_df.var_1.to_list())) + list(set(comphet_cross_df.var_2.to_list())))
    variant_chunks = [all_variants[round(i * len(all_variants) / num_threads):round((i+1) * len(all_variants) / num_threads)] for i in range(num_threads)]

    # score all variants
    with ThreadPoolExecutor() as executor:
        dicts_iterator = executor.map(thread_function_AutoCaSc_comp,
                                          zip(variant_chunks, [assembly] * len(variant_chunks)))
        instances_dict = {}
        for _dict in dicts_iterator:
            instances_dict.update(_dict)

    """with open(cache + "/comphet_variant_instances", "wb") as comphet_out:
        pickle.dump(instances_dict, comphet_out)

    with open(cache + "/comphet_variant_instances", "rb") as comphet_in:
        instances_dict = pickle.load(comphet_in)"""



    for i, row in comphet_cross_df.iterrows():
        var_1 = row["var_1"]
        var_2 = row["var_2"]
        instance_1 = instances_dict.get(var_1)
        instance_2 = instances_dict.get(var_2)

        if instance_1.status_code == 200 and instance_2.status_code == 200:
            _instance = copy.deepcopy(instance_1)
            _instance.other_autocasc_obj = copy.deepcopy(instance_2)
            _instance.calculate_candidate_score()
            comphet_cross_df.loc[i, "instance"] = _instance

    comphet_cross_df = comphet_cross_df.rename(columns={"var_1":"variant", "var_2":"other_variant"})
    return comphet_cross_df

def edit_java_script_functions(js_file, dp_filter=20, gq_filter=10):
    with open(js_file, "r") as original_file, open(js_file + ".temp", "w") as new_file:
        for line in original_file:
            if "function hq(sample)" in line:
                break
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
                    "if(sample.AB < 0.2 || sample.AB > 0.8) { return false; }\n"
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
    return df

def parent_affected(ped_file):
    pedigree = pd.read_csv(ped_file, sep="\t", header=None)
    # affected status = 2 means the individual is affected
    pedigree.columns = ["family_id", "individual_id", "paternal_id", "maternal_id", "sex", "affected_status"]
    for i in range(len(pedigree)):
        if pedigree.loc[i, "paternal_id"] not in pedigree.individual_id.to_list() \
                and pedigree.loc[i, "maternal_id"] not in pedigree.individual_id.to_list():
            if pedigree.loc[i, "affected_status"] == 2:
                print("parent affected")
                return True
    print("parent unaffected")
    return False

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
@click.option("--bed_file", "-b",
              type=click.Path(exists=True),
              help="Path to BED file.")
@click.option("--gnotate_file", "-g",
              required=False,
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
@click.option("--autosomal_af_max", "-aaf",
              default="0.0000001",
              help="Popmax AF in gnomad for autosomal dominant variants.")
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
@click.option("--deactivate_bed_filter", "-dbed",
              is_flag=True,
              default=False,
              help="Reference genome to use, either GRCh37 or GRCh38.")
@click.option("--skip_slivar", "-ssli",
              is_flag=True,
              default=False,
              help="Just for debugging purposes.")
@click.option("--mim_flag", "-mim",
              is_flag=True,
              default=True,
              help="Map MIM numbers.")
# todo change mim flag to path to morbid file
# Todo delete skip-slivar option
def score_vcf(vcf_file, ped_file, bed_file, gnotate_file, javascript_file, output_path,
              denovo_af_max, x_recessive_af_max, compf_af_max, autosomal_af_max, nhomalt,
              quality, gq_filter, dp_filter, cache, slivar_dir, assembly, deactivate_bed_filter,
              skip_slivar, mim_flag):
    trio_name = vcf_file.split("/")[-1]
    #ToDo noch mal genau mit regex überlegen, damit es nicht "family_XYZ.vcf.gz_filtered" heißt
    if not cache:
        cache = vcf_file.rstrip(trio_name) + "tmp"
        if not os.path.exists(cache):
            os.mkdir(cache)
    click.echo(f"Using cache directory: {cache}")

    if not bed_file:
        bed_file = str(sys.path[0]) + f"/data/BED/{assembly}.bed"
        # todo debug this line to see if it works

    # this part is creating a copy of the VCF containing only coding variants
    if not deactivate_bed_filter:
        with open(f"{cache}/vcf_filtered_temp", "wb") as vcf_filtered:
            bed_filter_command = f'bedtools intersect -wa -header -a {vcf_file} -b {bed_file}'
            bed_filter_process = subprocess.run(shlex.split(bed_filter_command),
                                            stdout=vcf_filtered,
                                            universal_newlines=True)
        if bed_filter_process.returncode == 0:
            vcf_file = f"{cache}/vcf_filtered_temp"
            click.echo("BED-filtering successfull!")
        else:
            click.echo("There has been a problem filtering for protein coding variants only! Continuing with all variants.")

    javascript_file = edit_java_script_functions(javascript_file, gq_filter=gq_filter, dp_filter=dp_filter)

    if mim_flag:
        omim_morbid = load_omim_morbid()

    if float(x_recessive_af_max) == 0.0:
        x_recessive_af_max = "0.000000001"
    if float(denovo_af_max) == 0.0:
        denovo_af_max = "0.000000001"
    if float(compf_af_max) == 0.0:
        compf_af_max = "0.000000001"

    if skip_slivar:
        slivar_noncomp_process = subprocess.run(shlex.split("echo test"))
        slivar_precomp_process = subprocess.run(shlex.split("echo test"))
        slivar_comp_process = subprocess.run(shlex.split("echo test"))
    else:
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

        if parent_affected(ped_file):
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
                                 f'&& trio_x_linked_recessive(kid, dad, mom)))" ' \
                                 f'--trio "autosomal_dominant:INFO.gnomad_popmax_af <= {autosomal_af_max} ' \
                                 '&& trio_autosomal_dominant(kid, mom, dad)" '
        else:
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



    if not all(c == 0 for c in [slivar_noncomp_process.returncode, slivar_precomp_process.returncode, slivar_comp_process.returncode]):
        click.echo("There has some error with the slivar subprocess! Discontinuing further actions!")
    else:
        print("Slivar subprocesses successfull, starting scoring!")
        comphet_variant_instances = score_comphets(f"{cache}/{trio_name}_comphets", cache, trio_name, assembly)
        print("comphets scored!")
        non_comphet_variant_instances = score_non_comphets(f"{cache}/{trio_name}_non_comphets", cache, trio_name, assembly)
        print("de_novos, x_linked & recessive scored!")


        with open(cache + "/non_comphet_variant_instances_scored", "wb") as non_comphet_out:
            pickle.dump(non_comphet_variant_instances, non_comphet_out)
        with open(cache + "/comphet_variant_instances_scored", "wb") as comphet_out:
            pickle.dump(comphet_variant_instances, comphet_out)

        with open(cache + "/non_comphet_variant_instances_scored", "rb") as non_comphet_out:
            non_comphet_variant_instances = pickle.load(non_comphet_out)
        with open(cache + "/comphet_variant_instances_scored", "rb") as comphet_out:
            comphet_variant_instances = pickle.load(comphet_out)

        merged_instances = pd.concat([non_comphet_variant_instances, comphet_variant_instances], ignore_index=True)
        merged_instances.fillna("", inplace=True)

        #merged_instances = comphet_variant_instances
        result_df = pd.DataFrame()

        for i, row in merged_instances.iterrows():
            _instance = row["instance"]
            if not _instance == "":
                result_df.loc[i, "variant"] = row["variant"]
                if row.other_variant != "":  # if it is not np.nan
                    try:
                        result_df.loc[i, "other_variant"] = _instance.other_autocasc_obj.__dict__.get("variant")
                    except AttributeError:
                        print("stop")
                result_df.loc[i, "gene_symbol"] = _instance.__dict__.get("gene_symbol")

                result_df.loc[i, "hgvsc"] = _instance.__dict__.get("hgvsc_change")
                result_df.loc[i, "hgvsp"] = _instance.__dict__.get("hgvsp_change")
                result_df.loc[i, "impact"] = _instance.__dict__.get("impact")
                result_df.loc[i, "inheritance"] = _instance.__dict__.get("inheritance")
                result_df.loc[i, "candidate_score_v1"] = _instance.__dict__.get("candidate_score_v1")
                result_df.loc[i, "candidate_score_v2"] = _instance.__dict__.get("candidate_score_v2")
                result_df.loc[i, "candidate_score_v3"] = _instance.__dict__.get("candidate_score_v3")
                result_df.loc[i, "transcript"] = _instance.__dict__.get("transcript")
                result_df.loc[i, "literature_score"] = _instance.__dict__.get("literature_score")
                result_df.loc[i, "CADD_phred"] = _instance.__dict__.get("cadd_phred") or 0
                result_df.loc[i, "status_code"] = _instance.__dict__.get("status_code")

                if _instance.__dict__.get("factors"):
                    try:
                        result_df.loc[i, "factors"] = " ".join([str(tup) for tup in _instance.__dict__.get("factors")])
                    except TypeError:
                        print("stop")

                if mim_flag:
                    result_df.loc[i, "mim_number"] = get_mim_number(_instance.__dict__.get("gene_symbol"), omim_morbid)

        # comphet_columns = list(comphet_results.columns)
        # comphet_results.columns = ["variant"] + comphet_columns[1:]
        # result_df = pd.concat([non_comphet_results, comphet_results.drop(columns="var_2")])
        # result_df = pd.concat([non_comphet_results, comphet_results])

        result_df.fillna("-", inplace=True)
        result_df = clean_up_duplicates(result_df)
        #result_df = result_df.sort_values("candidate_score", ascending=False)
        result_df.to_csv(f"{output_path}", index=False)

        shutil.rmtree(cache)
        os.remove(javascript_file)


if __name__ == "__main__":
    """score_vcf(shlex.split("-v /home/johann/Bernt_VCF/bernt_vcf_coding "
                      "-p /home/johann/Bernt_VCF/pedigrees/pedigree_MR144.ped "
                       "-g /home/johann/tools/slivar/gnotate/gnomad.hg38.v2.zip "
                       "-j /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/slivar-functions.js "
                      "-o /home/johann/Bernt_VCF/scored/MR144_scored.csv "
                      "-a GRCh38 "
                      "-s /home/johann/tools/slivar/slivar "
                      "-ssli "
                      "-dbed "
                      "-nc "
                      ))"""
    """score_vcf(shlex.split("-v /mnt/raid/users/johann/VCFs/vcf_filtered_temp "
                      "-p /home/johann/PycharmProjects/AutoCaSc_project_folder/sonstige/data/ped_files/L16-0467.ped "
                       "-g /home/johann/tools/slivar/gnotate/gnomad.hg38.v2.zip "
                       "-j /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/slivar-functions.js "
                      "-o /mnt/raid/users/johann/VCFs/scored/L16-0467.csv "
                      "-a GRCh38 "
                      "-s /home/johann/tools/slivar/slivar "
                      #"-ssli "
                      "-dbed "
                      ))"""
    """score_vcf(shlex.split("-v /home/johann/VCFs/modified_VCFs/annotated/ASH_sim01.vcf.gz "
                      "-p /home/johann/PEDs/ASH_a.ped "
                       "-g /home/johann/tools/slivar/gnotate/gnomad.hg37.zip "
                       "-j /home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/slivar-functions.js "
                      "-o /home/johann/trio_scoring_results/ASH_sim01.csv "
                      "-a GRCh37 "
                      "-s /home/johann/tools/slivar/slivar "
                      "-ssli "
                      "-dbed "
                      ))"""
    main(obj={})