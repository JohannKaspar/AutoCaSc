import copy
import os
import sys
import shlex
import subprocess
from io import StringIO
import click
from AutoCaSc import AutoCaSc, VERSION
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import re
import shutil
import pickle
import numpy as np


def create_bed_file(path=f"/home/johann/AutoCaSc/data/BED/Homo_sapiens.GRCh37.101.gtf"):
    print("create bed currently not working")
    """This function expands the exon positions by 50 to each side. Then one has to use the UNIX commands:
    sort -k 1,1 -k2,2n Homo_sapiens.GRCh38.bed | bedtools merge -i - > GRCh38.bed"""
    df = pd.read_csv(path, sep="\t", skiprows=5)
    df.columns = ["chr", "source", "type", "start", "stop", "unknown_1", "unknown_2", "unknown_3", "info"]
    exons = df.loc[df.type == "exon"]
    exons = exons.loc[exons["info"].str.contains("protein_coding")]
    exons = exons[["chr", "start", "stop"]]
    exons.start = exons.start - 50
    exons.stop = exons.stop + 50
    "chr" + exons.chr.astype(str)
    exons.to_csv(path + ".bed", sep="\t", index=False, header=False)


def thread_function_AutoCaSc_non_comp(params):
    vcf_chunk, assembly, ped_file, path_to_request_cache_dir = params
    vcf_chunk.reset_index(drop=True, inplace=True)
    return_df = pd.DataFrame(columns=["variant", "instance", "transcript_num"])

    # iterrating through the rows of the chunk, variant by variant
    for row_id, row in vcf_chunk.iterrows():
        try:
            chrom = row["#CHROM"]
        except pd.core.indexing.IndexingError:
            print("error")
        pos = row["POS"]
        ref = row["REF"]
        alt = row["ALT"]
        alt = alt.replace("*", "-")
        variant_vcf = ":".join(map(str, [chrom, pos, ref, alt]))
        if ped_file:
            quality_parameters = extract_quality_parameters(row, ped_file)
        else:
            quality_parameters = None

        if "comphet_side" in str(vcf_chunk.loc[row_id, "INFO"]):
            continue
        elif "de_novo" in str(vcf_chunk.loc[row_id, "INFO"]):
            inheritance = "de_novo"
        elif "homo" in str(vcf_chunk.loc[row_id, "INFO"]):
            inheritance = "homo"
        elif "x_linked" in str(vcf_chunk.loc[row_id, "INFO"]):
            inheritance = "x_linked"
        elif "autosomal_dominant" in str(vcf_chunk.loc[row_id, "INFO"]):
            inheritance = "ad_inherited"
        else:
            inheritance = "unknown"

        print(f"non_comph {row_id + 1} of {len(vcf_chunk)}: {variant_vcf}")

        autocasc_instance = AutoCaSc(variant=variant_vcf,
                                     inheritance=inheritance,
                                     assembly=assembly,
                                     path_to_request_cache_dir=path_to_request_cache_dir)

        if autocasc_instance.status_code == 200:
            autocasc_instance.calculate_candidate_score()

        i = len(return_df)
        return_df.loc[i, "variant"] = variant_vcf
        return_df.loc[i, "instance"] = autocasc_instance
        return_df.loc[i, "quality_parameters"] = quality_parameters
        return_df.loc[i, "transcript_num"] = 0

        # if there are multiple transcripts of same significance, calculate CaSc for all of them
        if len(autocasc_instance.high_priority_transcripts) > 1:
            for _transcript in autocasc_instance.high_priority_transcripts:
                alternative_instance = AutoCaSc(variant=variant_vcf,
                                                inheritance=inheritance,
                                                assembly=assembly,
                                                transcript=_transcript,
                                                path_to_request_cache_dir=path_to_request_cache_dir)
                if alternative_instance.status_code == 200:
                    alternative_instance.calculate_candidate_score()

                ix = len(return_df)
                return_df.loc[ix, "variant"] = variant_vcf + f" ({_transcript})"
                return_df.loc[ix, "instance"] = alternative_instance
                return_df.loc[ix, "transcript"] = _transcript
    return return_df


def thread_function_AutoCaSc_comp(params):
    # function that retrieves data
    variants, assembly, path_to_request_cache_dir = params
    return_dict = {}
    # iterrating through the variants
    for i, variant_vcf in enumerate(variants):
        print(f"comphet {i + 1} of {len(variants)}: {variant_vcf}")
        autocasc_instance = AutoCaSc(variant=variant_vcf,
                                     # inheritance="comphet",
                                     assembly=assembly,
                                     path_to_request_cache_dir=path_to_request_cache_dir)
        return_dict[variant_vcf] = autocasc_instance
    # multiple transcripts will be ignored by vcfAutoCaSc for compound heterozygous variants
    return return_dict


def score_non_comphets(filtered_vcf, cache, trio_name, assembly, ped_file, path_to_request_cache_dir, num_threads=2):
    # this loads the vcf containing all variants but compound heterozygous ones and converts it to a DataFrame
    with open(filtered_vcf, "r") as inp, open(
            f"{cache}/temp_{trio_name}.tsv",
            "w") as out:
        for row in inp:
            if "##" not in row:
                out.write(row)

    vcf_annotated = pd.read_csv(f"{cache}/temp_{trio_name}.tsv", sep="\t")

    vcf_chunks = [vcf_annotated.loc[round(i * len(vcf_annotated) / num_threads):
                                    round((i+1) * len(vcf_annotated) / num_threads), :] for i in range(num_threads)]

    # this starts annotation of all variants that are comphet
    with ProcessPoolExecutor() as executor:
        df_iterator = executor.map(thread_function_AutoCaSc_non_comp,
                                   zip(vcf_chunks,
                                       [assembly] * len(vcf_chunks),
                                       [ped_file] * len(vcf_chunks),
                                       [path_to_request_cache_dir] * len(vcf_chunks)))
        variant_df = pd.DataFrame()
        for _df in df_iterator:
            variant_df = pd.concat([variant_df, _df], ignore_index=True)
    return variant_df


def extract_quality_parameters(row, ped_file):
    """This function reads a ped_file, interprets it and extracts quality parameters from a row of a VCF file.
    """
    _, index_id, father_id, mother_id = interpret_pedigree(ped_file)
    quality_parameters = []
    if not re.findall(r"AC=[\w?]+(?=;)", row["INFO"]) == []:
        quality_parameters.append(re.findall(r"AC=[\w?]+(?=;)", row["INFO"])[0])
    if not re.findall(r"AF=[^;]+(?=;)", row["INFO"]) == []:
        quality_parameters.append(re.findall(r"AF=[^;]+(?=;)", row["INFO"])[0])
    quality_parameters.append(f'QUAL={row["QUAL"]}')
    if "GQ" in row["FORMAT"]:
        gq_position = row["FORMAT"].split("GQ")[0].count(":")
        quality_parameters.append(f'GQ_index={row[index_id].split(":")[gq_position]}')
    if "AD" in row["FORMAT"]:
        ad_position = row["FORMAT"].split("AD")[0].count(":")
        quality_parameters.append(f'AD_index={row[index_id].split(":")[ad_position]}'.replace(",",":"))
        quality_parameters.append(f'AD_father={row[father_id].split(":")[ad_position]}'.replace(",",":"))
        quality_parameters.append(f'AD_moth={row[mother_id].split(":")[ad_position]}'.replace(",",":"))
    if "DP" in row["FORMAT"]:
        dp_position = row["FORMAT"].split("DP")[0].count(":")
        quality_parameters.append(f'DP_index={row[index_id].split(":")[dp_position]}')
        quality_parameters.append(f'DP_father={row[father_id].split(":")[dp_position]}')
        quality_parameters.append(f'DP_moth={row[mother_id].split(":")[dp_position]}')
    return ";".join(quality_parameters)


def score_comphets(comphets_vcf, cache, trio_name, assembly, ped_file, path_to_request_cache_dir, num_threads=2):
    """This loads the vcf containing all compound heterozygous variants and converts it to a DataFrame.
    """
    with open(comphets_vcf, "r") as inp, open(f"{cache}/temp_{trio_name}.tsv", "w") as out:
        for row in inp:
            if "##" not in row:
                out.write(row)
    comphets_df = pd.read_csv(f"{cache}/temp_{trio_name}.tsv", sep="\t")

    # make comphet_cross_df
    comphet_cross_df = pd.DataFrame()
    for i, row in comphets_df.iterrows():
        chrom_1 = row["#CHROM"]
        pos_1 = row["POS"]
        ref_1 = row["REF"]
        alt_1 = row["ALT"]
        vcf_1 = ":".join(map(str, [chrom_1, pos_1, ref_1, alt_1]))

        if ped_file:
            quality_parameters = extract_quality_parameters(row, ped_file)
        else:
            quality_parameters = []

        info = row["INFO"]
        try:
            slivar_substrings = re.findall(r'(?:[^\/\,]+\/){6}[CTGA\-*]+', info.split("slivar_comphet=")[1])
            for substring in slivar_substrings:
                j = len(comphet_cross_df)
                chrom_2, pos_2, ref_2, alt_2 = substring.split("/")[-4:]
                vcf_2 = ":".join([chrom_2, pos_2, ref_2, alt_2])
                comphet_cross_df.loc[j, "var_1"] = vcf_1
                comphet_cross_df.loc[j, "var_2"] = vcf_2
                comphet_cross_df.loc[j, "quality_parameters"] = quality_parameters
        except IndexError:
            continue

    if len(comphet_cross_df) == 0:
        return comphet_cross_df

    all_variants = list(list(set(comphet_cross_df.var_1.to_list())) + list(set(comphet_cross_df.var_2.to_list())))
    variant_chunks = [all_variants[round(i * len(all_variants) / num_threads):round((i+1) * len(all_variants) / num_threads)] for i in range(num_threads)]

    with ProcessPoolExecutor() as executor:
        # retrieve data and inititate instances
        dicts_iterator = executor.map(thread_function_AutoCaSc_comp,
                                          zip(variant_chunks,
                                              [assembly] * len(variant_chunks),
                                              [path_to_request_cache_dir] * len(variant_chunks)))
        instances_dict = {}
        for _dict in dicts_iterator:
            instances_dict.update(_dict)

    for i, row in comphet_cross_df.iterrows():
        var_1 = row["var_1"]
        var_2 = row["var_2"]
        instance_1 = instances_dict.get(var_1)
        instance_2 = instances_dict.get(var_2)

        if instance_1.status_code == 200 and instance_2.status_code == 200:
            _instance = instance_1
            _instance.inheritance = "comphet"
            _other_instance = instance_2
            _other_instance.inheritance = "comphet"

            if _instance.transcript == _other_instance.transcript:
                _instance = copy.deepcopy(_instance)
                _instance.__dict__.pop("transcript_instances")
                _other_instance = copy.deepcopy(_other_instance)
                _other_instance.__dict__.pop("transcript_instances")
                _instance.other_autocasc_obj = _other_instance
                _instance.calculate_candidate_score()
            else:
                common_transcripts = list(set(_other_instance.affected_transcripts) & set(_instance.affected_transcripts))
                highest_score = 0
                for _transcript in common_transcripts:
                    _transcript_instance = copy.deepcopy(_instance)
                    try:
                        _transcript_instance.__dict__.pop("transcript_instances")
                    except KeyError:
                        pass
                    _other_transcript_instance = copy.deepcopy(_other_instance)
                    try:
                        _other_transcript_instance.__dict__.pop("transcript_instances")
                    except KeyError:
                        pass
                    _transcript_instance.other_autocasc_obj = _other_transcript_instance
                    _transcript_instance.calculate_candidate_score()
                    if _transcript_instance.candidate_score > highest_score:
                        _instance = copy.deepcopy(_transcript_instance)
                    elif _transcript_instance.candidate_score == highest_score:
                        if _transcript_instance.transcript in instance_1.canonical_transcripts:
                            _instance = copy.deepcopy(_transcript_instance)


            # elif _instance.transcript in _other_instance.affected_transcripts:
            #     _other_instance.transcript = _instance.transcript
            # elif _other_instance.transcript in _instance.affected_transcripts:
            #     _instance.transcript = _other_instance.transcript
            # else:
            #     for j in range(max([len(_instance.affected_transcripts), len(_other_instance.affected_transcripts)])):
            #         try:
            #             if _instance.affected_transcripts[j] in _other_instance.affected_transcripts:
            #                 _instance.transcript = _instance.affected_transcripts[j]
            #                 _other_instance.transcript = _instance.affected_transcripts[j]
            #                 break
            #         except IndexError:
            #             pass
            #         try:
            #             if _other_instance.affected_transcripts[j] in _instance.affected_transcripts:
            #                 _instance.transcript = _other_instance.affected_transcripts[j]
            #                 _other_instance.transcript = _other_instance.affected_transcripts[j]
            #                 break
            #         except IndexError:
            #             pass

        else:
            _instance = None
        comphet_cross_df.loc[i, "instance"] = _instance

    comphet_cross_df = comphet_cross_df.rename(columns={"var_1": "variant", "var_2": "other_variant"})
    return comphet_cross_df


def edit_java_script_functions(js_file, dp_filter=20, ab_filter=0.2, gq_filter=20):
    """This function reads the original slivar javascript file and edits the first line of code containing quality
    criteria.
    """
    with open(js_file, "r") as original_file, open(js_file + ".temp", "w") as new_file:
        new_file.write(f"var config = {{min_GQ: {gq_filter}, min_AB: {ab_filter}, min_DP: {dp_filter}}}")
        for line in original_file:
            if "var config" in line:
                continue
            new_file.write(line)
    js_file += ".temp"
    return js_file


def clean_up_duplicates(df):
    df = df.sort_values(by="variant", ascending=True)
    df = df.drop_duplicates(subset=list(df.columns)[1:], keep="first")
    return df


def interpret_pedigree(ped_file):
    """This function is used to check if parents are affected by the given phenotype and to extract the identifiers for
    both parents and the index."""
    pedigree = pd.read_csv(ped_file, sep="\t", header=None)
    # affected status = 2 means the individual is affected
    pedigree.columns = ["family_id", "individual_id", "paternal_id", "maternal_id", "sex", "affected_status"]

    individuals = pedigree.individual_id.to_list()
    parent_affected = False
    for i, row in pedigree.iterrows():
        if row["paternal_id"] in individuals and row["maternal_id"] in individuals:
            index_id = row["individual_id"]
        else:
            if row["sex"] == 1:
                father_id = row["individual_id"]
            else:
                mother_id = row["individual_id"]
            if row["affected_status"] == 2:
                parent_affected = True
    return parent_affected, index_id, father_id, mother_id


def execute_slivar(slivar_dir,
                   cache,
                   vcf_file,
                   bed_file,
                   ped_file,
                   gnotate_file,
                   javascript_file,
                   trio_name,
                   assembly,
                   quality,
                   x_recessive_af_max,
                   denovo_af_max,
                   autosomal_af_max,
                   comp_af_max,
                   ar_af_max,
                   gq_filter,
                   dp_filter,
                   ab_filter,
                   nhomalt,
                   pass_only,
                   deactivate_bed_filter):
    """This is a wrapper for slivar. It takes in different quality parameters and file locations and alters the command
    executing slivar accordingly.
    """
    if not bed_file:
        bed_file = str(sys.path[0]) + f"/data/BED/{assembly}.bed"

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

    javascript_file = edit_java_script_functions(javascript_file,
                                                 gq_filter=gq_filter,
                                                 dp_filter=dp_filter,
                                                 ab_filter=ab_filter)

    if float(x_recessive_af_max) == 0.0:
        x_recessive_af_max = "0.000000001"
    if float(denovo_af_max) == 0.0:
        denovo_af_max = "0.000000001"
    if float(comp_af_max) == 0.0:
        comp_af_max = "0.000000001"

    if pass_only:
        pass_expr = " && variant.FILTER == 'PASS'"
    else:
        pass_expr = ""
    slivar_noncomp_command = f'{slivar_dir} expr -v {vcf_file} -j {javascript_file} -p {ped_file} ' \
                                 f'--pass-only -g {gnotate_file} -o {cache}/{trio_name}_non_comphets ' \
                                 f'--info "INFO.gnomad_nhomalt <= {nhomalt} && INFO.impactful{pass_expr}" ' \
                                 f'--trio "homo:recessive(kid, mom, dad) && INFO.gnomad_popmax_af < {ar_af_max}" ' \
                                 f'--trio "x_linked_recessive:(variant.CHROM==\'X\' || variant.CHROM==\'chrX\') ' \
                                 f'&& INFO.gnomad_popmax_af <= {x_recessive_af_max} ' \
                                 f'&& x_recessive(kid, mom, dad)" ' \
                                 f'--trio "de_novo:INFO.gnomad_popmax_af <= {denovo_af_max} ' \
                                 '&& (denovo(kid, mom, dad) || ((variant.CHROM==\'X\' ' \
                                 '|| variant.CHROM==\'chrX\') ' \
                                 f'&& x_denovo(kid, mom, dad)))" ' \
                                 f'--trio "comphet_side:comphet_side(kid, mom, dad) ' \
                                 f'&& INFO.gnomad_popmax_af <= {comp_af_max}"'
    if interpret_pedigree(ped_file)[0]:
        slivar_noncomp_command += f' --trio "autosomal_dominant:INFO.gnomad_popmax_af <= {autosomal_af_max} ' \
                                 '&& trio_autosomal_dominant(kid, mom, dad)" '

    slivar_noncomp_process = subprocess.run(shlex.split(slivar_noncomp_command),
                                            stdout=subprocess.PIPE,
                                            universal_newlines=True)

    slivar_comp_command = f'{slivar_dir} compound-hets -v {cache}/{trio_name}_non_comphets ' \
                          f'-p {ped_file} -o {cache}/{trio_name}_comphets --sample-field comphet_side ' \
                          f'--sample-field de_novo'
    slivar_comp_process = subprocess.run(shlex.split(slivar_comp_command),
                                         stdout=subprocess.PIPE,
                                         universal_newlines=True)

    return all(c == 0 for c in [slivar_noncomp_process.returncode, slivar_comp_process.returncode])


def make_spreadsheet(merged_instances):
    """This function extracts informations from the scored variant instances and converts them to a clean dataframe.
    """
    result_df = pd.DataFrame()
    if not "other_variant" in merged_instances.columns:
        merged_instances.loc[:, "other_variant"] = ""
    for i, row in merged_instances.iterrows():
        _instance = row["instance"]
        if not _instance == "":
            result_df.loc[i, "variant"] = row["variant"]
            if row.other_variant != "":  # change this to not pd.isnull()?
                try:
                    result_df.loc[i, "other_variant"] = _instance.other_autocasc_obj.__dict__.get("variant")
                except AttributeError:
                    print("stop")
            result_df.loc[i, "gene_symbol"] = _instance.__dict__.get("gene_symbol")
            result_df.loc[i, "hgvsc"] = _instance.__dict__.get("hgvsc_change")
            result_df.loc[i, "hgvsp"] = _instance.__dict__.get("hgvsp_change")
            result_df.loc[i, "impact"] = _instance.__dict__.get("impact")
            result_df.loc[i, "inheritance"] = _instance.__dict__.get("inheritance")
            result_df.loc[i, "candidate_score"] = _instance.__dict__.get("candidate_score")
            result_df.loc[i, "transcript"] = _instance.__dict__.get("transcript")
            result_df.loc[i, "literature_score"] = _instance.__dict__.get("literature_score")
            result_df.loc[i, "CADD_phred"] = _instance.__dict__.get("cadd_phred") or 0
            result_df.loc[i, "status_code"] = _instance.__dict__.get("status_code")
            try:
                quality_parameters = row["quality_parameters"].split(";")
                for parameter in quality_parameters:
                    result_df.loc[i, parameter.split("=")[0]] = parameter.split("=")[1]
            except IndexError:
                pass

            if _instance.__dict__.get("factors"):
                try:
                    result_df.loc[i, "factors"] = " ".join([str(tup) for tup in _instance.__dict__.get("factors")])
                except TypeError:
                    print("stop")
    return result_df


def update_request_cache(path_to_request_cache_dir, assembly="GRCh37"):
    for api in ["gnomad", "vep"]:
        if os.path.isdir(f"{path_to_request_cache_dir}tmp/{api}"):
            with open(f"{path_to_request_cache_dir}{api}_requests_{assembly}", "rb") as requests_file:
                requests = pickle.load(requests_file)
                for entry in os.scandir(f"{path_to_request_cache_dir}tmp/{api}"):
                    with open(entry.path, "rb") as new_request_file:
                        try:
                            new_request = pickle.load(new_request_file)
                            requests.update(new_request)
                        except UnicodeDecodeError:
                            pass
            with open(f"{path_to_request_cache_dir}{api}_requests_{assembly}", "wb") as requests_file:
                pickle.dump(requests, requests_file)
            shutil.rmtree(f"{path_to_request_cache_dir}tmp/{api}")

def filter_xlinked_frequency(merged_instances, n_hemialt):
    drop_rows = []
    for i, row in merged_instances.iterrows():
        try:
            _instance = row["instance"]
            if _instance.inheritance == "x_linked":
                if _instance.male_count > n_hemialt:
                    drop_rows.append(i)
        except AttributeError:
            continue
    return merged_instances.drop(drop_rows).reset_index(drop=True)


def filter_blacklist(df, blacklist_path):
    with open(blacklist_path, "r") as file:
        blacklist = [line.strip() for line in file if not line[0] == "#"]
    for pattern in blacklist:
        df.loc[df.gene_symbol.str.contains(pattern, regex=True, na=False), "blacklist_filtered"] = True
        #df = df.loc[~df.gene_symbol.str.contains(pattern, regex=True, na=False)]
    return df.reset_index(drop=True)


def get_mim_number(gene, omim_morbid=None):
    omim_gene = omim_morbid.loc[omim_morbid.gene_symbol == gene]
    if omim_gene.empty:
        return ""
    else:
        return str(omim_gene.mim_number.to_list())


def mim_map(df, omim_morbid_path, column="gene_symbol"):
    # loading OMIM morbid dump
    with open(omim_morbid_path, "r") as raw_file:
        line_list = [line for line in raw_file if not line[0] in ["#", "[", "{", "?"]]

    omim_morbid = pd.read_csv(StringIO("\n".join(line_list)),
                              sep="\t",
                              header=None,
                              usecols=range(3))
    omim_morbid.columns = ["disease", "gene_symbol", "mim_gene_number"]
    omim_morbid.gene_symbol = omim_morbid.gene_symbol.apply(lambda x: x.split(",")[0] if "," in x else x)
    try:
        omim_morbid["mim_number"] = omim_morbid.disease.apply(lambda x: re.findall(", [0-9]{6}", x)[0].strip(", ") if re.findall(", [0-9]{6}", x) != [] else None)
    except IndexError:
        print("error indexing")
    omim_morbid = omim_morbid.dropna()

    for i in range(len(df)):
        _gene = df.loc[i, column]
        mim_number = get_mim_number(_gene, omim_morbid)
        df.loc[i, "mim_number"] = mim_number
    return df


def in_sysid(df, sysid_primary_path, sysid_candidates_path):
    sysid_primary = pd.read_csv(sysid_primary_path)["Gene symbol"].to_list()
    sysid_candidates = pd.read_csv(sysid_candidates_path)["Gene symbol"].to_list()
    for i, row in df.iterrows():
        gene_symbol = row["gene_symbol"]
        if gene_symbol in sysid_primary:
            df.loc[i, "sysid"] = "primary"
        elif gene_symbol in sysid_candidates:
            df.loc[i, "sysid"] = "candidates"
        else:
            df.loc[i, "sysid"] = ""
    return df


def filter_ac_mim(df, dp):
    if "sysid" in df.columns:
        temp = pd.concat([df.loc[df.sysid != ""], df.loc[df.mim_number == ""]]).drop_duplicates()
    else:
        temp = df.loc[df.mim_number == ""]
    temp.AC = pd.to_numeric(temp.AC, errors="coerce", downcast="unsigned").fillna(1)
    temp = temp.loc[~(temp.candidate_score == 0)]
    if "blacklist_filtered" in temp.columns:
        temp = temp.loc[temp.blacklist_filtered != True]

    temp = temp.loc[(pd.to_numeric(temp.DP_index, errors="coerce") > int(dp))
                    & (pd.to_numeric(temp.DP_moth, errors="coerce") > int(dp))
                    & (pd.to_numeric(temp.DP_father, errors="coerce") > int(dp))]

    temp = pd.concat(
        [
            temp.loc[(temp.inheritance == "de_novo") & (temp.AC == 1)],
            temp.loc[(temp.inheritance == "de_novo") & (temp.AC == 2) & (temp.variant.str[0].isin(["X", "Y"]))],
            temp.loc[(temp.inheritance == "homo") & (temp.AC <= 6)],  # homo --> AC == 4, +2 for being recessive
            temp.loc[(temp.inheritance == "comphet") & (temp.AC <= 4)],  # comphet --> AC == 2, +2 for being recessive
            temp.loc[(temp.inheritance == "x_linked") & (temp.AC <= 5)],  # x_linked --> AC == 3, +2 for being recessive
            temp.loc[(temp.inheritance == "ad_inherited") & (temp.AC == 2)],
        ],
        ignore_index=True)
    try:
        temp.loc[:, "autocasc_filter"] = "PASS"
    except ValueError:
        pass

    drop_indexes = []
    for i, row in temp.iterrows():
        if row.inheritance == "comphet":
            if row.other_variant not in temp.variant.to_list():
                drop_indexes.append(i)
    temp = temp.drop(temp.index[drop_indexes])
    temp.reset_index(drop=True, inplace=True)
    return temp


def add_ranks(df, dp):
    df[f"candidate_score"] = pd.to_numeric(df[f"candidate_score"],
                                                     errors="coerce")
    df.sort_values(f"candidate_score",
                   ascending=False,
                   inplace=True,
                   ignore_index=True)
    df.loc[:, f"rank"] = df.index
    df.loc[:, f"rank"] = df.loc[:, f"rank"].apply(lambda x: int(x+1))

    temp = filter_ac_mim(df, dp)

    if not temp.empty:
        temp.sort_values(f"candidate_score",
                         ascending=False,
                         inplace=True,
                         ignore_index=True)
        if not "other_variant" in temp.columns:
            temp.loc[:, f"rank_filtered"] = temp.index + 1
        else:
            temp.loc[:, f"rank_filtered"] = np.nan
            for i, row in temp.iterrows():
                if row["variant"] in temp.other_variant.to_list():
                    if pd.isnull(temp.loc[temp.other_variant == row["variant"], f"rank_filtered"].values[0]):
                        temp.loc[i, f"rank_filtered"] = temp[f"rank_filtered"].nunique() + 1
                    else:
                        temp.loc[i, f"rank_filtered"] = temp.loc[temp.other_variant == row["variant"],
                                                                           f"rank_filtered"].values[0]
                else:
                    temp.loc[i, f"rank_filtered"] = temp[f"rank_filtered"].nunique() + 1

        if "autocasc_filter" in temp.columns:
            merge_columns = ["rank", "rank_filtered", "autocasc_filter"]
        else:
            merge_columns = ["rank", "rank_filtered"]
        df = df.merge(temp[merge_columns],
                      on=f"rank",
                      how="left")
        df.sort_values("rank_filtered",
                          inplace=True)
    else:
        df["rank_filtered"] = np.nan
    df.loc[:, "rank_filtered"] = pd.to_numeric(df.loc[:, "rank_filtered"],
                                               downcast="unsigned",
                                               errors="ignore")
    return df


def post_scoring_polish(df,
                        omim_morbid_path,
                        dp=20,
                        blacklist_path=None,
                        sysid_primary_path=None,
                        sysid_candidates_path=None):
    if blacklist_path:
        df = filter_blacklist(df, blacklist_path)
    df = mim_map(df, omim_morbid_path)
    if sysid_primary_path and sysid_candidates_path:
        df = in_sysid(df, sysid_primary_path, sysid_candidates_path)
    df = add_ranks(df, dp)
    if "autocasc_filter" in df.columns:
        df.loc[df.autocasc_filter != "PASS", "autocasc_filter"] = "FAIL"
    else:
        df["autocasc_filter"] = "FAIL"
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
              required=False,
              type=click.Path(exists=True),
              help="Path to called vcf file.")
@click.option("--ped_file", "-p",
              type=click.Path(exists=True),
              help="Path to pedigree file.")
@click.option("--bed_file", "-b",
              type=click.Path(exists=True),
              help="Path to BED file.")
@click.option("--gnotate_file", "-g",
              type=click.Path(exists=True),
              help="Path to slivar gnotate file.")
@click.option("--javascript_file", "-j",
              type=click.Path(exists=True),
              help="Path to slivar javascript file.")
@click.option("--output_path", "-o",
              help="Output path for table with scored variants.")
@click.option("--denovo_af_max", "-daf",
              default="0.0000001",
              help="Popmax AF in gnomad for denovo variants.")
@click.option("--x_recessive_af_max", "-xaf",
              default="0.001",
              help="Popmax AF in gnomad for x-linked recessive variants.")
@click.option("--ar_af_max", "-raf",
              default="0.005",
              help="Popmax AF in gnomad for autosomal recessive variants.")
@click.option("--comp_af_max", "-caf",
              default="0.005",
              help="Popmax AF in gnomad for compound heterozygous variants.")
@click.option("--autosomal_af_max", "-aaf",
              default="0.0000001",
              help="Popmax AF in gnomad for autosomal dominant variants.")
@click.option("--nhomalt", "-nh",
              default="0",
              help="Max number of alternative sequence homozygotes in gnomad.")
@click.option("--n_hemialt", "-nhemi",
              default=0,
              help="Max number of alternative sequence hemizygotes in gnomad.")
@click.option("--quality", "-q",
              default="400",
              help="Minimum quality of variants.")
@click.option("--gq_filter", "-gq",
              default="10",
              help="Minimum GQ (quality) for variants.")
@click.option("--ab_filter", "-ab",
              default="0.2",
              help="Minimum allele balance for variants.")
@click.option("--dp_filter", "-dp",
              default="20",
              help="Minimum DP (reads) for variants.")
@click.option("--pass_only", "-pass",
              is_flag=True,
              default=False,
              help="Whether to use vants with VCF-filter 'PASS' only.")
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
              help="In case you want to skip BED-filtering.")
@click.option("--skip_slivar", "-ssli",
              is_flag=True,
              default=False,
              help="In case you dont want to run slivar yourself.")
@click.option("--vcf_comphet_path", "-vcf_ch",
              type=click.Path(exists=True),
              help="When skipping slivar, use this path to direct to VCf containing compound heterozygous variants.")
@click.option("--vcf_non_comphet_path", "-vcf_non_ch",
              type=click.Path(exists=True),
              help="When skipping slivar, use this path to direct to VCf containing NON compound heterozygous variants.")
@click.option("--path_to_request_cache_dir", "-req_cache",
              type=click.Path(exists=True),
              help="Path to directory for storing successfull API request responses for faster scoring in case of redundant requests.")
@click.option("--blacklist_path", "-blp",
              type=click.Path(exists=True),
              help="Path to list of Regex patterns. Genes matching the patterns will not be processed. (e.g. for excluding local artefacts)")
@click.option("--omim_morbid_path", "-omim",
              type=click.Path(exists=True),
              help="Path to OMIM morbid genes file downloaded from OMIM.")
@click.option("--sysid_primary_path", "-sys_prim",
              type=click.Path(exists=True),
              help="Path to .csv file from SysID containing list of NDD genes.")
@click.option("--sysid_candidates_path", "-sys_cand",
              type=click.Path(exists=True),
              help="Path to .csv file from SysID containing list of NDD candidates.")
@click.option("--trio_name", "-trio",
              help="Name of trio to use")
def score_vcf(vcf_file, ped_file, bed_file, gnotate_file, javascript_file, output_path,
              denovo_af_max, x_recessive_af_max, ar_af_max, comp_af_max, autosomal_af_max, nhomalt, n_hemialt,
              quality, gq_filter, ab_filter, dp_filter, pass_only, cache, slivar_dir, assembly, deactivate_bed_filter,
              skip_slivar, vcf_comphet_path, vcf_non_comphet_path, path_to_request_cache_dir,
              blacklist_path, omim_morbid_path, sysid_primary_path, sysid_candidates_path, trio_name):
    if not trio_name:
        trio_name = ped_file.split("/")[-1].split(".")[0]
    if skip_slivar:
        slivar_ok = True
    else:
        if not cache:
            cache = vcf_file.rstrip(trio_name) + "tmp"
            if not os.path.exists(cache):
                os.mkdir(cache)
        if not vcf_comphet_path:
            vcf_comphet_path = f"{cache}/{trio_name}_comphets"
        if not vcf_non_comphet_path:
            vcf_non_comphet_path = f"{cache}/{trio_name}_non_comphets"
        click.echo(f"Using cache directory for slivar filtered VCFs: {cache}")
        slivar_ok = execute_slivar(slivar_dir,
                                   cache,
                                   vcf_file,
                                   bed_file,
                                   ped_file,
                                   gnotate_file,
                                   javascript_file,
                                   trio_name,
                                   assembly,
                                   quality,
                                   x_recessive_af_max,
                                   denovo_af_max,
                                   autosomal_af_max,
                                   comp_af_max,
                                   ar_af_max,
                                   nhomalt,
                                   gq_filter,
                                   ab_filter,
                                   dp_filter,
                                   pass_only,
                                   deactivate_bed_filter)

    if not slivar_ok:
        click.echo("There has some error with the slivar subprocess! Discontinuing further actions!")
    else:
        print("Slivar subprocesses successfull, starting scoring!")
        if vcf_comphet_path:
            comphet_variant_instances = score_comphets(vcf_comphet_path,
                                                       cache,
                                                       trio_name,
                                                       assembly,
                                                       ped_file,
                                                       path_to_request_cache_dir)
            print("comphets scored!")
        else:
            comphet_variant_instances = pd.DataFrame()

        if vcf_non_comphet_path:
            non_comphet_variant_instances = score_non_comphets(vcf_non_comphet_path,
                                                               cache,
                                                               trio_name,
                                                               assembly,
                                                               ped_file,
                                                               path_to_request_cache_dir)
            print("non-comphet scored!")
        else:
            non_comphet_variant_instances = pd.DataFrame()

        merged_instances = pd.concat([non_comphet_variant_instances, comphet_variant_instances], ignore_index=True)
        merged_instances.fillna("", inplace=True)

        if merged_instances.empty:
            raise IOError("No variants left after filtering.")

        if n_hemialt is not None:
            merged_instances = filter_xlinked_frequency(merged_instances, n_hemialt)

        if merged_instances.empty:
            raise IOError("No variants left after filtering.")

        autocasc_df = make_spreadsheet(merged_instances)
        autocasc_df.fillna("-", inplace=True)
        autocasc_df = clean_up_duplicates(autocasc_df)

        if omim_morbid_path and sysid_primary_path and sysid_candidates_path:
            autocasc_df = post_scoring_polish(autocasc_df,
                                              omim_morbid_path,
                                              dp=dp_filter,
                                              blacklist_path=blacklist_path,
                                              sysid_primary_path=sysid_primary_path,
                                              sysid_candidates_path=sysid_candidates_path)
        autocasc_df["version"] = VERSION
        autocasc_df.to_csv(f"{output_path}",
                         index=False,
                         decimal=",",
                         sep="\t")
        #shutil.rmtree(cache)
        update_request_cache(path_to_request_cache_dir, assembly)
        #os.remove(javascript_file)


if __name__ == "__main__":
    # score_vcf(['-vcf_non_ch', '/home/johann/trio_scoring_results/synthetic_trios/2021-04-08/cache/ASH_sim01/ASH_a_non_comphets', '-vcf_ch', '/home/johann/trio_scoring_results/synthetic_trios/2021-04-08/cache/ASH_sim01/ASH_a_comphets', '-p', '/home/johann/PEDs/ASH_a.ped', '-g', '/home/johann/tools/slivar/gnotate/gnomad.hg37.zip', '-j', '/home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/slivar-functions.js', '-o', '/home/johann/trio_scoring_results/synthetic_trios/2021-04-08/ASH_sim01.csv', '-a', 'GRCh37', '-s', '/home/johann/tools/slivar/slivar', '-blp', '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/gene_blacklist.txt', '-omim', '/home/johann/PycharmProjects/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/OMIM_morbidmap.tsv', '-sys_prim', '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/sysid_primary_20210203.csv', '-sys_cand', '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/sysid_candidates_20210203.csv', '-ssli', '-dbed', '-req_cache', '/home/johann/PycharmProjects/AutoCaSc_project_folder/paper/data/', '--cache', '/home/johann/trio_scoring_results/synthetic_trios/2021-04-08/cache/ASH_sim01/', '-dp', '20', '-ab', '0.3', '-ssli']    )
    # score_vcf(['-vcf_non_ch', '/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/paper/data/tmp/ASH_a_non_comphets', '-vcf_ch', '/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/paper/data/tmp/ASH_a_comphets', '-p', '/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/paper/data/tmp/ASH_a.ped', '-o', '/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/paper/data/tmp/ASH_sim01.csv', '-a', 'GRCh37', '-blp', '/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/paper/data/gene_blacklist.txt', '-omim', '/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/OMIM_morbidmap.tsv', '-sys_prim', '/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/paper/data/sysid_primary_20210203.csv', '-sys_cand', '/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/paper/data/sysid_candidates_20210203.csv', '-ssli', '-dbed', '--cache', '/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/paper/data/tmp/', '-dp', '20', '-ab', '0.3', '-ssli']    )
    main(obj={})