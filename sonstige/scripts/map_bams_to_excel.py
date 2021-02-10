import pandas as pd
import xlrd
import re

bam_requests = pd.read_excel("/Users/johannkaspar/Documents/Promotion/AutoCaSc_Manuskript/family_sample_ids.xlsx",
                             sheet_name="bam_files")
bam_requests.loc[:, "mol_id"] = bam_requests.BAM_files.apply(lambda file_name: re.sub("\w{2}_S[0-9]+", "", file_name) if not "HUG" in file_name else file_name)
bam_requests.loc[:, "mol_id"] = bam_requests.mol_id.apply(lambda file_name: re.sub("-ready\.bam", "", file_name))
bam_requests.loc[:, "mol_id"] = bam_requests.mol_id.apply(lambda file_name: re.sub("-TWIST", "", file_name))

bam_requests = bam_requests[["BAM_files", "mol_id"]]

varvis_trio_results = pd.read_excel("/Users/johannkaspar/Documents/Promotion/AutoCaSc_Manuskript/TriosReal-List_2021-02-04.xlsx",
                                    sheet_name="TriosReal")

mol_ids = bam_requests.mol_id.to_list()
for i, row in varvis_trio_results.iterrows():
    try:
        mol_index = str(row['MOL-Nr. Index'])
        mol_mother = str(row['MOL-Nr. Mutter'])
        mol_father = str(row['MOL-Nr. Vater'])

        varvis_trio_results.loc[i, "index_bam"] = bam_requests.loc[bam_requests.BAM_files.str.contains(mol_index),
                                                                   "BAM_files"].values[0]
        varvis_trio_results.loc[i, "mother_bam"] = bam_requests.loc[bam_requests.BAM_files.str.contains(mol_mother),
                                                                   "BAM_files"].values[0]
        varvis_trio_results.loc[i, "father_bam"] = bam_requests.loc[bam_requests.BAM_files.str.contains(mol_father),
                                                                   "BAM_files"].values[0]

    except (TypeError, KeyError):
        pass

varvis_trio_results = varvis_trio_results[["Project_Identifier", "index_bam", "mother_bam", "father_bam"]]


varvis_trio_results = varvis_trio_results.set_index("Project_Identifier").T
varvis_trio_results.to_csv("/Users/johannkaspar/Documents/Promotion/AutoCaSc_Manuskript/TriosReal-List_bams_mapped.csv")
print("really done")
print("")