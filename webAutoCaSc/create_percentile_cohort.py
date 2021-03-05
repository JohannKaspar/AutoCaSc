from AutoCaSc_core.AutoCaSc import AutoCaSc
from AutoCaSc_core.AutoCaSc_vcf import mim_map, get_mim_number
import pandas as pd


candidate_table = pd.read_csv("/Users/johannkaspar/Documents/Promotion/AutoCaSc_analytics/data/candidate-scores_rerun_20210305.csv")

candidate_table = mim_map(candidate_table,
                          omim_morbid_path="/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/webAutoCaSc/AutoCaSc_core/data/OMIM_morbidmap.tsv",
                          column="HGNC_Symbol")

sysid = pd.concat([
    pd.read_csv("/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/AutoCaSc_maintenance/data/sysid/sysid_primary.csv"),
    pd.read_csv("/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/AutoCaSc_maintenance/data/sysid/sysid_candidates.csv"),
])
candidate_table = candidate_table.loc[~candidate_table.HGNC_Symbol.isin(sysid["Gene symbol"].to_list())]
print(candidate_table.candidate_score.dropna().to_list())