from AutoCaSc_core.vcfAutoCaSc import mim_map
import pandas as pd

# this was used to calculate the percentiles for webAutoCaSc
# "This variant scored better than XYZ% of all our candidates."
candidate_table = pd.read_csv("/Users/johannkaspar/Documents/Promotion/AutoCaSc_analytics/data/candidate-scores_rerun_20210305.csv")

candidate_table = mim_map(candidate_table,
                          omim_morbid_path="/AutoCaSc_core/data/OMIM_morbidmap.tsv",
                          column="HGNC_Symbol")

sysid = pd.concat([
    pd.read_csv("/update_data/data/sysid/sysid_primary.csv"),
    pd.read_csv("/update_data/data/sysid/sysid_candidates.csv"),
])
candidate_table = candidate_table.loc[~candidate_table.HGNC_Symbol.isin(sysid["Gene symbol"].to_list())]
print(candidate_table.candidate_score.dropna().to_list())