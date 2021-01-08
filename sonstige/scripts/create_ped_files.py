import pandas as pd
import xlrd


def load_sample_excel():
    sample_ids = pd.read_excel("/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/sonstige/data/family_sample_ids.xlsx",
                               sheet_name="Request_true_id")
    sample_ids.columns = ["_", "index_id", "mother_id", "father_id", "family_id"]
    sample_ids.family_id = sample_ids.family_id.apply(lambda x: x.rsplit("-", 1)[0])
    sample_ids.drop(columns=["_"], inplace=True)
    return sample_ids

def create_ped_files():
    for i, row in sample_ids.iterrows():
        index_id = str(row["index_id"])
        mother_id = str(row["mother_id"])
        father_id = str(row["father_id"])
        family_id = str(row["family_id"])

        with open("/Users/johannkaspar/Documents/Promotion/AutoCaSc_project_folder/sonstige/data/ped_files/" + f"{family_id}.ped", "w") as _ped_file:
            _ped_file.write("\n".join(["\t".join([family_id, index_id, father_id, mother_id, "other", "2"]),
                                       "\t".join([family_id, father_id, "0", "0", "1", "1"]),
                                       "\t".join([family_id, mother_id, "0", "0", "2", "1"])]))

sample_ids = load_sample_excel()
create_ped_files()
print("done")