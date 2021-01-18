import pandas as pd
import xlrd


def load_sample_excel():
    sample_ids = pd.read_excel("/home/johann/PycharmProjects/AutoCaSc_project_folder/sonstige/data/family_sample_ids.xlsx",
                               sheet_name="Request_true_id")
    sample_ids.columns = ["_", "index_id", "mother_id", "father_id", "family_id"]
    sample_ids.family_id = sample_ids.family_id.apply(lambda x: x.rsplit("-", 1)[0])
    sample_ids.drop(columns=["_"], inplace=True)
    return sample_ids

def create_ped_files():
    for i, row in sample_ids.iterrows():
        index_id = id_dict.get(str(row["index_id"]))
        mother_id = id_dict.get(str(row["mother_id"]))
        father_id = id_dict.get(str(row["father_id"]))
        family_id = str(row["family_id"])

        try:
            with open("/home/johann/PycharmProjects/AutoCaSc_project_folder/sonstige/data/ped_files/" + f"{family_id}.ped", "w") as _ped_file:
                _ped_file.write("\n".join(["\t".join([family_id, index_id, father_id, mother_id, "1", "2"]),
                                           "\t".join([family_id, father_id, "0", "0", "1", "1"]),
                                           "\t".join([family_id, mother_id, "0", "0", "2", "1"])]))
        except TypeError:
            print("stop")

def extract_id(input):
    if "23865" in input:
        print("stop")
    if "HUG" in input:
        return input
    if "_" not in input:
        return input
    part_1, part_2 = input.split("_")
    if len(part_1) >= len(part_2):
        if not "-" in part_1:
            return part_1
        else:
            return part_1.split("-")[0]
    else:
        if not "-" in part_2:
            return part_2
        else:
            return part_2.split("-")[0]

with open("/mnt/raid/users/johann/VCFs/old_sample_ids.txt", "r") as oldsample_ids_obj:
    old_line = oldsample_ids_obj.read()
    old_ids = "".join(old_line.split("\n")).split("\t", 9)[-1]

id_dict = {}
for old_id in old_ids.split("\t"):
    id_dict[extract_id(old_id)] = old_id

sample_ids = load_sample_excel()
create_ped_files()
print("done")