

def extract_id(input):
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

new_ids = []
for old_id in old_ids.split("\t"):
    new_ids.append(extract_id(old_id))

new_line = "\t".join("".join(old_line.split("\n")).split("\t", 9)[:-1] + new_ids) + "\n"