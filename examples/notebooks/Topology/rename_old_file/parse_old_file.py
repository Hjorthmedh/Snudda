import h5py
import csv

def collect_names(h5obj, path=""):
    rows = []

    for key in h5obj.keys():
        item = h5obj[key]
        full_path = f"{path}/{key}" if path else key

        # group or dataset
        obj_type = "group" if isinstance(item, h5py.Group) else "dataset"
        rows.append([obj_type, full_path, key, key])

        # attributes of this object
        for attr in item.attrs:
            rows.append(["attribute", full_path, attr, attr])

        # recurse into groups
        if isinstance(item, h5py.Group):
            rows.extend(collect_names(item, full_path))

    return rows


with h5py.File("../topology100/network-synapses.hdf5", "r") as f:
    rows = collect_names(f)

with open("rename_map.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["type", "path", "old_name", "new_name"])
    writer.writerows(rows)
