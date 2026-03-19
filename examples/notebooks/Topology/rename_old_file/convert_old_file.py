import h5py
import csv

def rename_from_csv(h5file, csv_path):
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # Sort by path depth descending so children rename before parents
    rows.sort(key=lambda r: r["path"].count("/"), reverse=True)

    for row in rows:
        if row["old_name"] == row["new_name"]:
            continue

        path = row["path"]
        parent_path = "/".join(path.split("/")[:-1])
        old_name = row["old_name"]
        new_name = row["new_name"]

        parent = h5file[parent_path] if parent_path else h5file

        if row["type"] in ("group", "dataset"):
            parent.move(old_name, new_name)

        elif row["type"] == "attribute":
            obj = h5file[path]
            obj.attrs[new_name] = obj.attrs[old_name]
            del obj.attrs[old_name]


with h5py.File("../topology100/network-synapses.hdf5", "r+") as f:
    rename_from_csv(f, "rename_map.csv")
