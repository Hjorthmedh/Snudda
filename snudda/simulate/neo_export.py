import h5py
import nixio

class NixExport:

    def __init__(self, snudda_file, nix_file: str):

        self.snudda_file = h5py.File(snudda_file, "r")
        self.nix_file = nixio.File.open(nix_file,
                                        nixio.FileMode.ReadWrite,
                                        compression=nixio.Compression.DeflateNormal)
        self.data_block = self.nix_file.blocks("Snudda Simulation")


    def convert(self):

        meta_section = self.nix_file.create_section("simulation_metadata", "nix.simulation_metadata")

        for data_name in self.snudda_file["meta_data"].keys():
            meta_section.create_property(data_name, self.snudda_file["meta_data"][data_name][()])


