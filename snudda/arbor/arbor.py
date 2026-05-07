#
# inspired by https://github.com/thorstenhater/cantata/blob/main/data/main.py
#

import arbor as A
from snudda import SnuddaLoad
import h5py

class recipe(A.recipe):
    def __init__(self, network_file, input_file, simulation_config):
        A.recipe.__init__(self)

        self.snudda_load = SnuddaLoad(network_file=network_file)

        self.input_data = h5py.File(input_file, "r")

    def cell_kind(self, gid):
        match self.snudda_load["neurons"][gid]["virtual_neuron"]:
            case True:
                return A.cell_kind.spike_source
            case False:
                return A.cell_kind.cable

    def num_cells(self):
        return self.snudda_load["num_neurons"]

    def connections_on(self, gid):

        synapses = self.snudda_load(post_id=gid)
        con_on = []

        for row in synapses:
            source = row[0]
            target = row[1]
            weight = row[11]
            delay = row[8] / self.axon_velocity + self.synapse_delay
            label = "detector"
            target = f"syn-{source}{row[2]}-{row[3]}-{row[4]}"

            con_on.append(
                A.connection((source, f"src-{label}"), target, weight, max(delay, self.dt))
            )

        return con_on


    def global_properties(self, kind):
        if kind == A.cell_kind.cable:
            return self.cable_props
        return None

    def cell_description(self, gid):
        kind = self.gid_to_kid[gid]
        if kind == 0:
            return self.make_cable_cell(gid)
        elif kind == 1:
            return self.make_lif_cell(gid)
        elif kind == 2:
            return self.make_vrt_cell(gid)
        else:
            raise RuntimeError("Unknown cell kind")

    def probes(self, gid):
        raise NotImplementedError("Not implemented yet")

    # https://github.com/thorstenhater/cantata/blob/main/data/main.py#L346
    def make_cable_cell(self, gid):
        raise NotImplementedError("Not implemented yet")

    def make_lif_cell(self, gid):
        raise NotImplementedError("Not implemented yet")

    def make_vrt_cell(self, gid):
        raise NotImplementedError("Not implemented yet")

    def load_cable_data(self, gid):
        raise NotImplementedError("Not implemented yet")


