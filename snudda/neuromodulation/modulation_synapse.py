import json
import os


class NeuromodulationSynapse:

    def __init__(self):
        self.synapse_modulation = dict()
        self.weight = None
        self.type = 'adaptive'

    def set_weight(self, weight):
        self.weight = weight

    def set_connection_type(self, neuromodulation_key, connector):
        """
                neurotransmitter_key is equivalent to the parameter used in the mod files, which marks the
                level, mod and maxMod parameters e.g. neurotransmitter_key, ACh, would have parameters levelACh,
                modACh and maxModACh in modulated mod files.
        """

        if connector in self.synapse_modulation.keys():
            raise KeyError('connector is already defined')

        self.synapse_modulation.update({neuromodulation_key: {'connector': connector,
                                                              'cells': dict()}})

    def add_cell_modulation(self, neuromodulation_key, cell, ion_channels=None, receptors=None, extrinsic=None,
                            type_connection=None):

        if ion_channels is None:
            ion_channels = dict(soma=list(), dendrite=list(), axon=list())
        self.synapse_modulation[neuromodulation_key]['cells'].update({cell: {
            'ion_channels': ion_channels,
            'receptors': receptors,
            'extrinsic': extrinsic,
            'type': type_connection}})

    def save(self, dir_path, name):

        cell_types = list()

        for neurotransmitter_type, modulation_information in self.synapse_modulation.items():

            for c in modulation_information["cells"]:
                if c not in cell_types:
                    cell_types.append(c)

        for neurotransmitter_type, modulation_information in self.synapse_modulation.items():
            for c in cell_types:
                if c not in modulation_information["cells"]:
                    modulation_information["cells"].update({c: {"ion_channels": dict(soma=list(), dendrite=list(), axon=list()),
                                                                "receptors": None,
                                                                "extrinsic": None,
                                                                "type": None}})
        temp = dict()
        temp.update({'type': self.type})
        temp.update({'description': self.synapse_modulation})
        temp.update({'weight': self.weight})
        with open(os.path.join(dir_path, name), 'w') as out_file:
            json.dump(temp, out_file, indent=4)


if __name__ == "__main__":
    sw = NeuromodulationSynapse()
    sw.set_weight(weight=0)

    # Acetylcholine

    sw.set_connection_type(connector="concACh", neuromodulation_key="ACh")

    sw.add_cell_modulation(neuromodulation_key="ACh",
                           cell="dSPN",
                           ion_channels={
                               "soma": ["kir_ms", "cal12_ms", "cal13_ms", "can_ms", "Im_ms"],
                               "dendrite": ["kir_ms", "cal12_ms", "cal13_ms"]},
                           type_connection="spiking-concentration")

    sw.add_cell_modulation(neuromodulation_key="ACh",
                           cell="iSPN",
                           ion_channels={
                               "soma": ["kir_ms", "cal12_ms", "cal13_ms", "can_ms"],
                               "dendrite": ["kir_ms", "cal12_ms", "cal13_ms"]},
                           type_connection="spiking-concentration")

    # Dopamine

    sw.set_connection_type(connector="concDA", neuromodulation_key="DA")

    sw.add_cell_modulation(neuromodulation_key="DA",
                           cell="dSPN",
                           ion_channels={
                               "soma": ["kas_ms", "kaf_ms", "can_ms"],
                               "dendrite": ["kaf_ms", "kas_ms"],
                               "axon": []},
                           receptors={"tmGabaA": {"maxMod": 0.8},
                                      "tmGlut": {"maxMod_AMPA": 1.2,
                                                 "maxMod_NMDA": 1.3,
                                                 "failRate": 0.7}},
                           extrinsic=["CorticalBase", "CorticalSignal", "Thalamic"],
                           type_connection="spiking-concentration")

    sw.save(dir_path="", name="DA-ACh-control.json")
