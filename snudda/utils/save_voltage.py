import os.path

import h5py
import numpy as np
from mpi4py import MPI  # This must be imported before neuron, to run parallel
from neuron import h  # , gui


class SnuddaSaveVoltage:

    def __init__(self, voltage_file, network_data=None):

        self.voltage_file = voltage_file
        self.network_data = network_data

        self.pc = h.ParallelContext()

    def write(self, t_save, v_save, v_key, voltage_file=None):

        if not voltage_file:
            voltage_file = self.voltage_file

        self.pc.barrier()

        if int(self.pc.id()) == 0:

            if not os.path.isdir(os.path.dirname(voltage_file)):
                os.mkdir(os.path.dirname(voltage_file))

            print(f"Writing voltage output to {voltage_file}")
            output_file = h5py.File(voltage_file, "w")

            meta_data = output_file.create_group("metaData")
            voltage_data = output_file.create_group("voltData")
            spike_data = output_file.create_group("spikeData")

            if self.network_data:
                neuron_id = np.array([x["neuronID"] for x in self.network_data["neurons"]])
                meta_data.create_dataset("ID", data=neuron_id)

                neuron_names = [x["name"] for x in self.network_data["neurons"]]
                str_type = 'S' + str(max(1, max([len(x) for x in neuron_names])))
                meta_data.create_dataset("name", (len(neuron_names),), str_type, neuron_names, compression="gzip")

                neuron_types = [x["type"] for x in self.network_data["neurons"]]
                str_type = 'S' + str(max(1, max([len(x) for x in neuron_types])))
                meta_data.create_dataset("type", (len(neuron_names),), str_type, neuron_names, compression="gzip")

                # TODO: Add these later
                # meta_data.create_dataset("parameterFile")
                # meta_data.create_dataset("metaFile")
                # meta_data.create_dataset("neuronMorphology")
                # meta_data.create_dataset("parameterKey")
                # meta_data.create_dataset("morphologyKey")
                # meta_data.create_dataset("modulationKey")

            voltage_data.create_dataset("time", data=t_save*1e-3, compression="gzip")
            output_file.close()

        for i in range(int(self.pc.nhost())):

            if i == int(self.pc.id()):
                output_file = h5py.File(voltage_file, "a")

                for neuron_id, voltage in zip(v_key, v_save):
                    output_file["voltData"].create_dataset(str(neuron_id), data=voltage*1e-3, compression="gzip")

                output_file.close()

            self.pc.barrier()
