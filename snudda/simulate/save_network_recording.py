import sys
import os.path
import h5py
import numpy as np
from mpi4py import MPI  # This must be imported before neuron, to run parallel
from neuron import h  # , gui

from snudda.neurons.NEURON_neuron_extractor import NEURONNeuronExtractor

"""

    Data hierarchy:
    
    SnuddaSaveNetworkRecordings - The Class which contains the data to be saved from the simulation
    
        Each Neuron has an instance of the NeuronRecordings class
        
            Within the NeuronRecordings class, the data can be saved as either:
            
                CompartmentData class, using self.register_compartment_data(**)
                    This contains all data within a neuron compartment (seg)
                
                SynapseData class, using self.register_synapse_data(**)
                    This contains all data within a synapse (point process)
                    
                SpikeData class, using self.register_spike_data(**)
                    This contains all spiking data (events) from a neuron
                    (the only data class which is not a time series)

"""

# TODO:

# Check if we can get synapse row number stored somehow...

# 4. Hur ska spikes sparas? (cell id, tid)
# 5. Gör en generell write funktion som tittar på vilken data som finns och anropar rätt sub-funktioner
# -- Skriv tutorial hur man lägger till en egen recording
#
# Förslag på namn: TimeSeriesData
#                  EventData

# NeuronActivity håller information om alla olika mätningar för en given neuron,
# Simulate har en lista med NeuronActivity
#
# I simulate måste vi:
# För varje neuron som ska mäta något på, skapa en NeuronActivity och lägg i en lista
# För varje measurement vi vill spara (ström, voltage, etc), använd add_measurements
# och sedan


class NeuronRecordings:
    """
        Container class for all recordings associated with a neuron.
    """

    def __init__(self, neuron_id):
        """
        Constructor

        Args:
            neuron_id  (int): Neuron ID
        """

        self.neuron_id = neuron_id

        self.data = dict()

    def register_compartment_data(self, data, data_type, sec_id, sec_x):
        """
        Adds a new measurement.
            
        Args:
            data (neuron.h.Vector): NEURON vector holding recording.
            data_type (str): Name of the tracked data.
            sec_id (int): Section ID
            sec_x (float): Section X (segment location)

        """

        if data_type not in self.data:
            self.data[data_type] = CompartmentData(neuron_id=self.neuron_id, data_type=data_type)

        self.data[data_type].append(data=data, sec_id=sec_id, sec_x=sec_x)

    def register_synapse_data(self, data, data_type, synapse_type, presynaptic_id, sec_id, sec_x, cond):

        if data_type not in self.data:
            self.data[data_type] = SynapseData(neuron_id=self.neuron_id, data_type=data_type)

        self.data[data_type].append(data=data, synapse_type=synapse_type, presynaptic_id=presynaptic_id,
                                    sec_id=sec_id, sec_x=sec_x, cond=cond)

    def register_spike_data(self, data, sec_id, sec_x):

        data_type = "spikes"
        if data_type not in self.data:
            self.data[data_type] = SpikeData(neuron_id=self.neuron_id)

        self.data[data_type].append(data=data, sec_id=sec_id, sec_x=sec_x)


class SynapseData:

    def __init__(self, neuron_id, data_type):
        self.neuron_id = neuron_id
        self.data = []
        self.data_type = data_type
        self.sec_id = []
        self.sec_x = []
        self.synapse_type = []
        self.presynaptic_id = []
        self.cond = []

    def append(self, data, synapse_type, presynaptic_id, sec_id, sec_x, cond):
        self.data.append(data)
        self.synapse_type.append(synapse_type)
        self.sec_id.append(sec_id)
        self.sec_x.append(sec_x)
        self.presynaptic_id.append(presynaptic_id)
        self.cond.append(cond)

    def to_numpy(self):
        """
            Returns:
                (np.ndarray): Data represented as np.ndarrays
        """
        return np.vstack([np.array(d) for d in self.data])


class CompartmentData:
    """
        Container class for recordings from a compartment. 
    """

    def __init__(self, neuron_id, data_type):
        """
        Constructor

        Args:
            neuron_id (int): Neuron ID
            data_type (str): Name of the tracked data.
        """
        self.neuron_id = neuron_id
        self.data_type = data_type
        self.data = []
        self.sec_id = []
        self.sec_x = []

    def append(self, data, sec_id, sec_x):
        # !!! Issue (?): The same compartment can hold several recordings now.

        """
        Appends a recording for a given compartment.

        Args:
            data (neuron.h.Vector): NEURON vector for the recording.
            sec_id (int): Section ID.
            sec_x (float): Section X (segment location). 

        """
        self.data.append(data)

        self.sec_id.append(sec_id)
        self.sec_x.append(sec_x)

    def to_numpy(self):  # Misnomer? (original data is preserved). (Instead "as_ndarray" (like in NEURON)?)
        """
            Returns:
                (np.ndarray): Data represented as np.ndarrays 
        """
        return np.vstack([np.array(d) if d.size() > 0 else np.array([]) for d in self.data])


class SpikeData(CompartmentData):

    def __init__(self, neuron_id):
        super().__init__(neuron_id, data_type="spikes")


class ExtracellularRecordings:

    def __init__(self, region_name, ecs):

        # Needs to identify with its index, and then we also track coordinates
        # we want to avoid duplicate data -- since a node has many coordinates inside it
        # EXAMPLE: https://neuron.yale.edu/neuron/docs/extracellular-diffusion

        self.ecs = ecs   # remove?
        self.region_name = region_name
        self.data = dict()
        self.coords = dict()

        self.xyz_min = np.array((self.ecs[self.region_name].compartments["ecs"]._xlo,
                                 self.ecs[self.region_name].compartments["ecs"]._ylo,
                                 self.ecs[self.region_name].compartments["ecs"]._zlo))

        self.xyz_max = np.array((self.ecs[self.region_name].compartments["ecs"]._xhi,
                                 self.ecs[self.region_name].compartments["ecs"]._yhi,
                                 self.ecs[self.region_name].compartments["ecs"]._zhi))

        self.n_voxels = np.array((self.ecs[self.region_name].compartments["ecs"]._nx,
                                  self.ecs[self.region_name].compartments["ecs"]._ny,
                                  self.ecs[self.region_name].compartments["ecs"]._nz))

    def get_coords(self, index_ijk):
        return np.vstack([self.coords[tuple(ijk)] for ijk in index_ijk])

    def register_extracellular_data(self, data, data_type, index_ijk):

        if data_type not in self.data:
            self.data[data_type] = ExtracellularNodeData(data_type=data_type)

        if index_ijk not in self.coords:
            delta_xyz = (self.xyz_max - self.xyz_min) / self.n_voxels
            self.coords[index_ijk] = self.xyz_min + delta_xyz * index_ijk

        self.data[data_type].append(data=data, index_ijk=index_ijk)

        # TODO: This object tracks the data in one "node",
        # The node knows its i,j,
        # It will have a dictionary, and there we need to know the species type, and the location
        # we will use i, j, k for index (but also store the centre of the node)
        # This node needs to know its own region, species, and location (x,y,z)
        #  We need to save the extracellular data, based on region, species, loc

# TODO: Check if RxD uses JIT, as they claim in Omar Awile et al, 2022.


class ExtracellularNodeData:

    def __init__(self, data_type):
        self.data_type = data_type
        self.data = dict()

    def append(self, data, index_ijk, unsafe=False):
        if not unsafe and index_ijk in self.data:
            raise ValueError(f"Index {index_ijk} already exists in the data structure.")

        self.data[index_ijk] = data

    def sort(self):
        sorted_items = sorted(self.data.items(), key=lambda x: x[0])
        self.data = dict(sorted_items)

    def to_numpy(self):
        return np.vstack([np.array(d) for d in self.data.values()])

    def get_data(self):
        self.sort()

        data = self.to_numpy()
        index = np.array(list(self.data.keys()))

        return index, data


class SnuddaSaveNetworkRecordings:

    # TODO: Add saving of simulation_config file (and experiment_config_file for pair recording)

    def __init__(self, output_file, network_data=None, sample_dt=None, node_id=0):

        # Only do this check on the first node
        if node_id == 0 and not self.check_file_available(output_file):
            sys.exit(-1)

        self.output_file = output_file
        self.network_data = network_data
        self.header_exists = False
        self.neuron_activities = dict()
        self.extracellular_data = dict()
        self.time = None
        self.sample_dt = sample_dt

        self.units = dict()

        self.pc = h.ParallelContext()

        self.save_geometry = False
        self.simulation = None

    def include_geometry(self, simulation):
        self.simulation = simulation
        self.save_geometry = True

    def set_new_output_file(self, output_file):
        self.output_file = output_file
        self.header_exists = False

    def add_unit(self, data_type, target_unit, conversion_factor):
        # Units reference chart: https://www.neuron.yale.edu/neuron/static/docs/units/unitchart.html
        self.units[data_type] = {"unit": target_unit, "conversion_factor": conversion_factor}

    def get_conversion(self, data_type):
        if data_type in self.units:
            return self.units[data_type]["conversion_factor"]
        else:
            return 1.0

    def register_compartment_data(self, data_type, neuron_id, data, sec_id, sec_x):
        if neuron_id not in self.neuron_activities:
            self.neuron_activities[neuron_id] = NeuronRecordings(neuron_id)

        self.neuron_activities[neuron_id].register_compartment_data(data=data, data_type=data_type,
                                                                    sec_id=sec_id, sec_x=sec_x)

    def register_synapse_data(self, data_type, neuron_id, data, synapse_type, presynaptic_id, sec_id, sec_x, cond):
        if neuron_id not in self.neuron_activities:
            self.neuron_activities[neuron_id] = NeuronRecordings(neuron_id)

        self.neuron_activities[neuron_id].register_synapse_data(data=data, data_type=data_type,
                                                                synapse_type=synapse_type,
                                                                presynaptic_id=presynaptic_id,
                                                                sec_id=sec_id, sec_x=sec_x, cond=cond)

    def register_spike_data(self, neuron_id, data, sec_id, sec_x):
        if neuron_id not in self.neuron_activities:
            self.neuron_activities[neuron_id] = NeuronRecordings(neuron_id)

        self.neuron_activities[neuron_id].register_spike_data(data=data, sec_id=sec_id, sec_x=sec_x)

    def register_extracellular_data(self, region_name, data_type, data, ecs, index_ijk):

        if region_name not in self.extracellular_data:
            self.extracellular_data[region_name] = ExtracellularRecordings(region_name=region_name, ecs=ecs)

        self.extracellular_data[region_name].register_extracellular_data(data=data, data_type=data_type,
                                                                         index_ijk=index_ijk)

    def register_time(self, time):
        self.time = time

    @staticmethod
    def spike_sort(t_spikes, id_spikes):

        spike_count = dict()

        for t, idx in zip(t_spikes, id_spikes):
            if idx in spike_count:
                spike_count[idx] += 1
            else:
                spike_count[idx] = 1

        spikes = dict()
        spike_ctr = dict()

        for idx, num in spike_count.items():
            spikes[idx] = np.full((num,), np.nan)  # Initialise to NaN to catch any unassigned values
            spike_ctr[idx] = 0

        for t, idx in zip(t_spikes, id_spikes):
            spikes[idx][spike_ctr[idx]] = t
            spike_ctr[idx] += 1

        # Internal consistency
        for idx in spikes:
            assert len(spikes[idx]) == 0 or not np.isnan(spikes[idx][-1])
            assert (idx in spike_ctr and spikes[idx].shape[0] == spike_ctr[idx]) \
                   or (idx not in spike_ctr and len(spikes[idx]) == 0)

        return spikes

    def write_string_meta_data(self, group, name):
        string_data = [x[name] if x[name] is not None else "" for x in self.network_data["neurons"]]
        max_len = max(1, max([len(x) for x in string_data]))
        str_type = f"S{max_len}"
        group.create_dataset(name, (len(string_data),), str_type, string_data, compression="gzip")

    def check_file_available(self, file_name):

        if os.path.isfile(file_name):
            # Try to open and close file, to make sure it is available
            try:
                f = h5py.File(file_name, "w")
                f.close()
            except BlockingIOError as e:
                print(f"Unable to create file {file_name}. Is some other program using the file?")
                return False

        return True

    def write_header(self):

        self.pc.barrier()

        if self.header_exists:
            return

        if int(self.pc.id()) == 0:

            if not os.path.isdir(os.path.dirname(self.output_file)):
                os.mkdir(os.path.dirname(self.output_file))

            print(f"Writing network output to {self.output_file}")

            try:
                out_file = h5py.File(self.output_file, "w")
            except Exception as e:

                print(e)
                print(f"Trying to recover, and save to different file name.")

                ctr = 1
                temp_name = f"{self.output_file}-{ctr}"

                while os.path.isfile(temp_name):
                    print(f"File exists: {temp_name}")
                    ctr += 1
                    temp_name = f"{self.output_file}-{ctr}"

                print(f"\n!!! Unable to create {self.output_file} (file locked?), using {temp_name} instead\n\n")
                self.output_file = temp_name

                out_file = h5py.File(temp_name, "w")

            meta_data = out_file.create_group("meta_data")

            meta_data.create_dataset("network_file", data=self.network_data["network_file"])
            out_file.create_group("neurons")

            if self.network_data:
                neuron_id = np.array([x["neuron_id"] for x in self.network_data["neurons"]])
                meta_data.create_dataset("id", data=neuron_id)

                for name in ["name", "type", "morphology", "parameter_key", "morphology_key", "modulation_key"]:
                    self.write_string_meta_data(group=meta_data, name=name)

                meta_data.create_dataset("population_unit", data=self.network_data["population_unit"], compression="gzip")
                meta_data.create_dataset("position", data=self.network_data["neuron_positions"], compression="gzip")
                meta_data.create_dataset("virtual_neuron", compression="gzip",
                                         data=[n["virtual_neuron"] for n in self.network_data["neurons"]])

                # TODO: Also write extraAxons info to file

            out_file.close()

        self.header_exists = True

        self.pc.barrier()

    def get_sample_step(self):

        if self.sample_dt is None:
            return None
        else:
            converted_time = np.array(self.time) * self.get_conversion("time")
            dt = converted_time[1] - converted_time[0]
            sample_step = int(np.round(self.sample_dt / dt))
            return sample_step

    def write(self):

        self.write_header()

        sample_step = self.get_sample_step()

        if int(self.pc.id()) == 0:

            # TODO: We need to save max time even if we do not save the soma voltage(!)
            out_file = h5py.File(self.output_file, "a")

            if "time" not in out_file and self.time is not None:
                print(f"Using sample dt = {self.sample_dt} (sample step size {sample_step})")
                try:
                    out_file.create_dataset("time", data=np.array(self.time)[::sample_step] * self.get_conversion("time"))
                    out_file.close()
                except:
                    import traceback
                    print(f"self.time={self.time}")
                    print(f"{traceback.format_exc()}", flush=True)
                    import pdb
                    pdb.set_trace()

        for i in range(int(self.pc.nhost())):

            self.pc.barrier()

            if int(self.pc.id()) == i:
                print(f"Worker {i+1}/{int(self.pc.nhost())} writing data to {self.output_file}")

                out_file = h5py.File(self.output_file, "a")

                # Make sure time data is written to file, in case node 0 only has virtual neurons
                if "time" not in out_file and self.time is not None:
                    print(f"Using sample dt = {self.sample_dt} (sample step size {sample_step})")
                    try:
                        out_file.create_dataset("time",
                                                data=np.array(self.time)[::sample_step] * self.get_conversion("time"))
                        out_file.close()
                    except:
                        import traceback
                        print(f"self.time={self.time}")
                        print(f"{traceback.format_exc()}", flush=True)
                        import pdb
                        pdb.set_trace()

                for na in self.neuron_activities.values():
                    neuron_id_str = str(na.neuron_id)

                    if neuron_id_str not in out_file["neurons"]:
                        out_file["neurons"].create_group(neuron_id_str)

                    for m in na.data.values():

                        conversion_factor = self.get_conversion(m.data_type)

                        if isinstance(m, SpikeData):
                            # Spike data is not a time series, and should never be downsampled
                            data_sample_step = None
                        else:
                            data_sample_step = sample_step

                        data_set = out_file["neurons"][neuron_id_str].create_dataset(m.data_type,
                                                                                     data=m.to_numpy()[:, ::data_sample_step]*conversion_factor,
                                                                                     compression="gzip")

                        data_set.attrs["sec_id"] = np.array(m.sec_id)
                        data_set.attrs["sec_x"] = np.array(m.sec_x)

                        if isinstance(m, SynapseData):
                            data_set.attrs["synapse_type"] = np.array(m.synapse_type)
                            data_set.attrs["presynaptic_id"] = np.array(m.presynaptic_id)
                            data_set.attrs["cond"] = np.array(m.cond)

                if int(self.pc.id()) == 0:
                    # RxD computes extracellular diffusion INDEPENDENTLY on ALL nodes (sigh!)
                    # thus we need to save only data from node 0

                    ecs_group = out_file.create_group("ecs")

                    for region_name, region_data in self.extracellular_data.items():
                        region_group = ecs_group.create_group(region_name)

                        for data_type, extracellular_data in region_data.data.items():
                            ecs_data_group = region_group.create_group(data_type)

                            index, data = extracellular_data.get_data()
                            coords = region_data.get_coords(index)

                            ecs_data_group.create_dataset("index", data=index)
                            ecs_data_group.create_dataset("xyz", data=coords)
                            ecs_data_group.create_dataset("data", data=data)

                if self.save_geometry:
                    nne = NEURONNeuronExtractor(sim=self.simulation)
                    nne.write_all_coordinates_to_hdf5(hdf5_file=out_file)

                out_file.close()

        self.pc.barrier()

