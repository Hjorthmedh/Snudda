#!/usr/bin/env python3

# A wrapper script for the touch detection algorithm
#
# Usage:
#
# snudda init <networkPath> --size XXX
# -- Creates an a json config file
#
# snudda place <networkPath> [--parallel]
# -- Cell placement within volumes specified
#
# snudda detect <networkPath> [--hvsize hyperVoxelSize] [--parallel]
# -- Touch detection of putative synapses
#
# snudda prune <networkPath> [--mergeonly] [--parallel]
# -- Prune the synapses
#
# snudda input <networkPath> [--input yourInputConfig] [--parallel]
#
# snudda export <networkPath>
# -- Export to SONATA format (optional)
#
# snudda simulate <networkPath>
#
# snudda analyse <networkPath>
#
#
# snudda help me

# Johannes Hjorth, Royal Institute of Technology (KTH)
# Human Brain Project 2019

#
# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union's Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
# and SGA2).
#


import os
import sys
import timeit
from collections import OrderedDict

# import pkg_resources
from importlib import resources

import json
import numpy as np

from snudda.utils import snudda_path
from snudda.utils.snudda_path import snudda_isfile, get_snudda_data


def get_data_file(*dirs):
    path = os.path.join("data", *dirs)

    try:
        return resources.files(__package__).joinpath(path).as_posix()
    except FileNotFoundError:
        raise FileNotFoundError(f"Data file '{path}' not found")

# def get_data_file(*dirs):
#     path = os.path.join("data", *dirs)
#     if not pkg_resources.resource_exists(__package__, path):
#         raise FileNotFoundError("Data file '{}' not found".format(path))
#     return pkg_resources.resource_filename(__package__, path)


class Snudda(object):

    """ Wrapper class, calls Snudda helper functions """

    def __init__(self, network_path, parallel=False, ipython_profile=None):

        """
        Instantiates Snudda
        :param network_path: Location of Snudda network
        """

        self.network_path = network_path

        self.d_view = None
        self.rc = None
        self.slurm_id = 0

        self.parallel = parallel
        self.ipython_profile = ipython_profile

        # Add current dir to python path
        sys.path.append(os.getcwd())

        self.start = timeit.default_timer()

    ############################################################################

    @staticmethod
    def help_info(args):
        """ Prints Snudda help """
        from snudda.help import snudda_help_text

    ############################################################################

    def init_config_wrapper(self, args):

        """
        Creates network-config.json in network_path.

        Args:
            args : command line arguments from argparse

        Example:
             snudda init -size 100 [-overwrite] [-randomseed 1234] [--profile] [--verbose] path
        """

        assert args.size is not None, "You need to specify --size when initialising config for the network"

        self.init_config(network_size=args.size,
                         snudda_data=args.snudda_data,
                         neurons_dir=args.neurons_dir,
                         connection_file=args.connection_file,
                         overwrite=args.overwrite,
                         random_seed=args.randomseed,
                         honor_stay_inside=args.stay_inside)

    def init_config(self,
                    network_size=None,
                    snudda_data=None,
                    struct_def=None,
                    neurons_dir=None,
                    connection_file=None,
                    honor_stay_inside=True,   # currently the cli.py defaults to sending False
                    overwrite=False,
                    random_seed=None):

        print(f"Legacy config creation.")

        # self.networkPath = args.path
        print("Creating config file")
        print(f"Network path: {self.network_path}")

        from snudda.init.init import SnuddaInit

        if struct_def is None:
            struct_def = {"Striatum": network_size,
                          "GPe": 0,
                          "GPi": 0,
                          "SNr": 0,
                          "STN": 0,
                          "Cortex": 0,  # Cortex and thalamus axons disabled right now, set to N > 0 to include a few
                          "Thalamus": 0}

        if not overwrite:
            assert not os.path.exists(self.network_path), \
                (f"Network path {self.network_path} already exists (aborting to prevent accidental overwriting)."
                 "\nCall snudda init with --overwrite to override and overwrite the old data.")

        self.make_dir_if_needed(self.network_path)

        config_file = os.path.join(self.network_path, "network-config.json")
        SnuddaInit(network_path=self.network_path,
                   struct_def=struct_def,
                   neurons_dir=neurons_dir,
                   snudda_data=snudda_data,
                   config_file=config_file,
                   honor_stay_inside=honor_stay_inside,
                   random_seed=random_seed,
                   connection_override_file=connection_file)

        if network_size is not None and network_size > 1e5:
            print(f"Make sure there is enough disk space in {self.network_path}")
            print("Large networks take up ALOT of space")

    ############################################################################

    def init_tiny(self, neuron_paths, neuron_names, number_of_neurons,
                  snudda_data=None,
                  morphology_key=None, parameter_key=None,
                  connection_config=None, random_seed=None, density=80500, d_min=15e-6):

        """
            network_path : Network path

        """

        from snudda.init.init import SnuddaInit
        from snudda.place import create_cube_mesh

        n_total = np.sum(number_of_neurons)

        si = SnuddaInit(network_path=self.network_path,
                        snudda_data=snudda_data,
                        random_seed=random_seed)

        si.define_structure(struct_name="Cube",
                            struct_mesh="cube",
                            d_min=d_min,
                            struct_centre=(0.0, 0.0, 0.0),
                            side_len=(n_total/density)**(1/3)*1e-3,
                            num_neurons=n_total,
                            n_putative_points=n_total*5)

        if connection_config is not None:
            si.replace_connectivity(connection_file=connection_config)

        if isinstance(neuron_paths, str):
            neuron_paths = [neuron_paths]

        if isinstance(neuron_names, str):
            neuron_names = [neuron_names]

        if isinstance(number_of_neurons, int):
            number_of_neurons = [int(number_of_neurons / len(neuron_paths)) for x in neuron_names]

        if isinstance(morphology_key, str):
            morphology_key = [morphology_key]

        if isinstance(parameter_key, str):
            parameter_key = [parameter_key]

        assert (morphology_key is None and parameter_key is None) or \
            len(morphology_key) == len(parameter_key) == len(neuron_paths)

        for idx, (path, name, cnt) in enumerate(zip(neuron_paths, neuron_names, number_of_neurons)):
            if name in si.network_data["regions"]["Cube"]["neurons"]:
                raise ValueError(f"neuron name {name} defined more than once")

            si.add_neurons(name=name, neuron_dir=path, region_name="Cube", num_neurons=cnt)

            if morphology_key is not None:
                si.network_data["regions"]["Cube"]["neurons"][name]["parameter_key"] = parameter_key[idx]
                si.network_data["regions"]["Cube"]["neurons"][name]["morphology_key"] = morphology_key[idx]

        si.write_json()

        return si

    ############################################################################

    def import_config_wrapper(self, args):

        self.import_config(network_config_file=args.config_file,
                           snudda_data=args.snudda_data,
                           overwrite=args.overwrite)

    def import_config(self, network_config_file, snudda_data=None, overwrite=False):

        from snudda.init.init_config import ConfigParser

        conf = ConfigParser(config_file=network_config_file, snudda_data=snudda_data)
        conf.parse_config()

        if not os.path.isdir(self.network_path):
            print(f"Creating directory {self.network_path}")
            os.makedirs(self.network_path)

        new_config_file = os.path.join(self.network_path, "network-config.json")

        if os.path.isfile(new_config_file) and not overwrite:
            print(f"\n!!! ERROR: File already exists: {new_config_file}\nSet 'overwrite' to overwrite old file.\n")
            exit(-1)

        conf.replace_network_path(network_path=os.path.abspath(self.network_path))
        conf.write_config(new_config_file)

    ############################################################################

    def create_network(self, honor_morphology_stay_inside=True):

        # This is a helper function, to create the full network
        self.place_neurons(honor_morphology_stay_inside=honor_morphology_stay_inside)
        self.detect_synapses()
        self.prune_synapses()

    ############################################################################

    def place_neurons_wrapper(self, args):
        """
        Places neurons in 3D space. Creates network-neuron-positions.hdf5 in network_path.

        Args:
            args : command line arguments from argparse

        Example:
            snudda place [--profile] [--verbose] [--h5legacy] [-parallel] path

        """

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        self.place_neurons(random_seed=args.randomseed,
                           parallel=args.parallel,
                           ipython_profile=args.ipython_profile,
                           ipython_timeout=args.ipython_timeout,
                           h5libver=h5libver,
                           verbose=args.verbose,
                           honor_morphology_stay_inside=args.stay_inside)

    def place_neurons(self,
                      random_seed=None,
                      parallel=None,
                      ipython_profile=None,
                      ipython_timeout=120,
                      h5libver="latest",
                      verbose=False,
                      honor_morphology_stay_inside=True):

        if parallel is None:
            parallel = self.parallel

        if ipython_profile is None:
            ipython_profile = self.ipython_profile

        # self.networkPath = args.path
        print("Placing neurons")
        print(f"Network path: {self.network_path}")

        log_file_name = os.path.join(self.network_path, "log", "place-neurons.txt")

        self.setup_log_file(log_file_name)  # sets self.logFile

        if parallel:
            self.setup_parallel(ipython_profile=ipython_profile, timeout=ipython_timeout)  # sets self.d_view

        from snudda.place.place import SnuddaPlace

        sp = SnuddaPlace(network_path=self.network_path,
                         log_file=self.logfile,
                         verbose=verbose,
                         d_view=self.d_view,
                         h5libver=h5libver,
                         random_seed=random_seed,
                         morphologies_stay_inside=honor_morphology_stay_inside)

        sp.place()

        self.cleanup_workers()

        self.stop_parallel()
        self.close_log_file()

        return sp

    ############################################################################

    def detect_synapses_wrapper(self, args):
        """
        Synapse touch detection. Writes results to network_path/voxels (one file per hypervoxel).
        Also adds synapse projections between structures.
        Results written to network-projection-synapses.hdf5 in network_path.

        Args:
            args : command line arguments from argparse

        Example:
            snudda detect [-cont] [-hvsize HVSIZE] [--volumeID VOLUMEID] [--profile] [--verbose] [--h5legacy] [-parallel] path

        """

        if args.hvsize is not None:
            hyper_voxel_size = int(args.hvsize)
        else:
            hyper_voxel_size = 100

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        self.detect_synapses(random_seed=args.randomseed,
                             parallel=args.parallel,
                             ipython_profile=args.ipython_profile,
                             ipython_timeout=args.ipython_timeout,
                             hyper_voxel_size=hyper_voxel_size,
                             volume_id=args.volumeID,
                             h5libver=h5libver,
                             verbose=args.verbose,
                             cont=args.cont)

    def detect_synapses(self,
                        random_seed=None,
                        parallel=None,
                        ipython_profile=None,
                        ipython_timeout=120,
                        hyper_voxel_size=100,
                        volume_id=None,
                        h5libver="latest",
                        verbose=False,
                        cont=False):

        if parallel is None:
            parallel = self.parallel

        if ipython_profile is None:
            ipython_profile = self.ipython_profile

        # self.networkPath = args.path
        print("Touch detection")
        print(f"Network path: {self.network_path}")

        log_dir = os.path.join(self.network_path, "log")
        if not os.path.exists(log_dir):
            print(f"Creating directory {log_dir}")
            os.makedirs(log_dir, exist_ok=True)

        config_file = os.path.join(self.network_path, "network-config.json")
        position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        log_filename = os.path.join(self.network_path, "log", "touch-detection.txt")
        save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")

        voxel_dir = os.path.join(self.network_path, "voxels")
        self.make_dir_if_needed(voxel_dir)

        self.setup_log_file(log_filename)  # sets self.logfile

        if parallel:
            self.setup_parallel(ipython_profile=ipython_profile, timeout=ipython_timeout)  # sets self.d_view

        from snudda.detect.detect import SnuddaDetect

        # You can now setup SnuddaDetect with only network_path and it will use default values
        # for config_file, position_file, logfile, save_file
        sd = SnuddaDetect(config_file=config_file,
                          position_file=position_file,
                          logfile=self.logfile,
                          save_file=save_file,
                          slurm_id=self.slurm_id,
                          volume_id=volume_id,
                          rc=self.rc,
                          hyper_voxel_size=hyper_voxel_size,
                          h5libver=h5libver,
                          random_seed=random_seed,
                          verbose=verbose)

        if cont:
            # Continue previous run
            print("Continuing previous touch detection")
            sd.detect(restart_detection_flag=False)
        else:
            sd.detect(restart_detection_flag=True)

        # Also run SnuddaProject to handle projections between volume

        from snudda.detect.project import SnuddaProject

        sp = SnuddaProject(network_path=self.network_path)
        sp.project()

        self.cleanup_workers()

        self.stop_parallel()
        self.close_log_file()

        return sd, sp

    ############################################################################

    def prune_synapses_wrapper(self, args):

        """
        Merges data from synapse detection, then does synapse pruning. Writes results to network-synapses.hdf5 in network_path.

        Args:
            args : command line arguments from argparse

        Example:
            snudda prune [--configFile CONFIG_FILE] [--profile] [--verbose] [--h5legacy] [--keepfiles] [-parallel] path
        """

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        self.prune_synapses(config_file=args.config_file,
                            random_seed=args.randomseed,
                            parallel=args.parallel,
                            ipython_profile=args.ipython_profile,
                            ipython_timeout=args.ipython_timeout,
                            verbose=args.verbose,
                            keep_files=args.keepfiles,
                            save_putative_synapses = args.savePutative)

    def prune_synapses(self,
                       config_file=None,
                       random_seed=None,
                       parallel=None,
                       ipython_profile=None,
                       ipython_timeout=120,
                       h5libver="latest",
                       verbose=False,
                       keep_files=False,
                       save_putative_synapses=False):

        if parallel is None:
            parallel = self.parallel

        if ipython_profile is None:
            ipython_profile = self.ipython_profile

        # self.networkPath = args.path
        print("Prune synapses")
        print(f"Network path: {self.network_path}")

        from snudda.detect.prune import SnuddaPrune

        log_filename = os.path.join(self.network_path, "log", "synapse-pruning.txt")

        self.setup_log_file(log_filename)  # sets self.logfile

        if parallel:
            self.setup_parallel(ipython_profile=ipython_profile, timeout=ipython_timeout)  # sets self.d_view

        # Optionally set this
        scratch_path = None

        sp = SnuddaPrune(network_path=self.network_path,
                         logfile=self.logfile,
                         logfile_name=log_filename,
                         config_file=config_file,
                         d_view=self.d_view,
                         scratch_path=scratch_path,
                         h5libver=h5libver,
                         random_seed=random_seed,
                         verbose=verbose,
                         keep_files=keep_files or save_putative_synapses)

        sp.prune()

        if save_putative_synapses:
            sp.save_putative_synapses()

        self.cleanup_workers()

        self.stop_parallel()
        self.close_log_file()

        return sp

    ############################################################################

    def setup_input_wrapper(self, args):
        """
        Creates synaptic input for network based on input.json, writes input-spikes.hdf5 in network_path.
        Args:
            args : command line arguments from argparse

        Example:
            snudda input [--input INPUT] [--inputFile INPUT_FILE] [--networkFile NETWORK_FILE] [--time TIME] [-randomseed 123] [--profile] [--verbose] [--h5legacy] [-parallel] path
        """

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        self.setup_input(network_file=args.network_file,
                         input_file=args.input_file,
                         input_config=args.input,
                         input_time=args.time,
                         random_seed=args.randomseed,
                         use_meta_input=not args.no_meta_input,
                         h5libver=h5libver,
                         parallel=args.parallel,
                         ipython_profile=args.ipython_profile,
                         ipython_timeout=args.ipython_timeout,
                         verbose=args.verbose)

    def setup_input(self,
                    network_file=None,
                    input_file=None,
                    input_config=None,
                    input_time=None,
                    use_meta_input=True,
                    random_seed=None,
                    h5libver="latest",
                    parallel=None,
                    ipython_profile=None,
                    ipython_timeout=120,
                    verbose=False):

        if parallel is None:
            parallel = self.parallel

        if ipython_profile is None:
            ipython_profile = self.ipython_profile

        print("Setting up inputs, assuming input.json exists")
        log_filename = os.path.join(self.network_path, "log", "setup-input.txt")
        self.setup_log_file(log_filename)  # sets self.logfile

        if parallel:
            self.setup_parallel(ipython_profile=ipython_profile, timeout=ipython_timeout)  # sets self.d_view

        from snudda.input.input import SnuddaInput

        if input_config is None:
            input_config = os.path.join(self.network_path, "input.json")

        snudda_data = get_snudda_data(network_path=self.network_path)
        if not snudda_isfile(input_config, snudda_data=snudda_data):
            print(f"Missing input config file: {input_config}")
            return

        if network_file is None:
            network_file = os.path.join(self.network_path, "network-synapses.hdf5")

        if input_file is None:
            input_file = os.path.join(self.network_path, "input-spikes.hdf5")

        print(f"Writing input spikes to {input_file}")

        si = SnuddaInput(input_config_file=input_config,
                         hdf5_network_file=network_file,
                         spike_data_filename=input_file,
                         time=input_time,
                         logfile=self.logfile,
                         rc=self.rc,
                         random_seed=random_seed,
                         h5libver=h5libver,
                         verbose=verbose,
                         use_meta_input=use_meta_input)
        si.generate()

        self.cleanup_workers()

        self.stop_parallel()
        self.close_log_file()

        return si

    ############################################################################

    def export_to_SONATA_wrapper(self, args):
        """
        Export network to SONATA files. (Currently not functioning)

        Args:
            args : command line arguments from argparse
        """

        self.export_to_SONATA(network_file=args.network_file,
                              input_file=args.input_file)

    def export_to_SONATA(self, network_file=None, input_file=None):

        from snudda.utils.export_sonata import ExportSonata

        print("Exporting to SONATA format")
        print(f"Network path: {self.network_path}")

        if network_file is None:
            network_file = os.path.join(self.network_path, "network-synapses.hdf5")

        if input_file is None:
            input_file = os.path.join(self.network_path, "input-spikes.hdf5")

        out_dir = os.path.join(self.network_path, "SONATA")

        cn = ExportSonata(network_file=network_file, input_file=input_file, out_dir=out_dir)

    ############################################################################

    @staticmethod
    def compile_mechanisms(mech_dir=None, snudda_data=None):

        if not mech_dir:
            mech_dir = os.path.realpath(snudda_path.snudda_parse_path(os.path.join("$DATA", "neurons", "mechanisms"),
                                                                      snudda_data=snudda_data))

        if not os.path.exists("x86_64") and not os.path.exists("nrnmech.dll")\
                and not os.path.exists("aarch64") and not os.path.exists("arm64"):

            from mpi4py import MPI  # This must be imported before neuron, to run parallel
            from neuron import h
            pc = h.ParallelContext()

            if pc.id() == 0:
                # Only run this on master node
                print(f"Running on master node:  nrnivmodl {mech_dir}")
                os.system(f"nrnivmodl {mech_dir}")
            else:
                print("Worker waiting for master node to compile NEURON modules.")

            pc.barrier()

            try:
                if os.path.exists("nrnmech.dll"):
                    mech_path = "nrnmech.dll"
                elif os.path.exists("x86_64"):
                    mech_path = "x86_64/.libs/libnrnmech.so"
                elif os.path.exists("aarch64"):
                    mech_path = "aarch64/.libs/libnrnmech.so"
                elif os.path.exists("arm64"):
                    mech_path = "arm64/.libs/libnrnmech.so"
                else:
                    print(f"Could not find compiled mechanisms. Compile using 'nrnivmodl {mech_dir}' "
                          f"and retry simulation.")
                    sys.exit(-1)

                print(f"Loading mechanisms from '{os.path.abspath(mech_path)}'")
                h.nrn_load_dll(mech_path)

            except:
                import traceback
                print(f"Error while loading mechanisms:\n{traceback.format_exc()}")


        else:
            print("NEURON mechanisms already compiled, make sure you have the correct version of NEURON modules."
                  "\nIf you delete x86_64, aarch64, arm64 directories (or nrnmech.dll) "
                  "then you will force a recompilation of the modules.")

    ############################################################################

    def simulate_wrapper(self, args):
        """
        Simulate network. Writes results to network_path/simulation.

        Args:
            args : command line arguments from argparse

        Example:
            snudda simulate [--networkFile NETWORK_FILE] [--inputFile INPUT_FILE] [--time TIME]
            [--spikesOut SPIKES_OUT] [--noVolt] [--disableGJ]
            [-mechdir MECH_DIR] [--profile] [--verbose] [--exportCoreNeuron] path
        """

        print(f"args: {args}")

        assert args.enable_rxd_neuromodulation is None or args.disable_rxd_neuromodulation is None, \
            "You can only specify enable_rxd_neuromodulation or disable_rxd_neuromodulation, not both"

        use_rxd_neuromodulation = None

        if args.enable_rxd_neuromodulation == True:
            use_rxd_neuromodulation = True

        if args.disable_rxd_neuromodulation == True:
            use_rxd_neuromodulation = False

        sim = self.simulate(network_file=args.network_file, input_file=args.input_file,
                            output_file=args.output_file, snudda_data=args.snudda_data,
                            time=args.time,
                            mech_dir=args.mech_dir,
                            # neuromodulation=args.neuromodulation,
                            disable_synapses=args.disable_synapses,
                            disable_gj=args.disable_gj,
                            record_volt=args.record_volt,
                            record_all=args.record_all,
                            simulation_config=args.simulation_config,
                            export_core_neuron=args.exportCoreNeuron,
                            use_rxd_neuromodulation=use_rxd_neuromodulation,
                            verbose=args.verbose)

        sim.clear_neuron()

    def simulate(self,
                 network_file=None,
                 input_file=None,
                 output_file=None,
                 snudda_data=None,
                 time=None,
                 mech_dir=None,
                 # neuromodulation=None,
                 disable_synapses=None,
                 disable_gj=False,
                 record_volt=True,
                 record_all=False,
                 sample_dt=None,
                 simulation_config=None,
                 export_core_neuron=False,
                 use_rxd_neuromodulation=None,
                 verbose=False):

        start = timeit.default_timer()

        from mpi4py import MPI  # This must be imported before neuron, to run parallel

        # Initialize MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        print(f"MPI Rank: {rank}, Size: {size}")

        from neuron import h
        pc = h.ParallelContext()

        if network_file is None:
            network_file = os.path.join(self.network_path, "network-synapses.hdf5")

        if input_file is None:
            putative_input_file = os.path.join(self.network_path, "input-spikes.hdf5")
            if os.path.isfile(putative_input_file):
                input_file = putative_input_file

        if output_file is None:
            output_file = os.path.join(self.network_path, "simulation", "output.hdf5")

        self.make_dir_if_needed(os.path.join(self.network_path, "simulation"))

        print(f"Using input file {input_file}")

        # nWorkers = args.ncores
        # print("Using " + str(nWorkers) + " workers for neuron")

        # Problems with nested symbolic links when the second one is a relative
        # path going beyond the original base path

        if mech_dir is None:
            # Take into account which SNUDDA_DATA the user wants to use
            from snudda.utils.snudda_path import get_snudda_data
            snudda_data = get_snudda_data(snudda_data=snudda_data,
                                          network_path=self.network_path)

            mech_dir = os.path.realpath(snudda_path.snudda_parse_path(os.path.join("$DATA", "neurons", "mechanisms"),
                                                                      snudda_data))

        self.compile_mechanisms(mech_dir=mech_dir)

        save_dir = os.path.join(os.path.dirname(network_file), "simulation")

        if not os.path.exists(save_dir):
            print(f"Creating directory {save_dir}")
            os.makedirs(save_dir, exist_ok=True)

        # Get the SlurmID, used in default file names
        slurm_id = os.getenv('SLURM_JOBID')

        if slurm_id is None:
            slurm_id = str(666)

        if disable_gj:
            print("!!! WE HAVE DISABLED GAP JUNCTIONS !!!")

        if disable_synapses:
            print("!!! SYNAPSES DISABLED")

        log_file = os.path.join(os.path.dirname(network_file), "log", "network-simulation-log.txt")

        log_dir = os.path.join(os.path.dirname(network_file), "log")
        if not os.path.exists(log_dir):
            print(f"Creating directory {log_dir}")
            os.makedirs(log_dir, exist_ok=True)

        from snudda.simulate.simulate import SnuddaSimulate

        # Simulate is deterministic, no random seed.
        sim = SnuddaSimulate(network_file=network_file,
                             input_file=input_file,
                             output_file=output_file,
                             snudda_data=snudda_data,
                             disable_gap_junctions=disable_gj,
                             disable_synapses=disable_synapses,
                             log_file=log_file,
                             simulation_config=simulation_config,
                             sample_dt=sample_dt,
                             use_rxd_neuromodulation=use_rxd_neuromodulation,
                             verbose=verbose)
        sim.setup()
        sim.add_external_input()

        sim.check_memory_status()

        if record_volt:
            # sim.add_volt_recording_all()
            # Either use record_volt or specify "record_all_soma" in simulation_config file
            sim.add_volt_recording_soma()
            # sim.addRecordingOfType("dSPN",5) # Side len let you record from a subset

        if record_all:
            record_cell_id = np.array([int(x) for x in record_all.split(",")])
            sim.add_volt_recording_all(cell_id=record_cell_id)
            sim.add_synapse_current_recording_all(record_cell_id)

        if time is not None:
            t_sim = time * 1000  # Convert from s to ms for Neuron simulator
        else:
            # Will attempt to read "time" from simulation_config
            t_sim = None

        if export_core_neuron:
            sim.export_to_core_neuron()
            return  # We do not run simulation when exporting to core neuron

        sim.check_memory_status()
        if t_sim is None or t_sim > 0:
            print(f"Running simulation for {t_sim} ms.")
            sim.run(t_sim)  # In milliseconds

            print("Simulation done, saving output")
            sim.write_output()
        else:
            print(f"Time set to {t_sim} ms. No simulation run.")

        stop = timeit.default_timer()
        if sim.pc.id() == 0:
            print(f"Program run time: {stop - start:.1f}s")

        # OBS! You want to do sim.clear_neuron() after the simulation if you need
        #      to setup another neuron simulation afterwards.

        # sim.plot()
        return sim

    ############################################################################

    def analyse(self, args):

        print("Add analysis code here, see Network_analyse.py")

    ############################################################################

    def setup_parallel(self, ipython_profile=None, timeout=120):
        """Setup ipyparallel workers."""

        self.slurm_id = os.getenv('SLURM_JOBID')

        if self.slurm_id is None:
            self.slurm_id = 0
        else:
            self.slurm_id = int(self.slurm_id)

        self.logfile.write(f"Using slurm_id: {self.slurm_id}")

        if ipython_profile is None:
            ipython_profile = os.getenv('IPYTHON_PROFILE')

        if not ipython_profile:
            ipython_profile = "default"

        ipython_dir = os.getenv('IPYTHONDIR')
        if not ipython_dir:
            ipython_dir = os.path.join(os.path.abspath(os.getcwd()), ".ipython")

        self.logfile.write('Creating ipyparallel client\n')

        from ipyparallel import Client

        u_file = os.path.join(ipython_dir, f"profile_{ipython_profile}", "security", "ipcontroller-client.json")
        print(f"Reading IPYPARALLEL connection info from {u_file}\n")
        self.logfile.write(f"Reading IPYPARALLEL connection info from {u_file}\n")
        # self.rc = Client(profile=ipython_profile, url_file=u_file, timeout=120, debug=False)
        self.rc = Client(profile=ipython_profile, connection_info=u_file, timeout=timeout, debug=False)

        self.logfile.write(f'Client IDs: {self.rc.ids}')

        # http://davidmasad.com/blog/simulation-with-ipyparallel/
        # http://people.duke.edu/~ccc14/sta-663-2016/19C_IPyParallel.html
        self.d_view = self.rc.direct_view(targets='all')  # rc[:] # Direct view into clients

        # Make sure SNUDDA_DATA is set on the workers, this might be needed if ipcluster
        # is started before SNUDDA_DATA is set
        if os.getenv('SNUDDA_DATA') is not None:
            print(f"Setting SNUDDA_DATA environment variable on workers to {os.getenv('SNUDDA_DATA')}")

            self.d_view.execute("import os")
            self.d_view.execute(f"os.environ['SNUDDA_DATA'] = '{os.getenv('SNUDDA_DATA')}'", block=True)

    ############################################################################

    def cleanup_workers(self):

        self.logfile.write(f"Calling cleanup on workers.")
        # Cleanup, and do garbage collection

        clean_cmd = ("sm = None\nsd = None\nsp = None\nspd = None\nnl = None\nsim = None"
                     "\ninner_mask = None\nmin_max = None\nneuron_hv_list = None"
                     "\nsyn_before = None\nsyn_after = None\ninpt = None"
                     "\nmerge_result_syn = None\nmerge_result_gj = None"
                     "\nimport gc\ngc.collect()")

        if self.d_view is not None:
            self.d_view.execute(clean_cmd, block=True)

    ############################################################################

    def stop_parallel(self):

        print("stop_parallel disabled, to keep pool running.")
        # Disable this function, keep the pool running for now
        return

        # if self.rc is not None:
        #    print("Stopping ipyparallel")
        #    self.rc.shutdown(hub=True)

    ############################################################################

    def setup_log_file(self, log_file_name):
        """
        Open log files for writing.

        Args:
            log_file_name (str) : Path to log file
        """

        data_dir = os.path.dirname(log_file_name)

        self.make_dir_if_needed(data_dir)

        try:
            self.logfile = open(log_file_name, 'w')
            self.logfile.write('Starting log file\n')
        except:
            print("Unable to set up log file " + str(log_file_name))

    ############################################################################

    def close_log_file(self):
        """ Close log file. """

        stop = timeit.default_timer()

        print(f"\nExecution time: {stop - self.start:.1f}s")

        self.logfile.write(f"Execution time: {stop - self.start:.1f}s")
        self.logfile.write("End of log. Closing file.")
        self.logfile.close()

    ############################################################################

    @staticmethod
    def make_dir_if_needed(dir_path):
        """ Creates directory if missing. """

        if not os.path.exists(dir_path):
            print("Creating missing directory " + dir_path)
            try:
                os.makedirs(dir_path)
                print("Created directory " + dir_path)
            except:
                print("Failed to create dir " + dir_path)

    ############################################################################

    @staticmethod
    def cleanup_workers_rc(rc):

        if rc is None:
            return

        d_view = rc.direct_view(targets='all')  # rc[:] # Direct view into clients

        clean_cmd = ("sm = None\nsd = None\nsp = None\nspd = None\nnl = None\nsim = None"
                     "\ninner_mask = None\nmin_max = None\nneuron_hv_list = None"
                     "\nsyn_before = None\nsyn_after = None\ninpt = None"
                     "\nmerge_result_syn = None\nmerge_result_gj = None"
                     "\nimport gc\ngc.collect()")

        d_view.execute(clean_cmd, block=True)


##############################################################################


if __name__ == "__main__":

    # This is fix to handle if user calles python from within neuron
    import sys

    if '-python' in sys.argv:
        print("Network_simulate.py called through nrniv, fixing arguments")
        pythonidx = sys.argv.index('-python')
        if len(sys.argv) > pythonidx:
            sys.argv = sys.argv[pythonidx + 1:]

    from .cli import snudda_cli

    snudda_cli()
