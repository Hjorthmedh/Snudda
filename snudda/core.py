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

import pkg_resources
import json

from snudda.utils import snudda_path
from snudda.utils.snudda_path import snudda_isfile


def get_data_file(*dirs):
    path = os.path.join("data", *dirs)
    if not pkg_resources.resource_exists(__package__, path):
        raise FileNotFoundError("Data file '{}' not found".format(path))
    return pkg_resources.resource_filename(__package__, path)


class Snudda(object):

    """ Wrapper class, calls Snudda helper functions """

    def __init__(self, network_path):

        """
        Instantiates Snudda
        :param network_path: Location of Snudda network
        """

        self.network_path = network_path
        self.d_view = None
        self.rc = None
        self.slurm_id = 0

        # Add current dir to python path
        sys.path.append(os.getcwd())

        self.start = timeit.default_timer()

    ############################################################################

    @staticmethod
    def help_info(args):
        """ Prints Snudda help """
        from snudda.help import snudda_help_text

    ############################################################################

    def init_config(self, args):
        """
        Creates network-config.json in network_path.

        Args:
            args : command line arguments from argparse

        Example:
             snudda init -size 100 [-overwrite] [-randomseed 1234] [--profile] [--verbose] path
        """
        # self.networkPath = args.path
        print("Creating config file")
        print(f"Network path: {self.network_path}")

        assert args.size is not None, "You need to specify --size when initialising config for the network"

        from snudda.init.init import SnuddaInit
        struct_def = {"Striatum": args.size,
                      "GPe": 0,
                      "GPi": 0,
                      "SNr": 0,
                      "STN": 0,
                      "Cortex": 0,
                      "Thalamus": 0}
        # Cortex and thalamus axons disabled right now, set to 1 to include one

        if not args.overwrite:
            assert not os.path.exists(self.network_path), \
                (f"Network path {self.network_path} already exists (aborting to prevent accidental overwriting)."
                 "\nCall snudda init with --overwrite to override and overwrite the old data.")

        self.make_dir_if_needed(self.network_path)

        random_seed = args.randomseed

        config_file = os.path.join(self.network_path, "network-config.json")
        SnuddaInit(struct_def=struct_def,
                   neurons_dir=args.neurons_dir,
                   config_file=config_file,
                   random_seed=random_seed,
                   connection_override_file=args.connectionFile)

        if args.size > 1e5:
            print(f"Make sure there is enough disk space in {self.network_path}")
            print("Large networks take up ALOT of space")

    ############################################################################

    def place_neurons(self, args):
        """
        Places neurons in 3D space. Creates network-neuron-positions.hdf5 in network_path.

        Args:
            args : command line arguments from argparse

        Example:
            snudda place [--raytraceBorders] [--profile] [--verbose] [--h5legacy] [-parallel] path

        """
        # self.networkPath = args.path
        print("Placing neurons")
        print(f"Network path: {self.network_path}")

        log_file_name = os.path.join(self.network_path, "log", "place-neurons.txt")

        random_seed = args.randomseed

        self.setup_log_file(log_file_name)  # sets self.logFile

        if args.parallel:
            self.setup_parallel()  # sets self.d_view

        from snudda.place.place import SnuddaPlace

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        sp = SnuddaPlace(network_path=self.network_path,
                         log_file=self.logfile,
                         verbose=args.verbose,
                         d_view=self.d_view,
                         h5libver=h5libver,
                         raytrace_borders=args.raytrace_borders,
                         random_seed=random_seed)

        sp.place()

        self.stop_parallel()
        self.close_log_file()

    ############################################################################

    def touch_detection(self, args):
        """
        Synapse touch detection. Writes results to network_path/voxels (one file per hypervoxel).
        Also adds synapse projections between structures.
        Results written to network-projection-synapses.hdf5 in network_path.

        Args:
            args : command line arguments from argparse

        Example:
            snudda detect [-cont] [-hvsize HVSIZE] [--volumeID VOLUMEID] [--profile] [--verbose] [--h5legacy] [-parallel] path

        """
        # self.networkPath = args.path
        print("Touch detection")
        print("Network path: " + str(self.network_path))

        if args.hvsize is not None:
            hyper_voxel_size = int(args.hvsize)
        else:
            hyper_voxel_size = 100

        if args.volumeID is not None:
            volume_id = args.volumeID
        else:
            volume_id = None

        log_dir = os.path.join(self.network_path, "log")
        if not os.path.exists(log_dir):
            print(f"Creating directory {log_dir}")
            os.makedirs(log_dir, exist_ok=True)

        config_file = os.path.join(self.network_path, "network-config.json")
        position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        log_filename = os.path.join(self.network_path, "log", "touch-detection.txt")
        save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")

        random_seed = args.randomseed

        voxel_dir = os.path.join(self.network_path, "voxels")
        self.make_dir_if_needed(voxel_dir)

        self.setup_log_file(log_filename)  # sets self.logfile

        if args.parallel:
            self.setup_parallel()  # sets self.d_view

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

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
                          verbose=args.verbose)

        if args.cont:
            # Continue previous run
            print("Continuing previous touch detection")
            sd.detect(restart_detection_flag=False)
        else:
            sd.detect(restart_detection_flag=True)

        # Also run SnuddaProject to handle projections between volume

        from snudda.detect.project import SnuddaProject

        sp = SnuddaProject(network_path=self.network_path)
        sp.project()

        self.stop_parallel()
        self.close_log_file()

    ############################################################################

    def prune_synapses(self, args):
        """
        Merges data from synapse detection, then does synapse pruning. Writes results to network-synapses.hdf5 in network_path.

        Args:
            args : command line arguments from argparse

        Example:
            snudda prune [--configFile CONFIG_FILE] [--profile] [--verbose] [--h5legacy] [--keepfiles] [-parallel] path
        """

        # self.networkPath = args.path
        print("Prune synapses")
        print("Network path: " + str(self.network_path))

        from snudda.detect.prune import SnuddaPrune

        log_filename = os.path.join(self.network_path, "log", "synapse-pruning.txt")

        random_seed = args.randomseed

        self.setup_log_file(log_filename)  # sets self.logfile

        if args.parallel:
            self.setup_parallel()  # sets self.d_view

        # Optionally set this
        scratch_path = None

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        sp = SnuddaPrune(network_path=self.network_path,
                         logfile=self.logfile,
                         logfile_name=log_filename,
                         config_file=args.config_file,
                         d_view=self.d_view,
                         scratch_path=scratch_path,
                         h5libver=h5libver,
                         random_seed=random_seed,
                         verbose=args.verbose,
                         keep_files=args.keepfiles or args.savePutative)

        sp.prune()

        if args.savePutative:
            sp.save_putative_synapses()

        self.stop_parallel()
        self.close_log_file()

    ############################################################################

    def setup_input(self, args):
        """
        Creates synaptic input for network based on input.json, writes input-spikes.hdf5 in network_path.
        Args:
            args : command line arguments from argparse

        Example:
            snudda input [--input INPUT] [--inputFile INPUT_FILE] [--networkFile NETWORK_FILE] [--time TIME] [-randomseed 123] [--profile] [--verbose] [--h5legacy] [-parallel] path
        """

        print("Setting up inputs, assuming input.json exists")
        log_filename = os.path.join(self.network_path, "log", "setup-input.txt")
        self.setup_log_file(log_filename)  # sets self.logfile

        if args.parallel:
            self.setup_parallel()  # sets self.d_view

        from snudda.input.input import SnuddaInput

        if "input" in args and args.input:
            input_config = args.input
        else:
            input_config = os.path.join(self.network_path, "input.json")

        if not snudda_isfile(input_config):
            print(f"Missing input config file: {input_config}")
            return

        if args.network_file:
            network_file = args.network_file
        else:
            network_file = os.path.join(self.network_path, "network-synapses.hdf5")

        if args.input_file:
            spike_file = args.input_file
        else:
            spike_file = os.path.join(self.network_path, "input-spikes.hdf5")

        if args.time:
            input_time = args.time
        else:
            input_time = None

        random_seed = args.randomseed

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        print(f"Writing input spikes to {spike_file}")

        si = SnuddaInput(input_config_file=input_config,
                         hdf5_network_file=network_file,
                         spike_data_filename=spike_file,
                         time=input_time,
                         logfile=self.logfile,
                         rc=self.rc,
                         random_seed=random_seed,
                         h5libver=h5libver,
                         verbose=args.verbose)
        si.generate()

        self.stop_parallel()
        self.close_log_file()

    ############################################################################

    def export_to_SONATA(self, args):
        """
        Export network to SONATA files. (Currently not functioning)

        Args:
            args : command line arguments from argparse
        """

        assert False, "Old export to SONATA borken, fixme!"
        # TODO: Fix this
        from snudda.utils.export_sonata import ExportSonata

        print("Exporting to SONATA format")
        print(f"Network path: {self.network_path}")

        if args.network_file:
            network_file = args.network_file
        else:
            network_file = os.path.join(self.network_path, "network-synapses.hdf5")

        if args.input_file:
            input_file = args.input_file
        else:
            input_file = os.path.join(self.network_path, "input-spikes.hdf5")

        out_dir = os.path.join(self.network_path, "SONATA")

        cn = ExportSonata(network_file=network_file, input_file=input_file, out_dir=out_dir)

    ############################################################################

    @staticmethod
    def compile_mechanisms(mech_dir=None):

        if not mech_dir:
            mech_dir = os.path.realpath(snudda_path.snudda_parse_path(os.path.join("$DATA", "neurons", "mechanisms")))

        if not os.path.exists("x86_64") and not os.path.exists("nrnmech.dll"):

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

            if os.path.exists("nrnmech.dll"):
                h.nrn_load_dll("nrnmech.dll")
            elif os.path.exists("x86_64"):
                h.nrn_load_dll("x86_64/.libs/libnrnmech.so")
            else:
                print(f"Could not find compiled mechanisms. Compile using 'nrnivmodl {mech_dir}' "
                      f"and retry simulation.")
                sys.exit(-1)

        else:
            print("NEURON mechanisms already compiled, make sure you have the correct version of NEURON modules."
                  "\nIf you delete x86_64 directory (or nrnmech.dll) "
                  "then you will force a recompilation of the modules.")

    ############################################################################

    def simulate(self, args):
        """
        Simulate network. Writes results to network_path/simulation.

        Args:
            args : command line arguments from argparse

        Example:
            snudda simulate [--networkFile NETWORK_FILE] [--inputFile INPUT_FILE] [--time TIME]
            [--spikesOut SPIKES_OUT] [--neuromodulation NEUROMODULATION] [--noVolt] [--disableGJ]
            [-mechdir MECH_DIR] [--profile] [--verbose] [--exportCoreNeuron] path
        """

        start = timeit.default_timer()

        from mpi4py import MPI  # This must be imported before neuron, to run parallel
        from neuron import h
        pc = h.ParallelContext()

        if args.network_file:
            network_file = args.network_file
        else:
            network_file = os.path.join(self.network_path, "network-synapses.hdf5")

        if args.input_file:
            input_file = args.input_file
        else:
            input_file = os.path.join(self.network_path, "input-spikes.hdf5")

        if args.output_file:
            output_file = args.output_file
        else:
            output_file = os.path.join(self.network_path, "simulation", "network-output.hdf5")

        self.make_dir_if_needed(os.path.join(self.network_path, "simulation"))

        print(f"Using input file {input_file}")

        # nWorkers = args.ncores
        # print("Using " + str(nWorkers) + " workers for neuron")

        # Problems with nested symbolic links when the second one is a relative
        # path going beyond the original base path

        if args.mech_dir:
            mech_dir = args.mech_dir
        else:
            # Take into account which SNUDDA_DATA the user wants to use
            mech_dir = os.path.realpath(snudda_path.snudda_parse_path(os.path.join("$DATA", "neurons", "mechanisms")))

            if args.neuromodulation is not None:
                # read neuromod file and determine if it is replay or adaptive, then if and import the correct one
                with open(args.neuromodulation, "r") as f:
                    neuromod_dict = json.load(f, object_pairs_hook=OrderedDict)

                if "adaptive" in neuromod_dict["type"]:
                    mech_dir = os.path.realpath(snudda_path.snudda_parse_path(os.path.join("$DATA", "neurons",
                                                                                           "mechanisms_ptr")))
        self.compile_mechanisms(mech_dir=mech_dir)

        save_dir = os.path.join(os.path.dirname(network_file), "simulation")

        if not os.path.exists(save_dir):
            print(f"Creating directory {save_dir}")
            os.makedirs(save_dir, exist_ok=True)

        # Get the SlurmID, used in default file names
        slurm_id = os.getenv('SLURM_JOBID')

        if slurm_id is None:
            slurm_id = str(666)

        print(f"args: {args}")

        disable_gj = args.disable_gj
        if disable_gj:
            print("!!! WE HAVE DISABLED GAP JUNCTIONS !!!")

        log_file = os.path.join(os.path.dirname(network_file), "log", "network-simulation-log.txt")

        log_dir = os.path.join(os.path.dirname(network_file), "log")
        if not os.path.exists(log_dir):
            print(f"Creating directory {log_dir}")
            os.makedirs(log_dir, exist_ok=True)

        if args.neuromodulation is not None:

            # read neuromod file and determine if it is replay or adaptive, then if and import the correct one

            with open(args.neuromodulation, 'r') as neuromod_f:
                neuromod_dict = json.load(neuromod_f, object_pairs_hook=OrderedDict)

            if 'type' not in neuromod_dict:
                print(f"Neuromodulation is not defined correctly in {args.neuromodulation} : 'type' is missing. Did you specify the correct file?")
                sys.exit(-1)

            elif 'replay' in neuromod_dict['type']:
                from snudda.neuromodulation.neuromodulation import SnuddaSimulateNeuromodulation

                sim = SnuddaSimulateNeuromodulation(network_file=network_file,
                                                    input_file=input_file,
                                                    disable_gap_junctions=disable_gj,
                                                    log_file=log_file,
                                                    verbose=args.verbose)

                sim.setup()
                sim.add_external_input()
                sim.apply_neuromodulation(neuromod_dict)
                sim.neuromodulation_network_wide()

            elif 'adaptive' in neuromod_dict['type']:
                from snudda.neuromodulation.neuromodulation_synapse import SnuddaSimulateNeuromodulationSynapse

                sim = SnuddaSimulateNeuromodulationSynapse(network_file=network_file,
                                                           input_file=input_file,
                                                           disable_gap_junctions=disable_gj,
                                                           log_file=log_file,
                                                           neuromodulator_description=neuromod_dict)

                sim.setup()
                sim.add_external_input()

        else:

            from snudda.simulate.simulate import SnuddaSimulate

            # Simulate is deterministic, no random seed.
            sim = SnuddaSimulate(network_file=network_file,
                                 input_file=input_file,
                                 disable_gap_junctions=disable_gj,
                                 log_file=log_file,
                                 verbose=args.verbose)
            sim.setup()
            sim.add_external_input()

        sim.check_memory_status()

        if args.record_volt:
            sim.add_recording(side_len=None)  # Side len let you record from a subset
            # sim.addRecordingOfType("dSPN",5) # Side len let you record from a subset

        t_sim = args.time * 1000  # Convert from s to ms for Neuron simulator

        if args.exportCoreNeuron:
            sim.export_to_core_neuron()
            return  # We do not run simulation when exporting to core neuron

        sim.check_memory_status()
        print(f"Running simulation for {t_sim} ms.")
        sim.run(t_sim)  # In milliseconds

        print("Simulation done, saving output")
        sim.write_output(output_file)

        stop = timeit.default_timer()
        if sim.pc.id() == 0:
            print(f"Program run time: {stop - start:.1f}s")

        # sim.plot()

    ############################################################################

    def analyse(self, args):

        print("Add analysis code here, see Network_analyse.py")

    ############################################################################

    def setup_parallel(self):
        """Setup ipyparallel workers."""

        self.slurm_id = os.getenv('SLURM_JOBID')

        if self.slurm_id is None:
            self.slurm_id = 0
        else:
            self.slurm_id = int(self.slurm_id)

        self.logfile.write(f"Using slurm_id: {self.slurm_id}")

        ipython_profile = os.getenv('IPYTHON_PROFILE')
        if not ipython_profile:
            ipython_profile = "default"

        ipython_dir = os.getenv('IPYTHONDIR')
        if not ipython_dir:
            ipython_dir = os.path.join(os.path.abspath(os.getcwd()), ".ipython")

        self.logfile.write('Creating ipyparallel client\n')

        from ipyparallel import Client

        u_file = os.path.join(ipython_dir, f"profile_{ipython_profile}", "security", "ipcontroller-client.json")
        self.rc = Client(url_file=u_file, timeout=120, debug=False)

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

    def stop_parallel(self):

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

        print(f"\nProgram run time: {stop - self.start:.1f}s")

        self.logfile.write(f"Program run time: {stop - self.start:.1f}s")
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
