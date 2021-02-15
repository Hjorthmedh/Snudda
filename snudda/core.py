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
import sys  # Used in __init__
import timeit
import numpy as np
import zmq
import pkg_resources


def get_data_file(*dirs):
    path = os.path.join("data", *dirs)
    if not pkg_resources.resource_exists(__package__, path):
        raise FileNotFoundError("Data file '{}' not found".format(path))
    return pkg_resources.resource_filename(__package__, path)


class Snudda(object):

    ############################################################################

    def __init__(self, network_path):

        self.network_path = network_path
        self.d_view = None
        self.lb_view = None
        self.rc = None
        self.slurm_id = 0

        # Add current dir to python path
        sys.path.append(os.getcwd())

        self.start = timeit.default_timer()

    ############################################################################

    @staticmethod
    def help_info(args):
        from snudda.help import snudda_help_text

    ############################################################################

    def init_config(self, args):
        # self.networkPath = args.path
        print("Creating config file")
        print(f"Network path: {self.network_path}")

        assert args.size is not None, "You need to specify --size when initialising config for the network"

        from .init import SnuddaInit
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
                "Network path {self.network_path} already exists (aborting to prevent accidental overwriting)"

        self.make_dir_if_needed(self.network_path)

        num_population_units = args.NumPopulationUnits
        population_unit_centres = args.PopulationUnitCentres
        population_unit_radius = args.PopulationUnitRadius
        random_seed = args.randomseed

        config_file = os.path.join(self.network_path, "network-config.json")
        SnuddaInit(struct_def=struct_def,
                   config_file=config_file,
                   num_population_units=num_population_units,
                   population_unit_centres=population_unit_centres,
                   population_unit_radius=population_unit_radius,
                   random_seed=random_seed)

        if args.size > 1e5:
            print(f"Make sure there is enough disk space in {self.network_path}")
            print("Large networks take up ALOT of space")

    ############################################################################

    def place_neurons(self, args):
        # self.networkPath = args.path
        print("Placing neurons")
        print(f"Network path: {self.network_path}")

        config_file = os.path.join(self.network_path, "network-config.json")
        position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        log_file_name = os.path.join(self.network_path, "log", "logFile-place-neurons.txt")

        random_seed = args.randomseed

        self.setup_log_file(log_file_name)  # sets self.logFile

        if args.parallel:
            self.setup_parallel()  # sets self.d_view and self.lb_view

        from .place import SnuddaPlace

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        sp = SnuddaPlace(config_file=config_file,
                         log_file=self.logfile,
                         verbose=True,
                         d_view=self.d_view,
                         h5libver=h5libver,
                         raytrace_borders=args.raytrace_borders,
                         random_seed=random_seed)

        sp.read_config()
        sp.write_data(position_file)

        self.stop_parallel()
        self.close_log_file()

    ############################################################################

    def touch_detection(self, args):
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
        log_filename = os.path.join(self.network_path, "log", "logFile-touch-detection.txt")
        save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")

        random_seed = args.randomseed

        voxel_dir = os.path.join(self.network_path, "voxels")
        self.make_dir_if_needed(voxel_dir)

        self.setup_log_file(log_filename)  # sets self.logfile

        if args.parallel:
            self.setup_parallel()  # sets self.d_view and self.lb_view

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        from .detect import SnuddaDetect

        sd = SnuddaDetect(config_file=config_file,
                          position_file=position_file,
                          logfile=self.logfile,
                          save_file=save_file,
                          slurm_id=self.slurm_id,
                          volume_id=volume_id,
                          rc=self.rc,
                          hyper_voxel_size=hyper_voxel_size,
                          h5libver=h5libver,
                          random_seed=random_seed)

        if args.cont:
            # Continue previous run
            print("Continuing previous touch detection")
            sd.detect(restart_detection_flag=False)
        else:
            sd.detect(restart_detection_flag=True)

        self.stop_parallel()
        self.close_log_file()

    ############################################################################

    def prune_synapses(self, args):
        # self.networkPath = args.path
        print("Prune synapses")
        print("Network path: " + str(self.network_path))

        from .prune import SnuddaPrune

        log_filename = os.path.join(self.network_path, "log", "logFile-synapse-pruning.txt")

        random_seed = args.randomseed

        self.setup_log_file(log_filename)  # sets self.logfile

        if args.parallel:
            self.setup_parallel()  # sets self.d_view and self.lb_view

        # Optionally set this
        scratch_path = None

        if args.merge_only:
            pre_merge_only = True
        else:
            pre_merge_only = False

        print(f"preMergeOnly : {pre_merge_only}")

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        sp = SnuddaPrune(network_path=self.network_path,
                         logfile=self.logfile,
                         logfile_name=log_filename,
                         config_file=args.config_file,
                         d_view=self.d_view, lb_view=self.lb_view,
                         scratch_path=scratch_path,
                         h5libver=h5libver,
                         random_seed=random_seed)

        sp.prune(pre_merge_only=pre_merge_only)

        self.stop_parallel()
        self.close_log_file()

    ############################################################################

    def setup_input(self, args):

        from .input import SnuddaInput

        print("Setting up inputs, assuming input.json exists")
        log_filename = os.path.join(self.network_path, "log", "logFile-setup-input.log")
        self.setup_log_file(log_filename)  # sets self.logfile

        if args.parallel:
            self.setup_parallel()  # sets self.d_view and self.lb_view

        if "input" in args:
            input_config = args.input
        else:
            input_config = os.path.join(self.network_path, "input.json")

        if not os.path.isfile(input_config):
            print(f"Missing input config file: {input_config}")
            return

        if args.network_file:
            network_file = args.network_file
        else:
            network_file = os.path.join(self.network_path, "network-pruned-synapses.hdf5")

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
                         h5libver=h5libver)
        si.generate()

        self.stop_parallel()
        self.close_log_file()

    ############################################################################

    def export_to_SONATA(self, args):

        from snudda.ConvertNetwork import ConvertNetwork

        print("Exporting to SONATA format")
        print(f"Network path: {self.network_path}")

        if args.networkFile:
            network_file = args.network_file
        else:
            network_file = os.path.join(self.network_path, "network-pruned-synapses.hdf5")

        if args.input_file:
            input_file = args.input_file
        else:
            input_file = os.path.join(self.network_path, "input-spikes.hdf5")

        out_dir = os.path.join(self.network_path, "SONATA")

        cn = ConvertNetwork(networkFile=network_file,
                            inputFile=input_file,
                            outDir=out_dir)

    ############################################################################

    def simulate(self, args):

        start = timeit.default_timer()

        from .simulate import SnuddaSimulate

        if args.network_file:
            network_file = args.network_file
        else:
            network_file = os.path.join(self.network_path, "network-pruned-synapses.hdf5")

        if args.input_file:
            input_file = args.input_file
        else:
            input_file = os.path.join(self.network_path, "input-spikes.hdf5")

        self.make_dir_if_needed(os.path.join(self.network_path, "simulation"))

        print(f"Using input file {input_file}")

        # nWorkers = args.ncores
        # print("Using " + str(nWorkers) + " workers for neuron")

        # Problems with nested symbolic links when the second one is a relative
        # path going beyond the original base path
        if args.mech_dir is None:
            # mech_dir = os.path.join(os.path.dirname(network_file), "mechanisms")

            # TODO!!! problem with paths, testing to create mechanism dir in current dir
            mech_dir = "mechanisms"

            if not os.path.exists(mech_dir):
                m_dir = os.path.join(os.path.dirname(__file__), "data", "cellspecs-v2", "mechanisms")
                os.symlink(m_dir, mech_dir)
        else:
            mech_dir = args.mech_dir

        # !!! These are saved in current directory x86_64
        # --- problem since nrnivmodl seems to want a relative path...

        make_mods_str = f"nrnivmodl {mech_dir}"
        if not os.path.exists('x86_64'):
            print(f"Please first run: {make_mods_str}")
            os.sys.exit(-1)
            # I was having problems when running nrnivmodl in the script, but
            # running it manually in bash works... WHY?!!

        # os.system(makeModsStr)

        save_dir = os.path.join(os.path.dirname(network_file), "simulation")

        if not os.path.exists(save_dir):
            print(f"Creating directory {save_dir}")
            os.makedirs(save_dir, exist_ok=True)

        # Get the SlurmID, used in default file names
        slurm_id = os.getenv('SLURM_JOBID')

        if slurm_id is None:
            slurm_id = str(666)

        print(f"args: {args}")

        if args.volt_out is not None:
            # Save neuron voltage
            if args.volt_out == "default":
                volt_file = os.path.join(save_dir, f"network-voltage-{slurm_id}.csv")
            else:
                volt_file = args.volt_out
        else:
            volt_file = None

        if args.spikes_out is None or args.spikes_out == "default":
            spikes_file = os.path.join(save_dir, f"network-output-spikes-{slurm_id}.txt")
        else:
            spikes_file = args.spikes_out

        disable_gj = args.disable_gj
        if disable_gj:
            print("!!! WE HAVE DISABLED GAP JUNCTIONS !!!")

        log_file = os.path.join(os.path.dirname(network_file), "log", "network-simulation-log.txt")

        log_dir = os.path.join(os.path.dirname(network_file), "log")
        if not os.path.exists(log_dir):
            print(f"Creating directory {log_dir}")
            os.makedirs(log_dir, exist_ok=True)

        from mpi4py import MPI  # This must be imported before neuron, to run parallel
        from neuron import h  # , gui

        pc = h.ParallelContext()

        # Simulate is deterministic, no random seed.
        sim = SnuddaSimulate(network_file=network_file,
                             input_file=input_file,
                             disable_gap_junctions=disable_gj,
                             log_file=log_file,
                             verbose=args.verbose)

        sim.add_external_input()

        sim.check_memory_status()

        if volt_file is not None:
            sim.add_recording(side_len=None)  # Side len let you record from a subset
            # sim.addRecordingOfType("dSPN",5) # Side len let you record from a subset

        t_sim = args.time * 1000  # Convert from s to ms for Neuron simulator

        sim.check_memory_status()
        print("Running simulation for " + str(t_sim) + " ms.")
        sim.run(t_sim)  # In milliseconds

        print("Simulation done, saving output")
        if spikes_file is not None:
            sim.write_spikes(spikes_file)

        if volt_file is not None:
            sim.write_voltage(volt_file)

        stop = timeit.default_timer()
        if sim.pc.id() == 0:
            print("Program run time: " + str(stop - start))

        # sim.plot()

    ############################################################################

    def analyse(self, args):

        print("Add analysis code here, see Network_analyse.py")

    ############################################################################

    def setup_parallel(self):

        self.slurm_id = os.getenv('SLURM_JOBID')

        if self.slurm_id is None:
            self.slurm_id = self.next_run_id()
        else:
            self.slurm_id = int(self.slurm_id)

        self.logfile.write(f"Using slurm_id: {self.slurm_id}")

        ipython_profile = os.getenv('IPYTHON_PROFILE')
        if not ipython_profile:
            ipython_profile = "Snudda"

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
        self.lb_view = self.rc.load_balanced_view(targets='all')

    ############################################################################

    def stop_parallel(self):

        # Disable this function, keep the pool running for now
        return

        # if self.rc is not None:
        #    print("Stopping ipyparallel")
        #    self.rc.shutdown(hub=True)

    ############################################################################

    def setup_log_file(self, log_file_name):
        data_dir = os.path.dirname(log_file_name)

        self.make_dir_if_needed(data_dir)

        try:
            self.logfile = open(log_file_name, 'w')
            self.logfile.write('Starting log file\n')
        except:
            print("Unable to set up log file " + str(log_file_name))

    ############################################################################

    def close_log_file(self):

        stop = timeit.default_timer()

        print("\nProgram run time: " + str(stop - self.start))

        self.logfile.write("Program run time: " + str(stop - self.start))
        self.logfile.write("End of log. Closing file.")
        self.logfile.close()

    ##############################################################################

    def next_run_id(self):

        import pickle

        run_id_file = ".runID.pickle"

        try:
            if os.path.isfile(run_id_file):

                with open(run_id_file, 'rb') as f:
                    run_id = pickle.load(f)
                    next_id = int(np.ceil(np.max(run_id)) + 1)

                run_id.append(next_id)

                with open(run_id_file, 'wb') as f:
                    pickle.dump(run_id, f, -1)

            else:

                with open(run_id_file, 'wb') as f:
                    next_id = 1
                    run_id = [1]
                    pickle.dump(run_id, f, -1)

        except Exception as e:
            import traceback
            tstr = traceback.format_exc()
            print(tstr)
            print("Problem reading .runID.pickle file, setting runID to 0")
            return 0

        print("Using runID = " + str(next_id))

        return next_id

    ############################################################################

    @staticmethod
    def make_dir_if_needed(dir_path):

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
