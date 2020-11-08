#!/usr/bin/env python3

# A wrapper script for the touch detection algorithm
#
# Usage:
#
# snudda init <networkPath> --size XXX
# -- Creates an a json config file
#
# snudda place <networkPath>
# -- Cell placement within volumes specified
#
# snudda detect <networkPath> [--hvsize hyperVoxelSize]
# -- Touch detection of putative synapses
#
# snudda prune <networkPath> [--mergeonly]
# -- Prune the synapses
#
# snudda input <networkPath> [--input yourInputConfig]
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
# the Human Brain Project, funded from the European Union’s Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
# and SGA2).
#


import os
import sys
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

    def __init__(self, networkPath):

        if networkPath[-1] == "/":
            self.network_path = networkPath[:-1]
        else:
            self.network_path = networkPath

        # Add current dir to python path
        sys.path.append(os.getcwd())

        self.start = timeit.default_timer()

    ############################################################################

    def help_info(self, args):
        from snudda.snudda_help import snudda_help_text
        print(snudda_help_text())

    ############################################################################

    def init_config(self, args):
        # self.networkPath = args.path
        print("Creating config file")
        print("Network path: " + str(self.network_path))

        assert args.size is not None, \
            "You need to speicfy --size when initialising config for network2"

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
                "Network path " + str(self.network_path) + " already exists" \
                + " (aborting to prevent accidental overwriting)"

        self.make_dir_if_needed(self.network_path)

        num_population_units = args.NumPopulationUnits
        population_unit_centres = args.PopulationUnitCentres
        population_unit_radius = args.PopulationUnitRadius

        config_file = self.network_path + "/network-config.json"
        SnuddaInit(struct_def=struct_def,
                   config_name=config_file,
                   num_population_units=num_population_units,
                   population_unit_centres=population_unit_centres,
                   population_unit_radius=population_unit_radius)

        if args.size > 1e5:
            print("Make sure there is enough disk space in " + str(self.network_path))
            print("Large networks take up ALOT of space")

    ############################################################################

    def place_neurons(self, args):
        # self.networkPath = args.path
        print("Placing neurons")
        print("Network path: " + str(self.network_path))

        config_file = self.network_path + "/network-config.json"
        position_file = self.network_path + "/network-neuron-positions.hdf5"
        log_file_name = self.network_path + "/log/logFile-place-neurons.txt"

        self.setup_log_file(log_file_name)  # sets self.logFile
        self.setup_parallel()  # sets self.dView and self.lbView

        from .place import SnuddaPlace

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        npn = SnuddaPlace(config_file=config_file,
                          log_file=self.logfile,
                          verbose=True,
                          d_view=self.d_view,
                          h5libver=h5libver)

        npn.write_data_HDF5(position_file)

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

        log_dir = self.network_path + "/log"

        config_file = self.network_path + "/network-config.json"
        position_file = self.network_path + "/network-neuron-positions.hdf5"
        log_filename = self.network_path + "/log/logFile-touch-detection.txt"
        save_file = self.network_path + "/voxels/network-putative-synapses.hdf5"

        voxel_dir = self.network_path + "/voxels"
        self.make_dir_if_needed(voxel_dir)

        self.setup_log_file(log_filename)  # sets self.logFile
        self.setup_parallel()  # sets self.dView and self.lbView

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        from .detect import SnuddaDetect

        if args.cont:
            # Continue previous run
            print("Continuing previous touch detection")

            ncv = SnuddaDetect(config_file=config_file,
                               position_file=position_file,
                               logfile=self.logfile,
                               save_file=save_file,
                               slurm_id=self.SlurmID,
                               volume_id=volume_id,
                               rc=self.rc,
                               hyper_voxel_size=hyper_voxel_size,
                               h5libver=h5libver,
                               restart_detection_flag=False)
        else:
            ncv = SnuddaDetect(config_file=config_file,
                               position_file=position_file,
                               logfile=self.logfile,
                               save_file=save_file,
                               slurm_id=self.SlurmID,
                               volume_id=volume_id,
                               rc=self.rc,
                               h5libver=h5libver,
                               hyper_voxel_size=hyper_voxel_size)

        self.stop_parallel()
        self.close_log_file()

    ############################################################################

    def prune_synapses(self, args):
        # self.networkPath = args.path
        print("Prune synapses")
        print("Network path: " + str(self.network_path))

        from .prune import SnuddaPrune

        log_filename = self.network_path + "/log/logFile-synapse-pruning.txt"

        work_log = self.network_path + "/log/network-detect-worklog.hdf5"

        self.setup_log_file(log_filename)  # sets self.logFile
        self.setup_parallel()  # sets self.dView and self.lbView

        # Optionally set this
        scratch_path = None

        if args.mergeonly:
            pre_merge_only = True
        else:
            pre_merge_only = False

        print("preMergeOnly : " + str(pre_merge_only))

        if args.h5legacy:
            h5libver = "earliest"
        else:
            h5libver = "latest"  # default

        ncvp = SnuddaPrune(work_history_file=work_log,
                           logfile=self.logfile,
                           logfile_name=log_filename,
                           d_view=self.d_view, lb_view=self.lb_view,
                           scratch_path=scratch_path,
                           h5libver=h5libver,
                           pre_merge_only=pre_merge_only)

        self.stop_parallel()
        self.close_log_file()

    ############################################################################

    def setup_input(self, args):

        from .input import SnuddaInput

        print("Setting up inputs, assuming input.json exists")
        log_filename = self.network_path + "/log/logFile-setup-input.log"
        self.setup_log_file(log_filename)  # sets self.logFile
        self.setup_parallel()  # sets self.dView and self.lbView

        if "input" in args:
            input_config = args.input
        else:
            input_config = self.network_path + "/input.json"

        if not os.path.isfile(input_config):
            print("Missing input config file: " + str(input_config))
            return

        if args.networkFile:
            network_file = args.networkFile
        else:
            network_file = self.network_path \
                          + "/network-pruned-synapses.hdf5"

        if args.inputFile:
            spike_file = args.inputFile
        else:
            spike_file = self.network_path + "/input-spikes.hdf5"

        if args.time:
            input_time = args.time

        print("Writing input spikes to " + spike_file)

        ni = SnuddaInput(input_config_file=input_config,
                         hdf5_network_file=network_file,
                         spike_data_filename=spike_file,
                         time=input_time,
                         logfile=self.logfile)

        self.stop_parallel()
        self.close_log_file()

    ############################################################################

    def export_to_SONATA(self, args):

        from snudda.ConvertNetwork import ConvertNetwork

        print("Exporting to SONATA format")
        print("Network path: " + str(self.network_path))

        if args.networkFile:
            network_file = args.networkFile
        else:
            network_file = self.network_path \
                          + "/network-pruned-synapses.hdf5"

        if args.input_file:
            input_file = args.input_file
        else:
            input_file = self.network_path + "/input-spikes.hdf5"

        out_dir = self.network_path + "/SONATA/"

        cn = ConvertNetwork(networkFile=network_file,
                            inputFile=input_file,
                            outDir=out_dir)

    ############################################################################

    def simulate(self, args):

        start = timeit.default_timer()

        from .simulate import SnuddaSimulate

        if args.networkFile:
            network_file = args.networkFile
        else:
            network_file = self.network_path \
                          + "/network-pruned-synapses.hdf5"

        if args.inputFile:
            input_file = args.inputFile
        else:
            input_file = self.network_path + "/input-spikes.hdf5"

        self.make_dir_if_needed(self.network_path + "/simulation")

        print("Using input file " + input_file)

        # nWorkers = args.ncores
        # print("Using " + str(nWorkers) + " workers for neuron")

        # Problems with nested symbolic links when the second one is a relative
        # path going beyond the original base path
        if args.mechDir is None:
            mech_dir = os.path.dirname(network_file) + "/mechanisms"

            # !!! problem with paths, testing to create mechanism dir in current dir
            mech_dir = "mechanisms"

            if not os.path.exists(mech_dir):
                m_dir = os.path.dirname(__file__) + "/data/cellspecs-v2/mechanisms"
                os.symlink(m_dir, mech_dir)
        else:
            mech_dir = args.mechDir

        # !!! These are saved in current directory x86_64
        # --- problem since nrnivmodl seems to want a relative path...

        make_mods_str = "nrnivmodl " + mech_dir
        if not os.path.exists('x86_64'):
            print("Please first run: " + make_mods_str)
            exit(-1)
            # I was having problems when running nrnivmodl in the script, but
            # running it manually in bash works... WHY?!!

        # os.system(makeModsStr)

        save_dir = os.path.dirname(network_file) + "/simulation/"

        if not os.path.exists(save_dir):
            print("Creating directory " + save_dir)
            os.makedirs(save_dir, exist_ok=True)

        # Get the SlurmID, used in default file names
        slurm_id = os.getenv('SLURM_JOBID')

        if slurm_id is None:
            slurm_id = str(666)

        print("args: " + str(args))

        if args.voltOut is not None:
            # Save neuron voltage
            if args.voltOut == "default":
                volt_file = save_dir + 'network-voltage-' + slurm_id + '.csv'
            else:
                volt_file = args.voltOut
        else:
            volt_file = None

        if args.spikesOut is None or args.spikesOut == "default":
            spikes_file = save_dir + 'network-output-spikes-' + slurm_id + '.txt'
        else:
            spikes_file = args.spikesOut

        disable_gj = args.disableGJ
        if disable_gj:
            print("!!! WE HAVE DISABLED GAP JUNCTIONS !!!")

        log_file = os.path.dirname(network_file) \
                  + "/log/network-simulation-log.txt"

        log_dir = os.path.dirname(network_file) + "/log"
        if not os.path.exists(log_dir):
            print("Creating directory " + log_dir)
            os.makedirs(log_dir, exist_ok=True)

        from mpi4py import MPI  # This must be imported before neuron, to run parallel
        from neuron import h  # , gui

        pc = h.ParallelContext()

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

        tSim = args.time * 1000  # Convert from s to ms for Neuron simulator

        sim.check_memory_status()
        print("Running simulation for " + str(tSim) + " ms.")
        sim.run(tSim)  # In milliseconds

        print("Simulation done, saving output")
        if spikes_file is not None:
            sim.write_spikes(spikes_file)

        if volt_file is not None:
            sim.write_voltage(volt_file)

        stop = timeit.default_timer()
        if sim.pc.id() == 0:
            print("Program run time: " + str(stop - start))

        # sim.plot()
        exit(0)

        # cmdStr = "nrnivmodl " + mechDir + " && mpiexec -n " + str(nWorkers) + " -map-by socket:OVERSUBSCRIBE python3 " + os.path.dirname(__file__) + " simulate.py " + networkFile + " " + inputFile + " --time " + str(args.time)

        # if(args.voltOut is not None):
        #  cmdStr += " --voltOut " + args.voltOut

        # os.system(cmdStr)

    ############################################################################

    def analyse(self, args):

        print("Add analysis code here, see Network_analyse.py")

    ############################################################################

    def setup_parallel(self):
        self.SlurmID = os.getenv('SLURM_JOBID')

        if self.SlurmID is None:
            self.SlurmID = self.next_run_id()
        else:
            self.SlurmID = int(self.SlurmID)

        self.logfile.write("Using SlurmID: " + str(self.SlurmID))

        if os.getenv('IPYTHON_PROFILE') is not None:

            self.logfile.write('Creating ipyparallel client\n')

            from ipyparallel import Client
            # self.rc = Client(profile=os.getenv('IPYTHON_PROFILE'),
            #            # sshserver='127.0.0.1',
            #            debug=False)

            u_file = os.getenv('IPYTHONDIR') + "/profile_" + os.getenv('IPYTHON_PROFILE') \
                      + "/security/ipcontroller-client.json"
            self.rc = Client(url_file=u_file, timeout=120, debug=False)

            self.logfile.write('Client IDs: ' + str(self.rc.ids))

            # http://davidmasad.com/blog/simulation-with-ipyparallel/
            # http://people.duke.edu/~ccc14/sta-663-2016/19C_IPyParallel.html
            self.d_view = self.rc.direct_view(targets='all')  # rc[:] # Direct view into clients
            self.lb_view = self.rc.load_balanced_view(targets='all')

            # Define nc globally
            # self.dView.execute("nc = None",block=True)
        else:
            self.logfile.write("No IPYTHON_PROFILE enviroment variable set, running in serial")
            self.d_view = None
            self.lb_view = None
            self.rc = None

    ############################################################################

    def stop_parallel(self):

        # Disable this function, keep the pool running for now
        return

        if self.rc is not None:
            print("Stopping ipyparallel")
            self.rc.shutdown(hub=True)

    ############################################################################

    def setup_log_file(self, logFileName):
        data_dir = os.path.dirname(logFileName)

        self.make_dir_if_needed(data_dir)

        try:
            self.logfile = open(logFileName, 'w')
            self.logfile.write('Starting log file\n')
        except:
            print("Unable to set up log file " + str(logFileName))

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
            import pdb
            pdb.set_trace()
            return 0

        print("Using runID = " + str(next_id))

        return next_id

    ############################################################################

    def make_dir_if_needed(self, dir_path):

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
