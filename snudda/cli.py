from argparse import ArgumentParser, RawTextHelpFormatter
from snudda.core import Snudda
from snudda.help import snudda_help_text
from snudda.utils.benchmark_logging import BenchmarkLogging
import os
import sys


def snudda_cli():
    """
    Command line parser for Snudda.
    Valid actions are: init, place, detect, prune, input, simulate, help
    """

    if '-python' in sys.argv:
        print("Snudda's cli.py called through nrniv, fixing arguments")
        pythonidx = sys.argv.index('-python')
        if len(sys.argv) > pythonidx:
            sys.argv = sys.argv[pythonidx + 1:]

    # print(f"Current working directory: {os.getcwd()}")
    
    parser = ArgumentParser(description=f"Snudda microcircuit generator\n\n{snudda_help_text()}",
                            formatter_class=RawTextHelpFormatter)

    sub_parsers = parser.add_subparsers(help="action", dest="action")
    sub_parsers.required = True

    # TODO: Remove the create_parser
    create_parser = sub_parsers.add_parser("create")
    create_parser.add_argument("path", help="Location of network")
    create_parser.add_argument("-overwrite", "--overwrite", help="Allow overwriting of old directory",
                               action="store_true")
    create_parser.add_argument("--profile", help="Run python cProfile", action="store_true")

    init_parser = sub_parsers.add_parser("init")
    init_parser.add_argument("path", help="Location of network")
    # init_parser.add_argument("size", type=int, help="Number of neurons in network", default=None)
    init_parser.add_argument("-size", "--size", dest="size",
                             type=int, help="Number of neurons in network", default=None)
    init_parser.add_argument("--snudda_data", "--SnuddaData", type=str, default=None, dest="snudda_data",
                             help="Path to SNUDDA_DATA")
    init_parser.add_argument("--neurons_dir", type=str, default=None,
                             help="Path to neurons_dir, default is $DATA/neurons (DEPRECATED, use --snudda_data instead")
    init_parser.add_argument("-overwrite", "--overwrite", help="Allow overwriting of old directory",
                             action="store_true")
    init_parser.add_argument("-connectionFile", "--connectionFile", default=None, dest="connection_file",
                             help="Use connectivity from user specified JSON file")
    init_parser.add_argument("--honorStayInside", "--stayInside", default=False, dest="stay_inside", action="store_true")
    init_parser.add_argument("-randomseed", "--randomseed", "--seed", default=None, help="Random seed", type=int)
    init_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    init_parser.add_argument("--verbose", action="store_true")

    import_parser = sub_parsers.add_parser("import")
    import_parser.add_argument("path", help="Location of network")
    import_parser.add_argument("config_file", help="Location of config_file to import")
    import_parser.add_argument("--snudda_data", help="Location of snudda_data", default=None)
    import_parser.add_argument("-overwrite", "--overwrite", action="store_true", help="Overwrite old config file")
    import_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    import_parser.add_argument("--verbose", action="store_true")

    place_parser = sub_parsers.add_parser("place")
    place_parser.add_argument("path", help="Location of network")
    place_parser.add_argument("-randomseed", "--randomseed", "--seed", default=None, help="Random seed", type=int)
    place_parser.add_argument("--honorStayInside", "--stayInside", dest="stay_inside", default=False, action="store_true")
    place_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    place_parser.add_argument("--verbose", action="store_true")
    place_parser.add_argument("--h5legacy", help="Use legacy hdf5 support", action="store_true")
    place_parser.add_argument("-parallel", "--parallel", action="store_true", default=False)
    place_parser.add_argument("-ipython_profile", "--ipython_profile", default=None)
    place_parser.add_argument("-ipython_timeout", "--ipython_timeout", default=120, type=int)


    detect_parser = sub_parsers.add_parser("detect")
    detect_parser.add_argument("path", help="Location of network")
    detect_parser.add_argument("-randomseed", "--randomseed", "--seed", default=None, help="Random seed", type=int)
    detect_parser.add_argument("-cont", "--cont", help="Continue partial touch detection", action="store_true")
    detect_parser.add_argument("-hvsize", "--hvsize", default=100,
                               help="Hyper voxel size, eg. 100 = 100x100x100 voxels in hypervoxel")
    detect_parser.add_argument("--volumeID", help="Specify volume ID for detection step")
    detect_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    detect_parser.add_argument("--verbose", action="store_true")
    detect_parser.add_argument("--h5legacy", help="Use legacy hdf5 support", action="store_true")
    detect_parser.add_argument("-parallel", "--parallel", action="store_true", default=False)
    detect_parser.add_argument("-ipython_profile", "--ipython_profile", default=None)
    detect_parser.add_argument("-ipython_timeout", "--ipython_timeout", default=120, type=int)


    prune_parser = sub_parsers.add_parser("prune")
    prune_parser.add_argument("path", help="Location of network")
    prune_parser.add_argument("-randomseed", "--randomseed", "--seed", default=None, help="Random seed", type=int)
    prune_parser.add_argument("--configFile", dest="config_file", default=None,
                              help="Prune using different network config file, useful when tuning pruning")
    prune_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    prune_parser.add_argument("--verbose", action="store_true")
    prune_parser.add_argument("--h5legacy", help="Use legacy hdf5 support", action="store_true")
    prune_parser.add_argument("--keepfiles", action="store_true",
                              help="Keep temp and voxel files after pruning (e.g. useful if you want to rerun pruning)")
    prune_parser.add_argument("--savePutative", action="store_true",
                              help="Also saved network-putative-synapses.hdf5 with unpruned network")
    prune_parser.add_argument("-parallel", "--parallel", action="store_true", default=False)
    prune_parser.add_argument("-ipython_profile", "--ipython_profile", default=None)
    prune_parser.add_argument("-ipython_timeout", "--ipython_timeout", default=120, type=int)


    input_parser = sub_parsers.add_parser("input")
    input_parser.add_argument("path", help="Location of network")
    input_parser.add_argument("--input", help="Input json config file (for input setup)", default=None)
    input_parser.add_argument("--inputFile", help="Input hdf5 file (for simulation)",
                              dest="input_file", default=None)
    input_parser.add_argument("--networkFile", help="Network file, if not network-synapses.hdf5",
                              dest="network_file")
    input_parser.add_argument("--time", type=float, default=None, help="Duration of simulation in seconds")
    input_parser.add_argument("-randomseed", "--randomseed", "--seed", default=None, help="Random seed", type=int)
    input_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    input_parser.add_argument("--verbose", action="store_true")
    input_parser.add_argument("--h5legacy", help="Use legacy hdf5 support", action="store_true")
    input_parser.add_argument("-parallel", "--parallel", action="store_true", default=False)
    input_parser.add_argument("-ipython_profile", "--ipython_profile", default=None)
    input_parser.add_argument("-ipython_timeout", "--ipython_timeout", default=120, type=int)
    input_parser.add_argument("-no_meta_input", "--no_meta_input", help="Do not use meta.json as stimulation input", action="store_true", default=False)

    simulate_parser = sub_parsers.add_parser("simulate")
    simulate_parser.add_argument("path", help="Location of network")
    simulate_parser.add_argument("--networkFile", "--network_file", help="Network file, if not network-synapses.hdf5",
                                 dest="network_file", default=None)
    simulate_parser.add_argument("--inputFile", "--input_file", help="Input hdf5 file (for simulation)",
                                 dest="input_file", default=None)
    simulate_parser.add_argument("--outputFile", "--output_file", help="Output hdf5 file (from simulation)",
                                 dest="output_file", default=None)
    simulate_parser.add_argument("--time", type=float, default=None, help="Duration of simulation in seconds")

    simulate_parser.add_argument("--snudda_data", "--SnuddaData", type=str, default=None, dest="snudda_data",
                                 help="Path to SNUDDA_DATA")
    simulate_parser.add_argument("--simulation_config", type=str, default=None)

    simulate_parser.add_argument("--noVolt", "--novolt", dest="record_volt", action="store_false",
                                 help="Exclude voltage data, to save time and space.")
    simulate_parser.add_argument("-randomseed", "--randomseed", "--seed", default=None, help="Random seed", type=int)

    simulate_parser.add_argument("--disableSyn", "--disableSynapses", action="store_true", dest="disable_synapses", default=None,
                                 help="Disable synapses")

    simulate_parser.add_argument("--disableGJ", "--disableGapJunctions", action="store_true", dest="disable_gj", default=None,
                                 help="Disable gap junctions")

    simulate_parser.add_argument("--mechdir", "--mechDir", "--mech_dir", dest="mech_dir",
                                 help="mechanism directory if not default", default=None)
    simulate_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    simulate_parser.add_argument("--verbose", action="store_true")
    simulate_parser.add_argument("--exportCoreNeuron", action="store_true")
    simulate_parser.add_argument("--recordALL", dest="record_all", type=str, default=None)
    simulate_parser.add_argument("--enable_rxd_neuromodulation", dest="enable_rxd_neuromodulation", action="store_true", default=None)
    simulate_parser.add_argument("--disable_rxd_neuromodulation", dest="disable_rxd_neuromodulation", action="store_true", default=None)

    export_parser = sub_parsers.add_parser("export")
    export_parser.add_argument("path", help="Location of network")
    export_parser.add_argument("--inputFile", "--input_file", help="Input hdf5 file (for simulation)",
                               dest="input_file")
    export_parser.add_argument("--networkFile", "--network_file", help="Network file, if not network-synapses.hdf5",
                               dest="network_file")
    export_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    export_parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    snudda = Snudda(args.path)

    actions = {"init": snudda.init_config_wrapper,
               "import": snudda.import_config_wrapper,
               "place": snudda.place_neurons_wrapper,
               "detect": snudda.detect_synapses_wrapper,
               "prune": snudda.prune_synapses_wrapper,
               "input": snudda.setup_input_wrapper,
               "export": snudda.export_to_SONATA_wrapper,
               "convert": snudda.export_to_SONATA_wrapper,
               "analyse": snudda.analyse,
               "simulate": snudda.simulate_wrapper,
               "help": snudda.help_info}

    if not hasattr(args, 'ipython_profile'):
        args.ipython_profile = None

    print(f"args.ipython_profile = {args.ipython_profile}")

    if args.profile:
        prof_file = f"profile-{args.action}.prof"
        print(f"Saving profile data to: {prof_file}")
        import cProfile
        cProfile.runctx("actions[args.action](args)", None, locals(), filename=prof_file)

        # To analyse profile data:
        import pstats
        from pstats import SortKey
        p = pstats.Stats(prof_file)
        p.strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats(100)

    else:
        if hasattr(args, "parallel"):
            run_parallel = args.parallel
        else:
            run_parallel = False

        running_neuron = (args.action == "simulate")

        bl = BenchmarkLogging(args.path, parallel_flag=run_parallel, running_neuron=running_neuron,
                              ipython_profile=args.ipython_profile)
        bl.start_timer(args.action)

        # Perform the requested action
        actions[args.action](args)

        bl.stop_timer(args.action)
        bl.write_log()


if __name__ == "__main__":
    snudda_cli()
