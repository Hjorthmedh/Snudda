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
    init_parser.add_argument("--neurons_dir", type=str, default=None,
                             help="Path to neurons_dir, default is $DATA/neurons")
    init_parser.add_argument("-overwrite", "--overwrite", help="Allow overwriting of old directory",
                             action="store_true")
    init_parser.add_argument("-connectionFile", "--connectionFile", default=None,
                             help="Use connectivity from user specified JSON file")
    init_parser.add_argument("-randomseed", "--randomseed", default=None, help="Random seed", type=int)
    init_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    init_parser.add_argument("--verbose", action="store_true")

    place_parser = sub_parsers.add_parser("place")
    place_parser.add_argument("path", help="Location of network")
    place_parser.add_argument("-randomseed", "--randomseed", default=None, help="Random seed", type=int)
    place_parser.add_argument("--raytraceBorders", help="Ray traces for more precise mesh edge detection",
                              action="store_true", dest="raytrace_borders", default=False)
    place_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    place_parser.add_argument("--verbose", action="store_true")
    place_parser.add_argument("--h5legacy", help="Use legacy hdf5 support", action="store_true")
    place_parser.add_argument("-parallel", "--parallel", action="store_true", default=False)

    detect_parser = sub_parsers.add_parser("detect")
    detect_parser.add_argument("path", help="Location of network")
    detect_parser.add_argument("-randomseed", "--randomseed", default=None, help="Random seed", type=int)
    detect_parser.add_argument("-cont", "--cont", help="Continue partial touch detection", action="store_true")
    detect_parser.add_argument("-hvsize", "--hvsize", default=100,
                               help="Hyper voxel size, eg. 100 = 100x100x100 voxels in hypervoxel")
    detect_parser.add_argument("--volumeID", help="Specify volume ID for detection step")
    detect_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    detect_parser.add_argument("--verbose", action="store_true")
    detect_parser.add_argument("--h5legacy", help="Use legacy hdf5 support", action="store_true")
    detect_parser.add_argument("-parallel", "--parallel", action="store_true", default=False)

    prune_parser = sub_parsers.add_parser("prune")
    prune_parser.add_argument("path", help="Location of network")
    prune_parser.add_argument("-randomseed", "--randomseed", default=None, help="Random seed", type=int)
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

    input_parser = sub_parsers.add_parser("input")
    input_parser.add_argument("path", help="Location of network")
    input_parser.add_argument("--input", help="Input json config file (for input setup)", default=None)
    input_parser.add_argument("--inputFile", help="Input hdf5 file (for simulation)",
                              dest="input_file", default=None)
    input_parser.add_argument("--networkFile", help="Network file, if not network-synapses.hdf5",
                              dest="network_file")
    input_parser.add_argument("--time", type=float, default=None, help="Duration of simulation in seconds")
    input_parser.add_argument("-randomseed", "--randomseed", default=None, help="Random seed", type=int)
    input_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    input_parser.add_argument("--verbose", action="store_true")
    input_parser.add_argument("--h5legacy", help="Use legacy hdf5 support", action="store_true")
    input_parser.add_argument("-parallel", "--parallel", action="store_true", default=False)

    simulate_parser = sub_parsers.add_parser("simulate")
    simulate_parser.add_argument("path", help="Location of network")
    simulate_parser.add_argument("--networkFile", help="Network file, if not network-synapses.hdf5",
                                 dest="network_file", default=None)
    simulate_parser.add_argument("--inputFile", help="Input hdf5 file (for simulation)",
                                 dest="input_file", default=None)
    simulate_parser.add_argument("--outputFile", help="Output hdf5 file (from simulation)",
                                 dest="output_file", default=None)

    simulate_parser.add_argument("--time", type=float, default=2.5, help="Duration of simulation in seconds")

    simulate_parser.add_argument("--noVolt", "--novolt", dest="record_volt", action="store_false",
                                 help="Exclude voltage data, to save time and space.")
    simulate_parser.add_argument("-randomseed", "--randomseed", default=None, help="Random seed", type=int)

    simulate_parser.add_argument("--neuromodulation", type=str, default=None,
                                 help=('replay plays back a vector of modulation level, '
                                       'adaptive sets modulation based on spiking activity'))

    simulate_parser.add_argument("--disableGJ", action="store_true", dest="disable_gj",
                                 help="Disable gap junctions")

    simulate_parser.add_argument("-mechdir", "--mechDir", dest="mech_dir",
                                 help="mechanism directory if not default", default=None)
    simulate_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    simulate_parser.add_argument("--verbose", action="store_true")
    simulate_parser.add_argument("--exportCoreNeuron", action="store_true")

    export_parser = sub_parsers.add_parser("export")
    export_parser.add_argument("path", help="Location of network")
    export_parser.add_argument("--inputFile", help="Input hdf5 file (for simulation)",
                               dest="input_file")
    export_parser.add_argument("--networkFile", help="Network file, if not network-synapses.hdf5",
                               dest="network_file")
    export_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    export_parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    snudda = Snudda(args.path)

    actions = {"init": snudda.init_config,
               "place": snudda.place_neurons,
               "detect": snudda.touch_detection,
               "prune": snudda.prune_synapses,
               "input": snudda.setup_input,
               "export": snudda.export_to_SONATA,
               "convert": snudda.export_to_SONATA,
               "analyse": snudda.analyse,
               "simulate": snudda.simulate,
               "help": snudda.help_info}

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

        bl = BenchmarkLogging(args.path, parallel_flag=run_parallel, running_neuron=running_neuron)
        bl.start_timer(args.action)

        # Perform the requested action
        actions[args.action](args)

        bl.stop_timer(args.action)
        bl.write_log()


if __name__ == "__main__":
    snudda_cli()
