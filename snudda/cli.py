from argparse import ArgumentParser, RawTextHelpFormatter
from snudda.core import Snudda, get_data_file
from snudda.help import snudda_help_text
import os


def snudda_cli():

    parser = ArgumentParser(description=f"Snudda microcircuit generator\n\n{snudda_help_text()}",
                            formatter_class=RawTextHelpFormatter)

    sub_parsers = parser.add_subparsers(help="action", dest="action")
    sub_parsers.required = True

    create_parser = sub_parsers.add_parser("create")
    create_parser.add_argument("path", help="Location of network")
    create_parser.add_argument("-overwrite", "--overwrite", help="Allow overwriting of old directory",
                               action="store_true")
    create_parser.add_argument("--profile", help="Run python cProfile", action="store_true")


    init_parser = sub_parsers.add_parser("init")
    init_parser.add_argument("path", help="Location of network")
    #init_parser.add_argument("size", type=int, help="Number of neurons in network", default=None)
    init_parser.add_argument("-size", "--size", dest="size",
                             type=int, help="Number of neurons in network", default=None)
    init_parser.add_argument("-overwrite", "--overwrite", help="Allow overwriting of old directory",
                             action="store_true")
    init_parser.add_argument("-randomseed", "--randomseed", default=None, help="Random seed", type=int)
    init_parser.add_argument("--NumPopulationUnits", type=int,
                             help="Number of Population Units in the structure, affects connectivity and input correlation",
                             default=1)
    init_parser.add_argument("--PopulationUnitCentres", help="A list which defines the population unit centres",
                             default="[[]]")
    init_parser.add_argument("--PopulationUnitRadius", type=float, help="Radius of population units", default=1)
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
    prune_parser.add_argument("--mergeonly", "--onlymerge", dest="merge_only",
                              help="Merge hyper voxel synapse files, skipping pruning step", action="store_true")
    prune_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    prune_parser.add_argument("--verbose", action="store_true")
    prune_parser.add_argument("--h5legacy", help="Use legacy hdf5 support", action="store_true")
    prune_parser.add_argument("-parallel", "--parallel", action="store_true", default=False)

    input_parser = sub_parsers.add_parser("input")
    input_parser.add_argument("path", help="Location of network")
    input_parser.add_argument("--input", help="Input json config file (for input setup)")
    input_parser.add_argument("--inputFile", help="Input hdf5 file (for simulation)",
                              dest="input_file")
    input_parser.add_argument("--networkFile", help="Network file, if not network-pruned-synapses.hdf5",
                              dest="network_file")
    input_parser.add_argument("--time", type=float, default=2.5, help="Duration of simulation in seconds")
    input_parser.add_argument("-randomseed", "--randomseed", default=None, help="Random seed", type=int)
    input_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    input_parser.add_argument("--verbose", action="store_true")
    input_parser.add_argument("--h5legacy", help="Use legacy hdf5 support", action="store_true")
    input_parser.add_argument("-parallel", "--parallel", action="store_true", default=False)

    simulate_parser = sub_parsers.add_parser("simulate")
    simulate_parser.add_argument("path", help="Location of network")
    simulate_parser.add_argument("--inputFile", help="Input hdf5 file (for simulation)",
                                 dest="input_file")
    simulate_parser.add_argument("--networkFile", help="Network file, if not network-pruned-synapses.hdf5",
                                 dest="network_file")

    simulate_parser.add_argument("--time", type=float, default=2.5, help="Duration of simulation in seconds")

    simulate_parser.add_argument("--voltOut", "--voltout", dest="volt_out", default=None,
                                 help="Name of voltage output file (csv)")
    simulate_parser.add_argument("--spikesOut", "--spikesout", dest="spikes_out", default=None,
                                 help="Name of spike output file (csv)")
    simulate_parser.add_argument("-randomseed", "--randomseed", default=None, help="Random seed", type=int)
    simulate_parser.add_argument("--disableGJ", action="store_true", dest="disable_gj",
                                 help="Disable gap junctions")

    simulate_parser.add_argument("-mechdir", "--mechDir", dest="mech_dir",
                                 help="mechanism directory if not default", default=None)
    simulate_parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    simulate_parser.add_argument("--verbose", action="store_true")

    export_parser = sub_parsers.add_parser("export")
    export_parser.add_argument("path", help="Location of network")
    export_parser.add_argument("--inputFile", help="Input hdf5 file (for simulation)",
                               dest="input_file")
    export_parser.add_argument("--networkFile", help="Network file, if not network-pruned-synapses.hdf5",
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
               "help": snudda.help_info,
               "create": create_project}

    if args.profile:
        prof_file = "profile-" + args.action + ".prof"
        print("Saving profile data to: " + prof_file)
        import cProfile
        cProfile.runctx("actions[args.action](args)", None, locals(), filename=prof_file)

        # To analyse profile data:
        import pstats
        from pstats import SortKey
        p = pstats.Stats(prof_file)
        p.strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats(30)

    else:
        # Performing the requested action
        actions[args.action](args)


def create_project(args):

    if not args.overwrite and os.path.exists(args.path):
        if input("Directory '{}' exists. Are you sure you wish to overwrite it? [y/n] ".format(
                args.path)).lower() != "y":
            print("Project creation aborted")
            return
        else:
            # Delete the existing folder
            import shutil
            shutil.rmtree(args.path)

    from distutils.dir_util import copy_tree

    # Copy the root data files folder to the specified path.
    # The root data folder is the "snudda/data" folder, containing config, synapses & cellspecs.
    copy_tree(get_data_file(), args.path)


if __name__ == "__main__":
    snudda_cli()
