from argparse import ArgumentParser, RawTextHelpFormatter
from snudda.core import Snudda, get_data_file
from snudda.help import snudda_help_text
import os


def snudda_cli():

    parser = ArgumentParser(description="Microcircuit generation\n\n" + snudda_help_text(),
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("action", choices=["init", "create", "place", "detect",
                                           "prune", "input", "export", "analyse", "convert", "simulate", "help"],
                        help="Action to do")
    parser.add_argument("path", help="Storage path for network files")
    parser.add_argument("--size", type=int, help="Number of neurons",
                        default=None)
    parser.add_argument("-cont", "--cont", help="Continue partial touch detection",
                        action="store_true")
    parser.add_argument("-hvsize", "--hvsize",
                        help="Hyper voxel size, 100 good value for full striatum, for small runs, use smaller values to more evenly distribute the workload between workers")
    parser.add_argument("--volumeID", help="Specify volume ID for detection step")
    parser.add_argument("--mergeonly", "--onlymerge",
                        help="Only merge synapses in hyper voxels into a big file. Pre-processing to pruning, normally run before. This allows the user to run this separately.",
                        action="store_true")
    parser.add_argument("--h5legacy", help="Use legacy hdf5 support", action="store_true")
    parser.add_argument("--profile", help="Run python cProfile", action="store_true")
    parser.add_argument("--NumPopulationUnits", type=int,
                        help="Number of Population Units in the structure, affects connectivity and input correlation",
                        default=1)
    parser.add_argument("--PopulationUnitCentres", help="A list which defines the population unit centres",
                        default="[[]]")
    parser.add_argument("--PopulationUnitRadius", type=float, help="Radius of population units", default=1)
    parser.add_argument("--input", help="Input json config file (for input setup)")
    parser.add_argument("--inputFile", help="Input hdf5 file (for simulation)")
    parser.add_argument("--networkFile", help="Network file, if not network-pruned-synapses.hdf5")
    parser.add_argument("--time", type=float, default=2.5,
                        help="Duration of simulation in seconds")
    parser.add_argument("--voltOut", "--voltout",
                        default=None,
                        help="Name of voltage output file (csv)")
    parser.add_argument("--spikesOut", "--spikesout",
                        default=None,
                        help="Name of spike output file (csv)")
    parser.add_argument("--disableGJ", action="store_true",
                        help="Disable gap junctions")

    # parser.add_argument("--ncores", default=12,
    #                    help="Number of cores used for simulation")
    parser.add_argument("--overwrite",
                        help="Skips check if network directory already exists",
                        action="store_true")
    parser.add_argument("--mechDir", help="mechanism directory if not default",
                        default=None)
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    if args.path is not None:
        if args.path[-1] == "/":
            args.path = args.path[:-1]

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
        cProfile.run("actions[args.action](args)", prof_file)

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
