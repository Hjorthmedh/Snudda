from argparse import ArgumentParser, RawTextHelpFormatter
from .core import Snudda
from .snudda_help import snudda_help_text

def snudda_cli():
      parser = ArgumentParser(description="Microcircuit generation\n\n" + snudda_help_text(), formatter_class=RawTextHelpFormatter)
      parser.add_argument("action", choices=["init","place","detect",
                                             "prune","input","export","analyse","convert","simulate","help"],
                          help="Action to do")
      parser.add_argument("path", help="Storage path for network files")
      parser.add_argument("--size",type=int,help="Number of neurons",
                          default=None)
      parser.add_argument("-cont", "--cont",help="Continue partial touch detection",
                          action="store_true")
      parser.add_argument("-hvsize", "--hvsize",help="Hyper voxel size, 100 good value for full striatum, for small runs, use smaller values to more evenly distribute the workload between workers")
      parser.add_argument("--volumeID", help="Specify volumeID for detection step")
      parser.add_argument("--mergeonly","--onlymerge", help="Only merge synapses in hyper voxels into a big file. Pre-processing to pruning, normally run before. This allows the user to run this separately.",action="store_true")
      parser.add_argument("--h5legacy",help="Use legacy hdf5 support",action="store_true")
      parser.add_argument("--profile",help="Run python cProfile",action="store_true")
      parser.add_argument("--nchannels",type=int,help="Number of functional channels in the structure, affects connectivity and input correlation",default=1)
      parser.add_argument("--input",help="Input json config file (for input setup)")
      parser.add_argument("--inputFile",help="Input hdf5 file (for simulation)")
      parser.add_argument("--networkFile", help="Network file, if not network-pruned-synapses.hdf5")
      parser.add_argument("--time",type=float,default=2.5,
                          help="Duration of simulation in seconds")
      parser.add_argument("--voltOut","--voltout",
                          default=None,
                          help="Name of voltage output file (csv)")
      parser.add_argument("--ncores", default=12,
                          help="Number of cores used for simulation")

      args = parser.parse_args()

      if(args.path is not None):
        if(args.path[-1] == "/"):
          args.path = args.path[:-1]

      snudda = Snudda(args.path)

      actions = { "init" : snudda.initConfig,
                  "place" : snudda.placeNeurons,
                  "detect" : snudda.touchDetection,
                  "prune" : snudda.pruneSynapses,
                  "input" : snudda.setupInput,
                  "export" : snudda.exportToSONATA,
                  "convert" : snudda.exportToSONATA,
                  "analyse" : snudda.analyse,
                  "simulate" : snudda.simulate,
                  "help" : snudda.helpInfo}



      if(args.profile):
        profFile = "profile-"+args.action+".prof"
        print("Saving profile data to: " + profFile)
        import cProfile
        cProfile.run("actions[args.action](args)", profFile)

        # To analyse profile data:
        import pstats
        from pstats import SortKey
        p = pstats.Stats(profFile)
        p.strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats(30)

      else:
        # Performing the requested action
        actions[args.action](args)
