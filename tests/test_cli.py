import argparse
import os
import sys
import time
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import snudda.cli
from snudda.init.init import SnuddaInit


# Duck punch the argument parser so it doesn't sys.exit
def on_argparse_error(self, message):
    raise argparse.ArgumentError(None, message)


argparse.ArgumentParser.error = on_argparse_error


def run_cli_command(command):
    argv = sys.argv
    sys.argv = command.split(" ")
    sys.argv.insert(0, "test_cli_command")
    result = snudda.cli.snudda_cli()
    sys.argv = argv
    return result


class TestCLI(unittest.TestCase):
    """
        Check if the CLI commands can be executed
    """

    def test_0_basics(self):
        self.assertRaises(argparse.ArgumentError, run_cli_command, "doesntexist")

    def test_workflow(self):

        #with self.subTest(stage="create"):
        #    run_cli_command("create test-project --overwrite")

        if True:
            network_path = "test-project"
            if os.path.exists(network_path):
                import shutil
                shutil.rmtree(network_path)
        
            os.mkdir(network_path)
            os.chdir(network_path)

        with self.subTest(stage="setup-parallel"):
            os.environ["IPYTHONDIR"] = os.path.join(os.path.abspath(os.getcwd()), ".ipython")
            os.environ["IPYTHON_PROFILE"] = "default"
            os.system("ipcluster start -n 4 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&")
            time.sleep(10)

        # with self.subTest(stage="init-parallel-BIG"):
        #     run_cli_command("init tiny_parallel --size 1000000 --overwrite")

        # with self.subTest(stage="place-parallel-BIG"):
        #     run_cli_command("place tiny_parallel --parallel")

        with self.subTest(stage="init-parallel"):
            run_cli_command("init tiny_parallel --size 100 --overwrite")

        # Lets reinit but a smaller network that contains all types of cells, to speed things up
        with self.subTest(stage="small-reinit-1"):
            config_name = os.path.join("tiny_parallel", "network-config.json")
            cnc = SnuddaInit(struct_def={}, config_file=config_name, random_seed=123456)
            cnc.define_striatum(num_dSPN=4, num_iSPN=4, num_FS=2, num_LTS=2, num_ChIN=2,
                                volume_type="cube", stay_inside=True)
            cnc.write_json(config_name)

        with self.subTest(stage="place-parallel"):
            run_cli_command("place tiny_parallel --parallel --stayInside")

        with self.subTest(stage="detect-parallel"):
            run_cli_command("detect tiny_parallel --parallel")

        with self.subTest(stage="prune-parallel"):
            run_cli_command("prune tiny_parallel --parallel")

        from shutil import copyfile
        print(f"listdir: {os.listdir()}")
        print(f"parent listdir: {os.listdir('..')}")
        input_file = os.path.join(os.path.dirname(__file__), os.path.pardir,
                                  "snudda", "data", "input_config", "input-v10-scaled.json")
        copyfile(input_file, os.path.join("tiny_parallel", "input.json"))

        with self.subTest(stage="input"):
            run_cli_command("input tiny_parallel --input tiny_parallel/input.json --parallel")

        # with self.subTest(stage="init-parallel-full"):
        #     run_cli_command("init large_parallel --size 1670000 --overwrite")

        # with self.subTest(stage="place-parallel-full"):
        #    run_cli_command("place large_parallel --parallel")

        with self.subTest(stage="parallel-stop"):
            os.system("ipcluster stop")

        #  Only serial tests below this line, we stopped ipcluster.

        with self.subTest(stage="simulate"):
            print("Running nrnivmodl:")
            # mech_dir = os.path.join(os.path.dirname(__file__), os.path.pardir,
            #                        "snudda", "data", "neurons", "mechanisms")

            mech_dir = os.path.join("..", "validation", "mechanisms")

            if not os.path.exists("mechanisms"):
                print("----> Copying mechanisms")
                # os.symlink(mech_dir, "mechanisms")
                from distutils.dir_util import copy_tree
                copy_tree(mech_dir, "mechanisms")
            else:
                print("------------->   !!! mechanisms already exists")

            eval_str = f"nrnivmodl mechanisms"  # f"nrnivmodl {mech_dir}
            print(f"Running: {eval_str}")
            os.system(eval_str)

            # For the unittest we for some reason need to load mechansism separately
            from mpi4py import MPI  # This must be imported before neuron, to run parallel
            from neuron import h  # , gui
            import neuron

            # For some reason we need to import modules manually
            # when running the unit test.
            if os.path.exists("x86_64/.libs/libnrnmech.so"):
                print("!!! Manually loading libraries")
                try:
                    h.nrn_load_dll("x86_64/.libs/libnrnmech.so")
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)

            if os.path.exists("aarch64/.libs/libnrnmech.so"):
                print("Manually loading libraries")
                try:
                    h.nrn_load_dll("aarch64/.libs/libnrnmech.so")
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)

            print("Time to run simulation...")
            # run_cli_command("simulate tiny_parallel --time 0.1")

            # Use the simulation config interface
            sim_config = os.path.join(os.path.dirname(__file__), "data", "sim-test-config.json")
            run_cli_command(f"simulate tiny_parallel --simulation_config {sim_config}")

        os.environ["SLURM_JOBID"] = "1234"

        with self.subTest(stage="init-serial"):
            # Remove the old folder if it exists
            if os.path.exists("tiny_serial"):
                import shutil
                shutil.rmtree("tiny_serial")

            run_cli_command("init tiny_serial --size 100 --profile")

        with self.subTest(stage="init-overwrite-fail"):
            # Should not allow overwriting of existing folder if --overwrite is not specified
            self.assertRaises(AssertionError, run_cli_command, "init tiny_serial --size 100")

        # Again, let us reinit to a smaller network to speed things up
        with self.subTest(stage="small-reinit-2"):
            config_name = os.path.join("tiny_serial", "network-config.json")
            cnc = SnuddaInit(struct_def={}, config_file=config_name, random_seed=1234)
            cnc.define_striatum(num_dSPN=3, num_iSPN=3, num_FS=2, num_LTS=2, num_ChIN=2,
                                volume_type="cube")
            cnc.write_json(config_name)

        with self.subTest(stage="place-serial"):
            run_cli_command("place tiny_serial --h5legacy")

        with self.subTest(stage="detect-serial"):
            run_cli_command("detect tiny_serial --volumeID Striatum --hvsize 120 --randomseed 123 --verbose --h5legacy")

        with self.subTest(stage="detect-serial-cont"):
            run_cli_command("detect tiny_serial --volumeID Striatum --hvsize 120 --cont --h5legacy")

        with self.subTest(stage="prune-serial"):
            run_cli_command("prune tiny_serial --h5legacy")

        input_file = os.path.join(os.path.dirname(__file__), os.path.pardir,
                                  "snudda", "data", "input_config", "input-v10-scaled.json")
        copyfile(input_file, "tiny_serial/input.json")
        with self.subTest(stage="input"):
            run_cli_command("input tiny_serial --time 1.0 --inputFile tiny_serial/input-spikes.hdf5")

    def tearDown(self) -> None:

        # Exit the test directory
        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))
            