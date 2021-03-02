import unittest, os, sys, argparse, time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import snudda.cli


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
        with self.subTest(stage="create"):
            run_cli_command("create test-project --overwrite")
        os.chdir('test-project')

        with self.subTest(stage="setup-parallel"):
            os.environ["IPYTHONDIR"] = os.path.join(os.path.abspath(os.getcwd()), ".ipython")
            os.environ["IPYTHON_PROFILE"] = "Snudda_local"
            os.system("ipcluster start -n 4 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&")
            time.sleep(10)

        with self.subTest(stage="init-parallel"):
            run_cli_command("init tiny_parallel --size 100 --overwrite")

        with self.subTest(stage="place-parallel"):
            run_cli_command("place tiny_parallel --parallel --raytraceBorders")

        with self.subTest(stage="detect-parallel"):
            run_cli_command("detect tiny_parallel --parallel")

        with self.subTest(stage="prune-parallel"):
            run_cli_command("prune tiny_parallel --parallel")

        from shutil import copyfile
        print(f"listdir: {os.listdir()}")
        print(f"parent listdir: {os.listdir('..')}")
        copyfile("../snudda/data/input_config/input-v10-scaled.json", "tiny_parallel/input.json")

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
            os.system("nrnivmodl ../snudda/data/neurons/mechanisms")
            print("Time to run simulation...")
            run_cli_command("simulate tiny_parallel --time 0.1 --voltOut default")

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

        with self.subTest(stage="place-serial"):
            run_cli_command("place tiny_serial --h5legacy")

        with self.subTest(stage="detect-serial"):
            run_cli_command("detect tiny_serial --volumeID Striatum --hvsize 120 --randomseed 123 --verbose --h5legacy")

        with self.subTest(stage="detect-serial-cont"):
            run_cli_command("detect tiny_serial --volumeID Striatum --hvsize 120 --cont --h5legacy")

        with self.subTest(stage="prune-serial-merge-only"):
            run_cli_command("prune tiny_serial --mergeonly --h5legacy")  # Testing the merge only code

        with self.subTest(stage="prune-serial"):
            run_cli_command("prune tiny_serial --h5legacy")

        copyfile("../snudda/data/input_config/input-v10-scaled.json", "tiny_serial/input.json")
        with self.subTest(stage="input"):
            run_cli_command("input tiny_serial --time 1.0 --inputFile tiny_serial/input-spikes.hdf5")
