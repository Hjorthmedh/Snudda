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

        with self.subTest(stage="init"):
            run_cli_command("init tiny_parallel --size 10 --overwrite")

        with self.subTest(stage="place"):
            run_cli_command("place tiny_parallel --parallel")

        with self.subTest(stage="detect"):
            run_cli_command("detect tiny_parallel --volumeID Striatum --parallel")

        with self.subTest(stage="prune"):
            run_cli_command("prune tiny_parallel --parallel")

        from shutil import copyfile
        print(f"listdir: {os.listdir()}")
        print(f"parent listdir: {os.listdir('..')}")
        copyfile("../snudda/data/config/input-v10-scaled.json", "tiny_parallel/input.json")

        with self.subTest(stage="input"):
            run_cli_command("input tiny_parallel --input tiny_parallel/input.json --time 1.0 --parallel")

        with self.subTest(stage="parallel-stop"):
            os.system("ipcluster stop")

        with self.subTest(stage="simulate"):
            print("Running nrnivmodl:")
            os.system("nrnivmodl ../snudda/data/cellspecs-v2/mechanisms")
            print("Time to run simulation...")
            run_cli_command("simulate tiny_parallel --time 0.1")

        with self.subTest(stage="init-serial"):
            run_cli_command("init tiny_serial --size 10 --overwrite")

        with self.subTest(stage="place-serial"):
            run_cli_command("place tiny_serial")

        with self.subTest(stage="detect-serial"):
            run_cli_command("detect tiny_serial --volumeID Striatum")

        with self.subTest(stage="prune-serial"):
            run_cli_command("prune tiny_serial")

        copyfile("../snudda/data/config/input-v10-scaled.json", "tiny_serial/input.json")
        with self.subTest(stage="input"):
            run_cli_command("input tiny_serial --input tiny_serial/input.json --time 1.0")
