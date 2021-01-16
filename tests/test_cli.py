import unittest, os, sys, argparse

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
        with self.subTest(stage="init"):
            run_cli_command("init tiny --size 10 --overwrite")
        with self.subTest(stage="place"):
            run_cli_command("place tiny")
        with self.subTest(stage="detect"):
            run_cli_command("detect tiny --volumeID Striatum")
        with self.subTest(stage="prune"):
            run_cli_command("prune tiny")
        from shutil import copyfile
        copyfile("data/config/input-tinytest-v9-freq-vectors.json", "tiny/input.json")
        with self.subTest(stage="input"):
            run_cli_command("input tiny --input tiny/input.json --time 1.0")
        with self.subTest(stage="simulate"):
            run_cli_command("simulate tiny --time 0.1")
