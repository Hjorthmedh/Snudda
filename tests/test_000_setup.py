import os
import unittest
import shutil
from pathlib import Path


def cleanup_neuron_folders(root_dir="."):
    """Delete NEURON special folders if they exist"""
    root_path = Path(root_dir)

    # Folders that NEURON creates
    neuron_folders = ['x86_64', 'aarch64', 'arm64', 'i686', 'powerpc', 'umac']

    deleted = []
    for folder_name in neuron_folders:
        folder = root_path / folder_name
        if folder.exists() and folder.is_dir():
            try:
                shutil.rmtree(folder)
                deleted.append(str(folder))
                print(f"Deleted: {folder}")
            except Exception as e:
                print(f"Error deleting {folder}: {e}")

    return deleted

def compile_mod_files(path):
    mech_dir = os.path.abspath(path)

    os.system(f"nrnivmodl {mech_dir}")

class CleanupTest(unittest.TestCase):

    def setUp(self):
        cleanup_neuron_folders(".")
        mech_path = os.path.join(os.path.dirname(__file__), "validation/mechanisms")
        compile_mod_files(mech_path)

    def test_cleanup_complete(self):
        """Dummy test - just confirms setup ran"""
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()
