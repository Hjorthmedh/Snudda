import os
import json
import h5py
import numpy as np
import functools


# We allow user to use $DATA to specify the Snudda data folder.
# Default is the Snudda/snudda/data folder, but the user can set the SNUDDA_DATA environment variable
#

# TODO: Add SNUDDA_PATH to network_config file

# snudda_parse_path is slow due to os calls, so cache results
# @functools.cache   # TODO: This is valid from python 3.9.2 -- change to this in the future, for now keep lru_cache
@functools.lru_cache(maxsize=None)
def snudda_parse_path(path, snudda_data):
    """ Parses a data path, replacing $DATA with the path to SNUDDA_DATA set by environment variable.

    Args:
        path (str) : Path to modify
        snudda_data (str) : Path to SNUDDA_DATA, this is optional, and if given overrides environment variable
        """

    if path and ("$DATA" in path or "$SNUDDA_DATA" in path):

        if snudda_data:
            data_path_str = snudda_data
        elif "SNUDDA_DATA" in os.environ:
            data_path_str = os.environ["SNUDDA_DATA"]
        else:
            data_path_str = os.path.join(os.path.dirname(__file__), os.pardir, "data")

        # Updated so both $DATA and $SNUDDA_DATA is possible to use for SNUDDA_DATA path
        p = path.replace("$DATA", data_path_str).replace("$SNUDDA_DATA", data_path_str)
        path = os.path.realpath(p)

    return path


def get_snudda_data(snudda_data=None, config_file=None, network_path=None, verbose=True):

    """ Note this function is slow, and should only be called once and then result stored.

        Args:
            snudda_data
            config_file
            network_path
            verbose
    """

    if config_file is None and network_path is not None:

        config_file2 = os.path.join(network_path, "network-config.json")
        if os.path.isfile(config_file2):
            config_file = config_file2

    if snudda_data is None and config_file is not None:
        snudda_data = read_snudda_data_from_config(config_file)
        if verbose:
            print(f"Reading SNUDDA_DATA={snudda_data} from {config_file}")

    if snudda_data is None and network_path is not None:
        network_file = os.path.join(network_path, "network-synapses.hdf5")
        if os.path.isfile(network_file):
            with h5py.File(network_file, "r") as f:
                if "meta" in f and "snuddaData" in f["meta"]:
                    snudda_data_str = f["meta/snuddaData"][()]
                    if type(snudda_data_str) in [bytes, np.bytes_]:
                        snudda_data = snudda_data_str.decode()
                    else:
                        snudda_data = snudda_data_str

        if snudda_data is not None and verbose:
            print(f"Reading SNUDDA_DATA={snudda_data} from {network_file}")

    if snudda_data is None and "SNUDDA_DATA" in os.environ:
        snudda_data = os.environ["SNUDDA_DATA"]
        if verbose:
            print(f"Reading SNUDDA_DATA={snudda_data} from environment variable $SNUDDA_DATA")

    if snudda_data is None:
        snudda_data = os.path.join(os.path.dirname(__file__), os.pardir, "data")

    assert os.path.isdir(snudda_data), f"SNUDDA_DATA = {snudda_data} DOES NOT EXIST"

    return snudda_data


def read_snudda_data_from_config(config_file):
    snudda_data = None

    with open(config_file, "rt") as f:
        config_data = json.load(f)

    if "SnuddaData" in config_data:
        snudda_data = config_data["SnuddaData"]
        assert os.path.isdir(snudda_data), \
            f"MISSING DIRECTORY SnuddaData = {snudda_data} specified in config file {config_file}"

    return snudda_data


def snudda_isfile(path, snudda_data):
    """ Checks if path is a file. """
    return os.path.isfile(snudda_parse_path(path, snudda_data))


def snudda_isdir(path, snudda_data):
    """ Checks if path is a directory. """
    return os.path.isdir(snudda_parse_path(path, snudda_data))


def snudda_path_exists(path, snudda_data):
    """ Checks if path exists. """
    return os.path.exists(snudda_parse_path(path, snudda_data))


# @functools.cache   # TODO: This is valid from python 3.9.2 -- change to this in the future, for now keep lru_cache
@functools.lru_cache(maxsize=None)
def snudda_simplify_path(path, snudda_data):
    """ Simplifies path, replacing any occurance of SNUDDA_DATA in the path with $SNUDDA_DATA.

    Args:
        path (str) : Path to be simplified
        snudda_data (str) : Path to SNUDDA_DATA, this is optional, and if given overrides environment variable

    """

    if snudda_data:
        data_path = os.path.realpath(snudda_data)
    else:
        data_path = snudda_parse_path("$SNUDDA_DATA", snudda_data=None)

    real_path = os.path.realpath(path)

    if path and data_path in real_path:
        path = real_path.replace(data_path, "$SNUDDA_DATA")

    return path
