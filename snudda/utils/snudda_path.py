import os

# We allow user to use $DATA to specify the Snudda data folder.
# Default is the Snudda/snudda/data folder, but the user can set the SNUDDA_DATA environment variable
#


def snudda_parse_path(path):

    """ Parses a data path, replacing $DATA with the path to SNUDDA_DATA set by environment variable.

    Args:
        path (str) : Path to modify
        """

    if path and ("$DATA" in path or "$SNUDDA_DATA" in path):

        if "SNUDDA_DATA" in os.environ:
            data_path_str = os.environ["SNUDDA_DATA"]
        else:
            data_path_str = os.path.join(os.path.dirname(__file__), os.pardir, "data")

        # Updated so both $DATA and $SNUDDA_DATA is possible to use for SNUDDA_DATA path
        p = path.replace("$DATA", data_path_str).replace("$SNUDDA_DATA", data_path_str)
        path = os.path.realpath(p)

    return path


def snudda_isfile(path):
    """ Checks if path is a file. """
    return os.path.isfile(snudda_parse_path(path))


def snudda_isdir(path):
    """ Checks if path is a directory. """
    return os.path.isdir(snudda_parse_path(path))


def snudda_path_exists(path):
    """ Checks if path exists. """
    return os.path.exists(snudda_parse_path(path))


def snudda_simplify_path(path):
    """ Simplifies path, replacing any occurance of SNUDDA_DATA in the path with $SNUDDA_DATA.

    Args:
        path (str) : Path to be simplified
    """
    data_path = snudda_parse_path("$SNUDDA_DATA")
    real_path = os.path.realpath(path)
    
    if path and data_path in real_path:
        path = real_path.replace(data_path, "$SNUDDA_DATA")

    return path
