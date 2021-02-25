import os

# We allow user to use $DATA to specify the Snudda data folder.
# Default is the Snudda/snudda/data folder, but the user can set the SNUDDA_DATA environment variable
#


def snudda_parse_path(path):

    if path and "$DATA" in path:

        if "SNUDDA_DATA" in os.environ:
            data_path_str = os.environ["SNUDDA_DATA"]
        else:
            data_path_str = os.path.join(os.path.dirname(__file__), os.pardir, "data")

        path = os.path.abspath(path.replace("$DATA", data_path_str))

    return path


def snudda_isfile(path):
    return os.path.isfile(snudda_parse_path(path))


def snudda_isdir(path):
    return os.path.isdir(snudda_parse_path(path))


def snudda_path_exists(path):
    return os.path.exists(snudda_parse_path(path))


def snudda_simplify_path(path):

    data_path = snudda_parse_path("$DATA")
    
    if path and data_path in path:
        path = path.replace(data_path, "$DATA")

    return path
