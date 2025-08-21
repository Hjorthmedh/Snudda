import setuptools, os

def get_version():
    version = {}
    with open(os.path.join("snudda", "__init__.py")) as f:
        for line in f:
            if line.startswith("__version__"):
                exec(line, version)
                break
    return version["__version__"]

if __name__ == "__main__":
    setuptools.setup(version=get_version())
