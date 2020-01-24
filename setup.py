import setuptools, sys, os
from snudda import __version__
from glob import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

# Collect all files recursively from the "snudda/data" folder
data_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), "snudda", "data"))
data_files = []
for (dirpath, dirnames, filenames) in os.walk(data_folder):
    rel_folder = os.path.relpath(dirpath, "snudda")
    if len(filenames) > 0:
        data_files.append(os.path.join(rel_folder, "*"))

setuptools.setup(
    name="snudda",
    version=__version__,
    author="Johannes Hjorth",
    author_email="hjorth@kth.se",
    description="Create realistic networks of neurons, synapses placed using touch detection between axons and dendrites.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Hjorthmedh/Snudda",
    license="GPLv3",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    package_data={
        "snudda": data_files,
    },
    entry_points={"console_scripts": ["snudda = snudda.cli:snudda_cli"]},
    install_requires=[
        "bluepyopt>=1.8.21",
        "h5py>=2.8.0",
        "ipyparallel>=6.2.3",
        "matplotlib>=3.0.2",
        "mpi4py>=3.0.1",
        "numpy>=1.15.4",
        "scipy>=1.2.0",
        "sonata>=0.0.1",
        "pyzmq>=18.0.0",
        "setuptools>=41.0.0"
    ],
    extras_require={"dev": ["sphinx"]},
)
