import setuptools, sys, os
from snudda import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="snudda",
    version=__version__,
    author="Johannes Hjorth",
    author_email="robingilbert.deschepper@unipv.it",
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
        "pyzmq>=18.0.0"
    ],
    extras_require={"dev": ["sphinx"]},
)
