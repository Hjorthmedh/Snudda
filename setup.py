import setuptools, sys, os
from glob import glob

def get_version():
    version = {}
    with open(os.path.join("snudda", "__init__.py")) as f:
        for line in f:
            if line.startswith("__version__"):
                exec(line, version)
                break
    return version["__version__"]

with open("README.md", "r") as fh:
    long_description = fh.read()


    
# Collect all files recursively from the "snudda/data" folder
data_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), "snudda", "data"))
data_files = []
for (dirpath, dirnames, filenames) in os.walk(data_folder):
    rel_folder = os.path.relpath(dirpath, "snudda")
    if len(filenames) > 0:
        data_files.append(os.path.join(rel_folder, "*"))

if os.environ.get('READTHEDOCS') == 'True':
    # We are in the readthedocs.org environment
    print("READTHEDOCS environment detected, clearing install_requires")
    install_requires = []
else:
    print(f"READTHEDOCS = {os.environ.get('READTHEDOCS')}") 
    install_requires = [
        "bluepyopt>=1.11.7",
        "h5py>=3.2.1",
        "ipyparallel>=6.3.0",
        "matplotlib>=3.3.4",
        "mpi4py>=3.0.3",
        "numpy>=1.26.4", 
        "scipy>=1.6.3",
        "sonata>=0.0.2",   # do we need libsonata?
        "pyzmq>=22.0.3",
        "setuptools",
        "psutil",
        "numexpr>=2.7.3",
        "numba>=0.56.4",
        "wheel",
        "open3d"
        # "igraph"
    ]
    
setuptools.setup(
    name="snudda",
    version=get_version(),
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
    entry_points={"console_scripts": ["snudda = snudda.cli:snudda_cli",
                                      "snudda_load = snudda.utils.load:snudda_load_cli",
                                      "snudda_load_simulation_data = snudda.utils.load_network_simulation:load_network_simulation_cli",
                                      "snudda_ablate_network = snudda.utils.ablate_network:snudda_ablate_network_cli",
                                      "snudda_plot_network = snudda.plotting.plot_network:snudda_plot_network_cli",
                                      "snudda_plot_trace = snudda.plotting.plot_traces:snudda_plot_traces_cli",
                                      "snudda_create_cube_mesh = snudda.place.create_cube_mesh:cube_cli"]},
    install_requires=install_requires,
    extras_require={"dev": ["sphinx"]},
    setup_requires=['wheel']
)
