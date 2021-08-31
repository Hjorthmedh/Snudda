# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys


sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('../snudda'))

# --- mock

import mock
 
MOCK_MODULES = ['mpi4py', 'NEURON', 'bluepyopt', 'bluepyopt.ephys', 'ipyparallel', 'matplotlib', 'matplotlib.pyplot', 'h5py', 'numpy', 'scipy', 'scipy.cluster', 'scipy.interpolate', 'sonata', 'pyzmq', 'numexpr', 'argparse', 'pyswarms', 'numba']

for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = mock.Mock()


# -- Project information -----------------------------------------------------

project = 'Snudda'
copyright = '2021, Johannes Hjorth'
author = 'Johannes Hjorth'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.autosummary','sphinx.ext.coverage', 'sphinx.ext.viewcode', 'sphinxarg.ext', 'sphinx.ext.napoleon', 'myst_parser']

autodoc_mock_imports = ['mpi4py', 'NEURON', 'bluepyopt', 'ipyparallel', 'matplotlib', 'h5py', 'numpy', 'scipy', 'sonata', 'pyzmq', 'numexpr', 'argparse', 'pyswarms', 'numba']

import snudda


napoleon_google_docstring = True


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

master_doc = 'index'
project = u"snudda"
version = snudda.__version__
release = snudda.__version__

pygments_style = 'sphinx'

autosummary_generate = True
autodoc_default_flags = ['show-inheritance']
autoclass_content = 'both'
tolerate_sphinx_warnings = True

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']




# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}
