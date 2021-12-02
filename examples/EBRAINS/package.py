from spack import *

class Snudda(PythonPackage):
    pypi = 'snudda-1.2.9-py3-none-any.whl'
    maintainers = ['hjorth']
    version('1.2.9', '248ec9320ad80dc4b05dca73e4deff78')

    depends_on('python@3.8:', type=('build','run'))

    depends_on('bluepyopt', type=('build','run'))
    depends_on('h5py', type=('build','run'))
    depends_on('ipyparallel', type=('build','run'))
    depends_on('matplotlib', type=('build','run'))
    depends_on('mpi4py', type=('build','run'))
    depends_on('numpy', type=('build','run'))
    depends_on('scipy', type=('build','run'))
    depends_on('sonata', type=('build','run'))
    depends_on('pyzmq', type=('build','run'))
    depends_on('numexpr', type=('build','run'))
    depends_on('argparse', type=('build','run'))
    depends_on('NEURON', type=('build','run'))
    depends_on('pyswarms', type=('build','run'))
    depends_on('setuptools', type=('build','run'))
    depends_on('psutil', type=('build','run'))
    depends_on('cython', type=('build','run'))
    depends_on('numba', type=('build','run'))
    
