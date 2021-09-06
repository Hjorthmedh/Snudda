import os
import json
import timeit
from collections import OrderedDict


class BenchmarkLogging:

    """ Saves benchmark logging when running snudda from command line.

    Trivial example below:

    import time
    bl = BenchmarkLogging("/home/hjorth/HBP/Snudda/snudda")
    bl.start_timer("test")
    time.sleep(1)
    bl.start_timer("test2")
    time.sleep(2)
    bl.stop_timer("test2")
    bl.stop_timer("test")
    bl.write_log()

    """

    def __init__(self, network_path, parallel_flag=False, log_file=None, running_neuron=False):

        """
        Constructor.

        Args:
            network_path (str): Path to network
            parallel_flag (bool): Running in parallel, should we determine number of workers?
            log_file (str) : Log file to save text to
            running_neuron (bool) : Are we running NEURON? (Sets method for determining number of workers)
        """

        if log_file:
            self.log_file = log_file
        else:
            self.log_file = os.path.join(network_path, "benchmark_log.json")

        self.network_name = self.get_network_name(network_path)

        self.start_time = dict()
        self.end_time = dict()
        self.pc = None  # Used if running neuron

        if parallel_flag or running_neuron:
            self.num_workers = self.get_number_of_workers(running_neuron=running_neuron)
        else:
            self.num_workers = 1

    def get_number_of_workers(self, running_neuron):

        """
        Returns number of workers.

        Args:
            running_neuron (bool) : Are we running NEURON? Used when determining number of workers).

        """

        if running_neuron:
            # We are running neuron, different way to detect number of workers
            from mpi4py import MPI
            from neuron import h
            self.pc = h.ParallelContext()
            return self.pc.nhost()

            # Is there a simpler way to get the number of workers?

        ipython_profile = os.getenv('IPYTHON_PROFILE')
        if not ipython_profile:
            ipython_profile = "default"

        ipython_dir = os.getenv('IPYTHONDIR')
        if not ipython_dir:
            ipython_dir = os.path.join(os.path.abspath(os.getcwd()), ".ipython")

        import ipyparallel
        u_file = os.path.join(ipython_dir, f"profile_{ipython_profile}", "security", "ipcontroller-client.json")
        rc = ipyparallel.Client(url_file=u_file, timeout=120, debug=False)
        d_view = rc.direct_view(targets='all')  # rc[:] # Direct view into clients

        return len(d_view) + 1  # We also include the master node

    @staticmethod
    def get_network_name(network_path):

        """ Returns network name based on network_path. """

        if os.path.basename(network_path):
            network_name = os.path.basename(network_path)
        else:
            network_name = os.path.basename(os.path.dirname(network_path))

        return network_name

    def start_timer(self, item_name):

        """ Start benchmark timer for item_name. """

        self.start_time[item_name] = timeit.default_timer()

    def stop_timer(self, item_name):
        """ Stops benchmark timer for item_name. """
        self.end_time[item_name] = timeit.default_timer()

    def write_log(self):

        """ Writes to log. """

        if self.pc and self.pc.id() != 0:
            # If this is true we are running NEURON, and with id != 0 we are not master node, just return
            return

        if os.path.exists(self.log_file):
            try:
                with open(self.log_file, "r") as fr:
                    data = json.load(fr, object_pairs_hook=OrderedDict)
            except:
                print(f"Error loading {self.log_file}, creating new benchmark log.")
                # Start with fresh data, this will overwrite old data
                data = OrderedDict()

        else:
            data = OrderedDict()

        if self.network_name not in data:
            data[self.network_name] = OrderedDict()

        del_keys = []

        for item_name in self.end_time.keys():
            assert item_name in self.start_time, f"Benchmark logging was not started for {item_name}"
            duration = self.end_time[item_name] - self.start_time[item_name]

            if item_name in data[self.network_name]:
                data[self.network_name][item_name].append([duration, self.num_workers])
            else:
                data[self.network_name][item_name] = [[duration, self.num_workers]]

            del_keys.append(item_name)

        for key in del_keys:
            # Remove the old start and end times
            del self.start_time[key]
            del self.end_time[key]

        with open(self.log_file, "w") as fw:
            json.dump(data, fw, indent=4)

