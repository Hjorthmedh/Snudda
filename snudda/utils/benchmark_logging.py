import os
import json
import timeit
from collections import OrderedDict


class BenchmarkLogging:

    def __init__(self, network_path, parallel_flag=False, log_file=None):

        if log_file:
            self.log_file = log_file
        else:
            self.log_file = os.path.join(network_path, "benchmark_log.json")

        self.network_name = self.get_network_name(network_path)

        self.start_time = dict()
        self.end_time = dict()

        if parallel_flag:
            self.num_workers = self.get_number_of_workers()
        else:
            self.num_workers = 1

    @staticmethod
    def get_number_of_workers():
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
        if os.path.basename(network_path):
            network_name = os.path.basename(network_path)
        else:
            network_name = os.path.basename(os.path.dirname(network_path))

        return network_name

    def start_timer(self, item_name):
        self.start_time[item_name] = timeit.default_timer()

    def stop_timer(self, item_name):
        self.end_time[item_name] = timeit.default_timer()

    def write_log(self):

        if os.path.exists(self.log_file):
            try:
                with open(self.log_file, "r") as fr:
                    data = json.load(fr, object_pairs_hook=OrderedDict)
            except:
                try:
                    import datetime
                    cur_time = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
                    backup_file = f"{self.log_file}-old-{cur_time}"
                    os.rename(self.log_file, backup_file)
                    self.write_log(f"Renamed corrupt {self.log_file} as {backup_file}")
                except:
                    self.write_log(f"Failed to create {backup_file}")

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


if __name__ == "__main__":
    import time
    bl = BenchmarkLogging("/home/hjorth/HBP/Snudda/snudda")
    bl.start_timer("test")
    time.sleep(1)
    bl.start_timer("test2")
    time.sleep(2)
    bl.stop_timer("test2")
    bl.stop_timer("test")
    bl.write_log()
