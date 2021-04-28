import os
import json
import timeit
from collections import OrderedDict


class BenchmarkLogging:

    def __init__(self, network_path, log_file=None):

        if log_file:
            self.log_file = log_file
        else:
            self.log_file = os.path.join(network_path, "benchmark_log.json")

        self.network_name = self.get_network_name(network_path)

        self.start_time = dict()
        self.end_time = dict()

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
            with open(self.log_file, "r") as fr:
                data = json.load(fr, object_pairs_hook=OrderedDict)
        else:
            data = OrderedDict()

        if self.network_name not in data:
            data[self.network_name] = OrderedDict()

        del_keys = []

        for item_name in self.end_time.keys():
            assert item_name in self.start_time, f"Benchmark logging was not started for {item_name}"
            duration = self.end_time[item_name] - self.start_time[item_name]

            if item_name in data[self.network_name]:
                data[self.network_name][item_name].append(duration)
            else:
                data[self.network_name][item_name] = [duration]

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
