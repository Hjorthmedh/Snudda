from collections import OrderedDict
import os
import heapq
import numpy as np
import json
import datetime

from snudda.utils.numpy_encoder import  NumpyEncoder


class ParameterBookkeeper:

    def __init__(self, n_max=10, old_book=None, old_book_file=None):

        self.n_max = n_max

        if old_book:
            self.book = old_book.copy()
            heapq.heapify(self.book)
        else:
            self.book = []

        if old_book_file and os.path.exists(old_book_file):
            self.load(old_book_file)

    def add_parameters(self, parameter_set, section_id, section_x, error, dt=None, volt=None):

        now = datetime.datetime.now().timestamp()

        data = OrderedDict()
        data["parameters"] = parameter_set
        data["section_id"] = section_id
        data["section_x"] = section_x
        data["error"] = error

        if dt:
            data["dt"] = dt
        if volt:
            data["volt"] = volt

        if len(self.book) >= self.n_max:
            worst_error = self.book[0][0]
            if error > np.abs(worst_error):
                # The new parameter is worse than the current worst one, don't add since we are already full
                return

        try:
            heapq.heappush(self.book, (-error, now, data))
        except:
            import traceback
            tstr = traceback.format_exc()
            print(tstr)
            exit(-1)

        if len(self.book) > self.n_max:
            heapq.heappop(self.book)  # Throw away largest error

    def merge(self, *other_books):

        heapq.merge(self.book, *other_books)
        while len(self.book) > self.n_max:
            heapq.heappop(self.book)

    def get_dictionary(self):
        data_dict = OrderedDict()
        data_list = []
        for i in range(0, len(self.book)):
            _, _, data = heapq.heappop(self.book)
            data_list.insert(0, data)

        for idx, data in enumerate(data_list):
            data_dict[idx] = data

        return data_dict

    def get_best_parameterset(self):
        if len(self.book) == 0:
            return None
        else:
            sorted_data = sorted(self.book, key=lambda x: -x[0])
            best_data = sorted_data[0]
            return best_data[2]["parameters"]  # First is error

    def get_best_dataset(self):
        if len(self.book) == 0:
            return None
        else:
            sorted_data = sorted(self.book, key=lambda x: -x[0])
            best_data = sorted_data[0]
            return best_data[2]  # First is error, second datetime (tiebreaker)

    def clear(self):
        self.book = []

    def load(self, file_name, clear=False):

        if clear:
            self.clear()

        with open(file_name, "r") as f:
            data = json.load(f)

        for d in data.values():
            self.add_parameters(parameter_set=np.array(d["parameters"]),
                                section_id=np.array(d["section_id"]).astype(int),
                                section_x=np.array(d["section_x"]),
                                error=d["error"],
                                # dt=d["dt"],
                                # volt=np.array(d["volt"])
                                )

    def save(self, file_name):
        print(f"Writing parameter data to {file_name}")
        data_dict = self.get_dictionary()

        with open(file_name, "w") as f:
            json.dump(data_dict, f, indent=4, cls=NumpyEncoder)

    def check_integrity(self):

        if len(self.book) <= 1:
            return

        param_len = np.array([len(d[2]["parameters"]) for d in self.book])
        assert (param_len == param_len[0]).all(), "Parameters should all be same length"

        section_id = [d[2]["section_id"] for d in self.book]
        section_x = [d[2]["section_x"] for d in self.book]

        assert np.array([s == section_id[0] for s in section_id]).all(), \
            "section_id should match for all parametersets"

        assert np.array([s == section_x[0] for s in section_x]).all(), \
            "section_x should match for all parametersets"





