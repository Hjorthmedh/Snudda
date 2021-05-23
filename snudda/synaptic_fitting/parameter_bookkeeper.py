from collections import OrderedDict
import heapq


class ParameterBookkeeper:

    def __init__(self, n_max=10, old_book=None):

        self.n_max = n_max

        if old_book:
            self.book = old_book.copy()
            heapq.heapify(self.book)
        else:
            self.book = []

    def add_parameters(self, parameter_set, error, dt, volt):

        data = OrderedDict()
        data["params"] = parameter_set
        data["error"] = error
        data["dt"] = dt
        data["volt"] = volt

        heapq.push(self.book, (-error, data))

        if len(self.book) > self.n_max:
            heapq.pop(self.book)  # Throw away largest error

    def merge(self, *other_books):

        heapq.merge(self.book, *other_books)
        while len(self.book) > self.n_max:
            heapq.pop(self.book)

    def get_dictionary(self):
        param_dict = OrderedDict()
        param_list = []
        for i in range(0, len(self.book)):
            _, param = heapq.pop(self.book)
            param_list.insert(0, param)

        for idx, param in enumerate(param_list):
            param_dict[idx] = param

        return param_dict
