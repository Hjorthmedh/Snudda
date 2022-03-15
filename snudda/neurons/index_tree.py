import igraph


class IndexTree(igraph.Graph):

    def __init__(self, links, coords):

        super().__init__(directed=True)

        self.add_vertices(coords.shape[0])
        self.add_edges(links)

        self.coords = coords

    # Kom på ett bra namn.
    def walk(self, distance, direction):

        """ Returns all compartments within distance in direction = None, Proximal, Distal """

        pass

    def get_path(self):

        pass

    def verify(self):
        # Check SWC file is not circular
        pass

    def find_parent(self):

        pass

    def find_children(self):

        pass



    # TODO: Funktionalitet vi vill ha
    #
    # -


class IndexLookupError(Exception):
    """Must set a unique value in a BijectiveMap."""

    def __init__(self, value):
        self.value = value
        msg = 'The value "{}" is already in the mapping.'
        super().__init__(msg.format(value))


#  Credit to : https://stackoverflow.com/a/34460187

class IndexLookup(dict):
    """Invertible map."""

    def __init__(self, inverse=None):
        if inverse is None:
            inverse = self.__class__(inverse=self)
        self.inverse = inverse

    def __setitem__(self, key, value):
        if value in self.inverse:
            raise IndexLookupError(value)

        self.inverse._set_item(value, key)
        self._set_item(key, value)

    def __delitem__(self, key):
        self.inverse._del_item(self[key])
        self._del_item(key)

    def _del_item(self, key):
        super().__delitem__(key)

    def _set_item(self, key, value):
        super().__setitem__(key, value)


