# TODO : Consolidate logging to one class, right now separate function everywhere

from psutil import virtual_memory


def memory_status():

    """ Returns memory status. Tuple (memory available, memory total) """

    mem = virtual_memory()
    return mem.available, mem.total
