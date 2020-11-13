# TODO : Consolidate logging to one class, right now separate function everywhere

from psutil import virtual_memory


def low_memory(threshold=0.1):

    mem = virtual_memory()
    memory_ratio = mem.available / mem.total

    return memory_ratio < threshold


def memory_status():

    mem = virtual_memory()
    return mem.available, mem.total
