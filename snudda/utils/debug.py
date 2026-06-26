import sys

def can_debug():
    return sys.stdin.isatty() and sys.stdout.isatty()