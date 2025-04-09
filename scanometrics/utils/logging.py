"""
Implements logging functions throughout ScanOMetrics package
"""


from os import EX_USAGE  # "command used incorrectly" exit code
from sys import exit
__VERBOSE_LEVEL__ = 0

def set_verbose(verbose):
    """Set global verbose level.

    :param verbose: verbose level (0: hide everything, 1: warning/errors, 2: notes/prints
    :type verbose: int
    :return: no return as sets global variable __VERBOSE_LEVEL__
    """
    global __VERBOSE_LEVEL__
    if verbose not in [0, 1, 2]:
        raise TypeError('Verbose level must be equal to 0 (hide everything), 1 (warning/errors), or 2 (notes/prints)')
    __VERBOSE_LEVEL__ = verbose

def get_verbose():
    """Return global verbose level"""
    return __VERBOSE_LEVEL__

def PRINT(*args, **kwargs):
    if __VERBOSE_LEVEL__ >= 2:
        print(*args, **kwargs)

def LOG(msg, prefix=''):
    if __VERBOSE_LEVEL__ >= 2:
        print(prefix+"\033[0;32m%s\033[0m" % msg)

def NOTE(msg, prefix=''):
    if __VERBOSE_LEVEL__ == 2:
        print(prefix+"\033[0;30;44m[ NOTE ]\033[0;34m %s\033[0m" % msg)

def WARNING(msg, prefix=''):
    if __VERBOSE_LEVEL__ >= 1:
        print(prefix+"\033[0;30;43m[ WARNING ]\033[0;33m %s\033[0m" % msg)

def ERROR(msg, prefix=''):
    if __VERBOSE_LEVEL__ >= 1:
        msg = prefix+"\033[0;30;41m[ ERROR ]\033[0;31m %s\033[0m\n" % msg
    print(msg)  # Otherwise message isn't printed by exit() in some platforms
    exit(Exception(msg))
