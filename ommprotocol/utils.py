#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ommprotocol: A command line application to launch
#              MD protocols with OpenMM
# By Jaime RGP <@jaimergp>

"""
ommprotocol.utils
-----------------

Collection of miscelaneous functions
"""

from __future__ import print_function
from contextlib import contextmanager
import argparse
import os
import random
import string
import sys
try:
    import thread
except ImportError:
    import _thread as thread
import threading


if sys.version_info.major == 3:
    basestring = str
    raw_input = input


def assert_not_exists(path, sep='.'):
    """
    If path exists, modify to add a counter in the filename. Useful
    for preventing accidental overrides. For example, if `file.txt`
    exists, check if `file.1.txt` also exists. Repeat until we find
    a non-existing version, such as `file.12.txt`.

    Parameters
    ----------
    path : str
        Path to be checked

    Returns
    -------
    newpath : str
        A modified version of path with a counter right before the extension.
    """
    name, ext = os.path.splitext(path)
    i = 1
    while os.path.exists(path):
        path = '{}{}{}{}'.format(name, sep, i, ext)
        i += 1
    return path


def assertinstance(obj, types):
    """
    Make sure object `obj` is of type `types`. Else, raise TypeError.
    """
    if isinstance(obj, types):
        return obj
    raise TypeError('{} must be instance of {}'.format(obj, types))


def available_platforms():
    from simtk import openmm as mm
    names = []
    for i in range(mm.Platform.getNumPlatforms()):
        platform = mm.Platform.getPlatform(i)
        names.append(platform.getName())
    return names


def available_platforms_properties():
    from simtk import openmm as mm
    for i in range(1, mm.Platform.getNumPlatforms()):
        platform = mm.Platform.getPlatform(i)
        name = platform.getName()
        print('{}\n{}'.format(name, '-'*len(name)))
        for prop in platform.getPropertyNames():
            value = platform.getPropertyDefaultValue(prop)
            print(prop, '(default={!r})'.format(value) if value else '')
        print()


def extant_file(path):
    """
    Check if file exists with argparse
    """
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError("{} does not exist".format(path))
    return path


@contextmanager
def ignored_exceptions(*exceptions):
    try:
        yield
    except exceptions:
        pass


def random_string(length=5):
    return ''.join(random.choice(string.ascii_letters) for _ in range(length))


def sanitize_args_for_file(args, config_file):
    if isinstance(args, basestring):
        args = (args,)
    path = sanitize_path_for_file(args[0], config_file)
    return (path,) + tuple(args[1:])


def sanitize_path_for_file(path, config_file):
    basepath = os.path.dirname(config_file)
    path = os.path.expanduser(path)
    return os.path.join(basepath, path)


def sort_key_for_numeric_suffixes(path, sep='.', suffix_index=-2):
    """
    Sort files taking into account potentially absent suffixes like
        somefile.dcd
        somefile.1000.dcd
        somefile.2000.dcd

    To be used with sorted(..., key=callable).
    """
    chunks = path.split(sep)
    # Remove suffix from path and convert to int
    if chunks[suffix_index].isdigit():
        return sep.join(chunks[:suffix_index] + chunks[suffix_index+1:]), int(chunks[suffix_index])
    return path, 0


def timed_input(prompt, timeout=300.0):
    print(prompt, end='')
    timer = threading.Timer(timeout, thread.interrupt_main)
    astring = None
    try:
        timer.start()
        astring = raw_input()
    except KeyboardInterrupt:
        pass
    timer.cancel()
    return astring


def warned_getattr(name, attr, default):
    if attr is None:
        return default
    try:
        return getattr(name, attr)
    except AttributeError:
        print('! Value for `{}` not found. Using default `{}`'.format(attr, default))
        return default
