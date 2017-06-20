#!/usr/bin/env python
# -*- coding: utf-8 -*-

#################################################
#           insiliChem OpenMM launcher          #
# --------------------------------------------- #
# By Jaime RGP <jaime@insilichem.com> @ 2016    #
#################################################

"""
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

from simtk import openmm as mm

if sys.version_info.major == 3:
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
    names = []
    for i in range(mm.Platform.getNumPlatforms()):
        platform = mm.Platform.getPlatform(i)
        names.append(platform.getName())
    return names


def extant_file(path):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(path):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
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


def sanitize_path_for_file(path, config_file):
    basepath = os.path.dirname(config_file)
    path = os.path.expanduser(path)
    return os.path.join(basepath, path)


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
    try:
        return getattr(name, attr)
    except AttributeError:
        print('! Value for `{}` not found. Using default `{}`'.format(attr, default))
        return default
