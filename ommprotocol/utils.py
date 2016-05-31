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

import os
import random
import string

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


def random_string(length=5):
    return ''.join(random.choice(string.ascii_letters) for _ in range(length))


def sanitize_path_for_file(path, config_file):
    basepath = os.path.dirname(config_file)
    path = os.path.expanduser(path)
    return os.path.join(basepath, path)