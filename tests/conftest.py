#!/usr/bin/env python
# -*- coding: utf-8 -*-

#################################################
#           insiliChem OpenMM launcher          #
# --------------------------------------------- #
# By Jaime RGP <jaime@insilichem.com> @ 2016    #
#################################################

import os
import sys
import pytest

def get_file(path):
    prefix = os.environ.get('CONDA_PREFIX', '/etc')
    return os.path.join(prefix, 'share', 'openmm', 'examples', path)