#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ommprotocol: A command line application to launch
#              MD protocols with OpenMM
# By Jaime RGP <@jaimergp>

import os
import sys
import pytest

def get_file(path):
    prefix = os.environ.get('CONDA_PREFIX', '/etc')
    return os.path.join(prefix, 'share', 'openmm', 'examples', path)