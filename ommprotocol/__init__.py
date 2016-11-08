#!/usr/bin/env python
# -*- coding: utf-8 -*-

#################################################
#           insiliChem OpenMM launcher          #
# --------------------------------------------- #
# By Jaime RGP <jaime@insilichem.com> @ 2016    #
#################################################

import sys as _sys
from ommprotocol.io import prepare_input, statexml2pdb
from ommprotocol.md import protocol
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

def run_protocol():
    handler, cfg = prepare_input()
    protocol(handler, cfg)

def state_to_pdb():
    if len(_sys.argv) in (3, 4):
        statexml2pdb(*_sys.argv[1:])
    else:
        _sys.exit("StateXML2PDB usage: state2pdb <topology> <xmlstate> [<output.pdb>]")

if __name__ == '__main__':
    run_protocol()


