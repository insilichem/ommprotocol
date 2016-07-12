#!/usr/bin/env python
# -*- coding: utf-8 -*-

#################################################
#           insiliChem OpenMM launcher          #
# --------------------------------------------- #
# By Jaime RGP <jaime@insilichem.com> @ 2016    #
#################################################

from ommprotocol.io import prepare_input
from ommprotocol.md import protocol


def run_protocol():
    handler, cfg = prepare_input()
    protocol(handler, cfg)

if __name__ == '__main__':
    run_protocol()
