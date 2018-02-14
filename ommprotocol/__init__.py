#!/usr/bin/env python
# -*- coding: utf-8 -*-

#################################################
#           insiliChem OpenMM launcher          #
# --------------------------------------------- #
# By Jaime RGP <jaime@insilichem.com> @ 2016    #
#################################################

import sys as _sys
import logging
import os

logger = logging.getLogger(__name__)
if not os.environ.get('OMMPROTOCOL_SLAVE'):
    class CustomFormatter(logging.Formatter):

        CUSTOM_FORMATS = {
            logging.DEBUG: "DEBUG: %(module)s: %(lineno)d: %(message)s",
            logging.INFO: "%(message)s",
            logging.WARNING: "Warning: %(message)s",
            logging.ERROR: "[!] %(message)s",
            logging.CRITICAL: "CRITICAL: %(message)s",
            100: "%(message)s"
        }

        def format(self, record):
            format_orig = self._fmt
            self._fmt = self.CUSTOM_FORMATS.get(record.levelno, format_orig)
            result = logging.Formatter.format(self, record)
            self._fmt = format_orig
            return result

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    formatter = CustomFormatter()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

from ommprotocol.io import prepare_input, statexml2pdb, export_frame_coordinates
from ommprotocol.md import protocol
from ._version import get_versions
__version__ = get_versions()['version']
__copyright__ = "InsiliChem OMMProtocol v{}".format(__version__)
del get_versions

_banner = """
     ____    __      __  __      __  ______
  __    __  ____  ____  ____  ____  __    __
 __    __  __  __  __  __  __  __  ______
__    __  __      __  __      __  __
 ____    __      __  __      __  __
"""


def run_protocol():
    print(_banner[1:])
    print(__copyright__)
    print('-' * max([len(l) for l in _banner.splitlines() + [__copyright__]]))
    handler, cfg, args = prepare_input()
    if not args.check:
        protocol(handler, cfg)


def state_to_pdb():
    if len(_sys.argv) in (3, 4):
        statexml2pdb(*_sys.argv[1:])
    else:
        _sys.exit("StateXML2PDB usage: state2pdb <topology> <xmlstate> [<output.pdb>]")


def export_frame():
    if len(_sys.argv) in (4, 5):
        export_frame_coordinates(*_sys.argv[1:])
    else:
        _sys.exit("exportframe usage: crdfromframe <topology> <trajectory> <nframe> [<output.crd>]")


if __name__ == '__main__':
    run_protocol()
