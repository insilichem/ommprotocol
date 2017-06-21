#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ommprotocol: A command line application to launch
#              MD protocols with OpenMM
# By Jaime RGP <@jaimergp>

import pytest
from ommprotocol.io import SystemHandler, Positions
from ommprotocol.md import Stage
from conftest import get_file

STAGE_OPTIONS = dict(steps=100, 
                     minimization=False, 
                     barostat=False, 
                     verbose=False,
                     trajectory_every=0,
                     restart_every=0,
                     report=False,
                     report_every=0,
                     save_state_at_end=False,
                     attempt_rescue=False)

@pytest.mark.parametrize("top, pos", [
    ('input.pdb', None),
    ('input.prmtop', 'input.inpcrd'),
    ('input.top', 'input.gro')
])
def test_formats(top, pos):
    top = get_file(top)
    kw = dict(positions=Positions.load(get_file(pos))) if pos else {}
    handler = SystemHandler.load(top, **kw)
    stage = Stage(handler, **STAGE_OPTIONS)
    stage.run()
