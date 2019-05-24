#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ommprotocol: A command line application to launch
#              MD protocols with OpenMM
# By Jaime RGP <@jaimergp>

import os
from distutils.spawn import find_executable

import pytest
import numpy as np

from ommprotocol.io import SystemHandler, Positions, Restart, prepare_handler
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

HERE = os.path.dirname(os.path.abspath(__file__))


inputs = [('input.pdb', None), ('input.prmtop', 'input.inpcrd')]
if find_executable('mdrun'):
    inputs.append(('input.top', 'input.gro'))
@pytest.mark.parametrize("top, pos", inputs)
def test_formats(top, pos):
    top = get_file(top)
    kw = dict(positions=Positions.load(get_file(pos))) if pos else {}
    handler = SystemHandler.load(top, **kw)
    stage = Stage(handler, **STAGE_OPTIONS)
    stage.run()


def test_prepare_handler():
    cfg = {'topology': get_file('input.prmtop'),
           'positions': get_file('input.inpcrd')}
    handler = prepare_handler(cfg)
    assert handler.topology.getNumAtoms() == len(handler.positions) == 2269
    stage = Stage(handler, **STAGE_OPTIONS)
    stage.run()


def test_prepare_handler_tuples():
    cfg = {'topology': get_file('input.prmtop'),
           'positions': (get_file('input.inpcrd'),) }
    handler = prepare_handler(cfg)
    assert handler.topology.getNumAtoms() == len(handler.positions) == 2269
    stage = Stage(handler, **STAGE_OPTIONS)
    stage.run()


def test_random_velocities_with_checkpoint():
    """
    Fix issue reported in https://github.com/insilichem/ommprotocol/issues/13
    """
    cfg = {
        'topology': os.path.join(HERE, 'data', 'issue-13', 'input.prmtop'),
        # you can also use input_1_md.state below
        'checkpoint': os.path.join(HERE, 'data', 'issue-13', 'input_1_md.rs'),
        'velocities': None
    }
    checkpoint = Restart.load(cfg['checkpoint'])
    handler = prepare_handler(cfg)
    assert handler.velocities != checkpoint.velocities
    all_velocities = []
    for _ in range(3):
        stage = Stage(handler, **STAGE_OPTIONS)
        state = stage.simulation.context.getState(getVelocities=True)
        all_velocities.append(np.around(state.getVelocities()._value, decimals=3))
    assert not np.isclose(all_velocities[0], all_velocities[1]).all()
    assert not np.isclose(all_velocities[0], all_velocities[2]).all()
    assert not np.isclose(all_velocities[1], all_velocities[2]).all()

