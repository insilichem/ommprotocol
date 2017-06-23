#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ommprotocol: A command line application to launch
#              MD protocols with OpenMM
# By Jaime RGP <@jaimergp>

"""
ommprotocol.analyze
-------------------

Trajectory analysis routines
"""

import argparse
import sys
import matplotlib.pyplot as plt
plt.ioff()
from .utils import extant_file


def plot_log(path):
    import pandas as pd
    df = pd.read_csv(path, sep="\t", index_col=1)
    df.ix[:, 1:5].plot(subplots=True, layout=(2,2))
    plt.show()


def plot_rmsd(trajectory, topology):
    import mdtraj
    t = mdtraj.load(trajectory, top=topology)
    rmsd = mdtraj.rmsd(t, t)*10.0 # nm->Angstroms!
    plt.plot(rmsd)
    fig = plt.gca()
    fig.set_title(trajectory)
    fig.set_xlabel('Frames')
    fig.set_ylabel('RMSD (A)')
    plt.show()


def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers()
    
    # 'log' subcommand
    p_log = sp.add_parser('log', help='Plot .log reports')
    p_log.add_argument('path', type=extant_file)
    p_log.set_defaults(func=plot_log)

    # 'rmsd' subcommand
    p_rmsd = sp.add_parser('rmsd', help='Plot RMSD of a trajectory')
    p_rmsd.add_argument('trajectory', help='Trajectory file', type=extant_file)
    p_rmsd.add_argument('topology', help='Topology file', type=extant_file)
    p_rmsd.set_defaults(func=plot_rmsd)

    # Autocall each requested method
    args = p.parse_args()
    if not hasattr(args, 'func'):
        p.print_help()
        sys.exit()
    fcode = args.func.__code__
    args.func(*[getattr(args, arg) for arg in fcode.co_varnames[:fcode.co_argcount]])


if __name__ == "__main__":
    main()