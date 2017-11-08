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


def plot_log(paths):
    import pandas as pd
    multi_csv = (pd.read_csv(path, sep="\t", index_col=1) for path in paths)
    df = pd.concat(multi_csv, ignore_index=True)
    df.ix[:, 1:5].plot(subplots=True, layout=(2,2))
    plt.show()


def plot_rmsd(topology, trajectory):
    import mdtraj
    t = mdtraj.load(trajectory, discard_overlapping_frames=True, top=topology)
    rmsd = mdtraj.rmsd(t, t)*10.0 # nm->Angstroms!
    plt.plot(rmsd)
    fig = plt.gca()
    fig.set_title(', '.join(trajectory))
    fig.set_xlabel('Frames')
    fig.set_ylabel('RMSD (A)')
    plt.show()


def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers()

    # 'log' subcommand
    p_log = sp.add_parser('log', help='Plot .log reports')
    p_log.add_argument('paths', type=extant_file, nargs='+')
    p_log.set_defaults(func=plot_log)

    # 'rmsd' subcommand
    p_rmsd = sp.add_parser('rmsd', help='Plot RMSD of a trajectory')
    p_rmsd.add_argument('topology', help='Topology file', type=extant_file)
    p_rmsd.add_argument('trajectory', help='Trajectory file', type=extant_file,
                        nargs='+')
    p_rmsd.set_defaults(func=plot_rmsd)

    # Autocall each requested method
    cli_args = p.parse_args()
    if not hasattr(cli_args, 'func'):
        p.print_help()
        sys.exit()
    f_code = cli_args.func.__code__
    f_args = [getattr(cli_args, a) for a in f_code.co_varnames[:f_code.co_argcount]]
    cli_args.func(*f_args)


if __name__ == "__main__":
    main()