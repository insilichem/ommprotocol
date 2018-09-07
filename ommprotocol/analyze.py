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
import os
import sys
import matplotlib.pyplot as plt
plt.ioff()
from .utils import extant_file, sort_key_for_numeric_suffixes


def plot_log(paths):
    import pandas as pd
    multi_csv = (pd.read_csv(path, sep="\t", index_col=1) for path in paths)
    df = pd.concat(multi_csv, ignore_index=True)
    df.ix[:, 1:5].plot(subplots=True, layout=(2,2))
    plt.show()


def plot_rmsd(trajectories, topology=None, subset=None, output='rmsd.dat', chunksize=100,
              reimage=False):
    import mdtraj
    import numpy as np
    from tqdm import tqdm
    if topology:
        topology = mdtraj.load_topology(topology)
    if subset:
        subset = topology.select(subset)
    trajectories = sorted(trajectories, key=sort_key_for_numeric_suffixes)
    first_frame = mdtraj.load_frame(trajectories[0], 0, top=topology)
    frame_size = first_frame.xyz[0].nbytes
    if reimage:
        first_frame.image_molecules(inplace=True)
    rmsds = []
    for trajectory in tqdm(trajectories, unit='file'):
        _, ext = os.path.splitext(trajectory)
        total, unit_scale = None, None
        if ext.lower() == '.dcd':
            n_frames = round(os.path.getsize(trajectory) / frame_size, -1 * len(str(chunksize)[1:]))
            total = int(n_frames/chunksize)
            unit_scale = chunksize
        itertraj = mdtraj.iterload(trajectory, top=topology, chunk=chunksize)
        tqdm_kwargs = {'total': total, 'unit': 'frames', 'unit_scale': unit_scale,
                       'postfix': {'traj': trajectory}}
        for chunk in tqdm(itertraj, **tqdm_kwargs):
            if reimage:
                chunk.image_molecules(inplace=True)
            rmsd = mdtraj.rmsd(chunk, first_frame, atom_indices=subset) * 10.0 # nm->A
            rmsds.append(rmsd)

    rmsds = np.concatenate(rmsds)
    with open(output, 'w') as f:
        f.write('\n'.join(map(str, rmsds)))
    print('\nWrote RMSD values to', output)
    print('Plotting results...')
    plt.plot(rmsds)
    fig = plt.gca()
    fig.set_title('{}{}'.format(trajectories[0], ' and {} more'.format(len(trajectories[1:])
                                                 if len(trajectories) > 1 else '')))
    fig.set_xlabel('Frames')
    fig.set_ylabel('RMSD (A)')
    plt.show()


def inspect_topology(topology, subset=None):
    def usage():
        print('(i) Subset queries must be either a valid MDTraj DSL query (visit URL\n'
              '    http://mdtraj.org/latest/atom_selection.html) OR a list of atom indices:')
        print('    - resname UNK or resname GLN')
        print('    - [0, 1, 2, 3, 4]')
    import mdtraj
    from ommprotocol.io import SystemHandler
    import numpy as np
    ommtop = SystemHandler.load(topology, strict=False).topology
    top = mdtraj.Topology.from_openmm(ommtop)
    print('Topology', topology)
    print('***')
    print('Contents:')
    print('-', top.n_chains, 'chains')
    print('-', top.n_residues, 'residues')
    print('-', top.n_atoms, 'atoms')
    print('-', top.n_bonds, 'bonds')
    print('***')
    if subset is not None:
        indices = np.empty((0,))
        try:
            indices = top.select(subset)
        except Exception as e1:
            try:
                eval_subset = eval(subset)
                if isinstance(eval_subset, (list, tuple)) and \
                all(isinstance(i, int) for i in eval_subset):
                    indices = np.array(eval_subset)
            except Exception as e2:
                print('!!! Could not parse subset `{}`'.format(subset), 'due to:')
                print('    - MDTraj select raised {}:'.format(e1.__class__.__name__))
                print('     ', str(e1).splitlines()[0])
                print('    - Python eval raised {}:'.format(e2.__class__.__name__))
                print('     ', e2)
                usage()
                return
        print('\nSubset `{}` will select {} atoms.'.format(subset, indices.size))
        if indices.size:
            df, dfbonds = top.to_dataframe()
            report = df.iloc[indices].drop(['serial'], axis=1).to_string()
            if indices.size < 5000:
                print(report)
            else:
                print('Dumping table to file {}.subset...'.format(topology))
                with open('{}.subset'.format(topology), 'w') as f:
                    f.write(report)
        else:
            usage()


def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers()

    # 'log' subcommand
    p_log = sp.add_parser('log', help='Plot .log reports')
    p_log.add_argument('paths', type=extant_file, nargs='+')
    p_log.set_defaults(func=plot_log)

    # 'rmsd' subcommand
    p_rmsd = sp.add_parser('rmsd', help='Plot RMSD of one or more (sorted) trajectories')
    p_rmsd.add_argument('trajectories', help='Trajectory file(s)', type=extant_file,
                        nargs='+')
    p_rmsd.add_argument('-t', '--topology', help='Topology file', type=extant_file)
    p_rmsd.add_argument('-s', '--subset', help='DSL query to select atoms in Topology',
                        type=str, default=None)
    p_rmsd.add_argument('-o', '--output', help='Plain-text file where to write RMSD results',
                        type=str, default='rmsd.dat')
    p_rmsd.add_argument('-c', '--chunksize', help='Frames to be loaded at once',
                        type=int, default=100)
    p_rmsd.add_argument('-i', '--reimage', help='Wrap coordinates in central PBC image (SLOW!)',
                        default=False, action='store_true')
    p_rmsd.set_defaults(func=plot_rmsd)

    # 'top' subcommand
    p_top = sp.add_parser('top', help='Inspect topologies and subsets')
    p_top.add_argument('topology', type=extant_file)
    p_top.add_argument('-s', '--subset', type=str, default=None,
        help='DSL query to select atoms in Topology. Uses MDTraj parser '
             '(http://mdtraj.org/latest/atom_selection.html).')
    p_top.set_defaults(func=inspect_topology)

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