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

import sys
import pandas as pd
import matplotlib.pyplot as plt
plt.ioff()


def plot_logfile(path):
    df = pd.read_csv(path, sep="\t", index_col=1)
    df.ix[:, 1:5].plot(subplots=True, layout=(2,2))
    plt.show()

def main():
    fname = sys.argv[1]
    plot_logfile(fname)

if __name__ == "__main__":
    main()