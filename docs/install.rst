.. _install:

============
Installation
============

First method: Standalone installer
----------------------------------

If you haven't used Anaconda or Miniconda before (a Python distribution with a cool package manager), your best bet is to simply download the installer for the latest release, which includes everything you need.

1. Go to the `OMMProtocol releases page <https://github.com/insilichem/ommprotocol/releases/latest>`_ and download the latest installer for your platform.
2. Run the installer and follow the instructions!
    a. In Linux, open the terminal and run ``bash ~/Downloads/ommprotocol*.sh`` or whatever path the file got saved.
    b. In Windows, double click on the downloaded ``ommprotocol*.exe``.
3. The installer will create, by default, a new directory called ``ommprotocol`` in your ``$HOME``. Under ``ommprotocol/bin``(Linux) or ``ommprotocol/Scripts`` (Windows) you will find the ``ommprotocol`` executable.

Second method: Conda package
----------------------------

OMMProtocol is also distributed as a separate ``conda`` package. If you already have Anaconda/Miniconda installed, you won't probably want to download a separate Python distribution. In that case, skip to step 2.

1. Download and install `Miniconda <http://conda.pydata.org/miniconda.html>`_, a tiny Python distribution with a cool package manager and installer. Check `its webpage <http://conda.pydata.org/docs/>`_ for more info.

::

    For Linux:

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3*.sh

    For Windows, download the EXE installer:
    https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe


2. Install it in the default environment...

::

    conda install -c omnia -c insilichem ommprotocol


3. ... or use a new, separate environment (optional):

::

    conda create -n ommprotocol -c omnia -c insilichem ommprotocol
    conda activate ommprotocol


4. If everything is OK, these sould run correctly.

::

        ommprotocol -h
        ommanalyze -h

Third method: From source
-------------------------

If there's no package for your platform, install the dependencies with ``conda`` and then install ``ommprotocol`` from pip or source.

::

    conda create -n ommprotocol -c omnia openmm ruamel_yaml parmed openmoltools mdtraj netcdf4 jinja2 pdbfixer
    conda activate ommprotocol
    # stable version
    pip install ommprotocol
    # dev version
    pip install https://github.com/insilichem/ommprotocol/archive/master.zip


Updating OMMProtocol
--------------------

Depending on the installation method, updating OMMProtocol involves different steps.

First method: Standalone installer
----------------------------------

Just `download the installer for the new version <https://github.com/insilichem/ommprotocol/releases>`_ and run it. In Linux/MacOS you will need to append the ``-U`` flag to the installer. In Windows, just follow the wizard.

Second method: Conda package
----------------------------

Within the activated environment, run ``conda update -c insilichem -c omnia ommprotocol``. That's it.

Third method: From source
-------------------------

Simply pass the ``-U`` flag to pip: ``pip install -U ommprotocol`` or, for development version, ``pip install -U https://github.com/insilichem/ommprotocol/archive/master.zip``.