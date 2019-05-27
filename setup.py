#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import io
import os
import versioneer

VERSION = versioneer.get_version()


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(os.path.join(os.path.dirname(__file__), filename), encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.rst')

setup(
    name='ommprotocol',
    version=VERSION,
    cmdclass=versioneer.get_cmdclass(),
    url='https://github.com/insilichem/ommprotocol',
    download_url='https://github.com/insilichem/ommprotocol/tarball/v' + VERSION,
    license='LGPL',
    author="Jaime Rodríguez-Guerra",
    author_email='jaime.rogue@gmail.com',
    description='Easy to deploy MD protocols for OpenMM',
    long_description=long_description,
    packages=find_packages(),
    package_data={'': ['../examples/*.yaml']},
    platforms='any',
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 3 - Alpha',
        'Natural Language :: English',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    entry_points='''
        [console_scripts]
        ommprotocol=ommprotocol:run_protocol
        ommanalyze=ommprotocol.analyze:main
        state2pdb=ommprotocol:state_to_pdb
        exportframe=ommprotocol:export_frame
    ''',
)
