#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
import os

VERSION = "0.1.1"


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='ommprotocol',
    version=VERSION,
    url='https://github.com/insilichem/ommprotocol',
    download_url='https://github.com/insilichem/ommprotocol/tarball/v' + VERSION,
    license='LGPL',
    author="Jaime Rodríguez-Guerra",
    author_email='jaime.rogue@gmail.com',
    description='Easy to deploy MD protocols for OpenMM',
    long_description=read('README.md'),
    packages=['ommprotocol'],
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
    ''',
)
