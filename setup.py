#!/usr/bin/env python3
import os
from distutils.core import setup
from setuptools import find_packages

classifiers = """
Development Status :: 3 - Alpha
Environment :: Console
License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python
Programming Language :: Python :: 3.6
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

exec(open('mob_suite/version.py').read())

setup(
    name='mob_suite',
    version='1.4.9',

    python_requires='>3.5.2',
    packages=find_packages(exclude=['tests']),
    url='https://github.com/phac-nml/mob-suite',
    license='GPLv3',
    author='James Robertson',
    author_email='james.robertson@canada.ca',
    description=('mob_suite is a set of tools for finding, typing and reconstruction of plasmids from draft and complete genome assemblies.'),
    keywords='Plasmids finding typing reconstruction',
    classifiers=classifiers,
    package_dir={'mob_suite':'mob_suite'},

    install_requires=[
        'numpy>=1.11.1',
        'pandas>=0.18.1',
        'tables>=3.3.0',
        'pandas>=0.22.0',
        'biopython>=1.70',
        'pycurl>=7.43.0',
        'scipy>=1.1.0'
        
    ],

    entry_points={
        'console_scripts': [
            'mob_init=mob_suite.mob_init:main',
            'mob_recon=mob_suite.mob_recon:main',
            'mob_cluster=mob_suite.mob_cluster:main',
            'mob_typer=mob_suite.mob_typer:main',
            'best_blast_hits=mob_suite.blast_best_hits:main',
        ],
    },
)
