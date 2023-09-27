#!/usr/bin/env python3
import os
from distutils.core import setup
from setuptools import find_packages
from mob_suite.version import __version__
author = 'James Robertson, Kyrylo Bessonov'

classifiers = """
Development Status :: 3 - Alpha
Environment :: Console
License :: OSI Approved :: Apache Software License
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python
Programming Language :: Python :: 3.7
Programming Language :: Python :: 3.8
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


exec(open('mob_suite/version.py').read())

setup(
    name='mob_suite',
    include_package_data=True,
    version=__version__,
    python_requires='>=3.7.0,<4',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    packages=find_packages(exclude=['tests', 'databases']),
    url='https://github.com/phac-nml/mob-suite',
    license='GPLv3',
    author='James Robertson',
    author_email='james.robertson@canada.ca',
    description=(
        'MOB-suite is a set of tools for finding, typing and reconstruction of plasmids from draft and complete genome assemblies.'),
    keywords='Plasmids finding typing reconstruction',
    classifiers=classifiers,
    package_dir={'mob_suite': 'mob_suite'},
    package_data={'mob_suite': ['config.json']},

    install_requires=[
        'numpy>=1.11.1,<1.23.5',
        'tables>=3.3.0,<4',
        'pandas>=0.22.0,<=1.0.5',
        'biopython>=1.8,<2',
        'pycurl>=7.43.0,<8',
        'scipy>=1.1.0,<2',
        'ete3>=3.1.3,<4',
        'six>=1.10,<2',
    ],

    entry_points={
        'console_scripts': [
            'mob_init=mob_suite.mob_init:main',
            'mob_recon=mob_suite.mob_recon:main',
            'mob_cluster=mob_suite.mob_cluster:main',
            'mob_typer=mob_suite.mob_typer:main',
        ],
    },
)
