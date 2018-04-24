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
Programming Language :: Python :: 2.7
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

exec(open('mob_suite/version.py').read())

setup(
    name='mob_suite',
    version=__version__,
    packages=find_packages(exclude=['tests']),
    url='https://github.com/jrober84/mob_suite',
    license='GPLv3',
    author='James Robertson',
    author_email='james.robertson@canada.com',
    description=('mob_suite is a set of tools for finding, typing and reconstruction of plasmids from draft and complete genome assemblies.'),
    keywords='Plasmids finding typing reconstruction',
    classifiers=classifiers,
    package_dir={'mob_suite':'mob_suite'},
    package_data={'mob_suite': ['databases/*.msh',
                            'databases/*.fas',
                            'databases/*.faa',
                            'databases/*.nhr',
                            'databases/*.nin',
                            'databases/*.nsq',
                            ]},
    install_requires=[
        'numpy>=1.11.1',
        'pandas>=0.18.1',
        'tables>=3.3.0',
        'pandas==0.22.0',
        
    ],

    entry_points={
        'console_scripts': [
            'mob_recon=mob_suite.mob_recon:main',
            'mob_cluster=mob_suite.mob_cluster:main',
            'mob_typer=mob_suite.mob_typer:main',
        ],
    },
)
