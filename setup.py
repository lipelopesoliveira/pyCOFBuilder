# -*- coding: utf-8 -*-

from setuptools import setup

with open('LICENSE') as f:
    license = f.read()


VERSION = '0.0.8.1'
DESCRIPTION = 'A package for Covalent Organic Frameworks sturcture assembly.'

setup(
    name='pyCOFBuilder',
    version=VERSION,
    description=DESCRIPTION,
    long_description=open('README.md').read(),
    author='Felipe Lopes Oliveira',
    author_email='felipe.lopes@nano.ufrj.br',
    url='https://github.com/lipelopesoliveira/pyCOFBuilder',
    license=license,
    package_dir={'': 'src'},
    package_data={'pycofbuilder': ['pycofbuilder/data/*']},
    install_requires=['simplejson',
                      'numpy>=1.6.3',
                      'scipy>=1.2',
                      'pymatgen>=2022.0.8',
                      'jupyter',
                      'jupyterlab',
                      'pandas',
                      'tqdm',
                      'gemmi',
                      'ase',],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
