# -*- coding: utf-8 -*-

from setuptools import setup

with open('LICENSE') as f:
    license = f.read()


VERSION = '0.0.5'
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
    package_dir={'': 'pycofbuilder'},
    install_requires=['os',
                      'math',
                      'simplejason',
                      'numpy>=1.2',
                      'scipy>=1.6.3',
                      'pymatgen>=2022.0.8'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        "Intended Audience :: Science/Research",
        "License :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
