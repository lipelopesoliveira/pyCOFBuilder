# -*- coding: utf-8 -*-

from setuptools import setup


with open('LICENSE') as f:
    license = f.read()


VERSION = '0.0.1'
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
    package_dir={'':'pycofbuilder'},
    install_requires=['os', 'numpy', 'scipy', 'pymatgen'], 
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
)

