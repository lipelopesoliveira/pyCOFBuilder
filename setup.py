# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()


VERSION = '0.0.1'
DESCRIPTION = 'A package for Covalent Organic Frameworks sturcture assembly.'

setup(
    name='pyCOFBuilder',
    version=VERSION,
    description=DESCRIPTION,
    long_description=readme,
    author='Felipe Lopes Oliveira',
    author_email='felipe.lopes@nano.ufrj.br',
    url='https://github.com/kennethreitz/samplemod',
    license=license,
    package_dir={'':'pycofbuilder'},
    install_requires=['os', 'numpy', 'scipy', 'pymatgen'], 
    #packages=find_packages(exclude=('tests', 'docs', 'data'))
)

