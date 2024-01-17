Installing pyCOFBuilder
=======================

Setting the environment
-----------------------

The pyCOFBuilder package is written in Python and requires a few dependencies to be installed. The easiest way to install 
these dependencies is to use the conda package manager. If you do not have conda installed, we recommend installing the 
`anaconda`_ distribution. Once conda is installed, you can create a new environment 
with the required dependencies using the following command:

.. code-block:: bash

    conda env create --file environment.yml

Requirements
------------

0. Python >= 3.10
1. pymatgen >= 2022.0.0
2. numpy >= 1.2
3. scipy >= 1.6.3
4. simplejson
5. ase
6. gemmi

.. important::
   Don't forget to activate the conda environment before installing or running pyCOFBuilder. To do so, run the following command:
   
   .. code-block:: bash
    
    conda activate pycofbuilder

Installation
------------

Option 1: Using pip
~~~~~~~~~~~~~~~~~~~

The easiest way to install ``pyCOFBuilder`` is using ``pip``. To do so, run the following command:

.. code-block:: bash

    pip install pycofbuilder

Option 2: Manually importing the module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to just download the package and manually import ``pyCOFBuilder`` it using the ``sys`` module, as exemplified below:

.. code-block:: python
    
    # importing module
    import sys
    
    # appending a path
    sys.path.append('{PATH_TO_PYCOFBUILDER}/pyCOFBuilder/src')

    import pycofbuilder as pcb

Just remember to change the ``{PATH_TO_PYCOFBUILDER}`` to the directory where you download the pyCOFBuilder package.

.. _`anaconda`: https://www.anaconda.com/download