.. pycofbuilder documentation master file, created by
   sphinx-quickstart on Tue Jan 16 12:53:54 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: img/header-logo.png
   :width: 450
   :align: center
   :alt: pyCOFBuilder logo

Welcome to pyCOFBuilder's documentation!
========================================

The Python-based Covalent Organic Framework builder (pyCOFBuilder) is an open-source tool to generate 
and manipulate COF structures based on the reticular approach. The package provides a simple, lightweight, 
and fast method for creating structures with different topologies, building blocks, connection chemistry.
The package also provides  a set of tools to manipulate the generated structures, such as adding 
functional groups or changing  the stacking pattern (for 2D structures) and interpenetration degree (for 3D structures). 
The package is designed to be easily integrated with other Python-based packages for further analysis and simulation of COFs.

Here you will find tutorials and examples on how to install and use pyCOFBuilder to generate and manipulate COF 
structures. These are aimed for new users and people with more experience on molecular modeling.

First steps
-----------

Some very basic information. If you are not familiar with Python, COFs or molecular simulations, 
maybe it is best to start here, otherwise you can just skip and go to the next sections.

* :doc:`first_steps/install`

* :doc:`first_steps/simple_use`

* :doc:`first_steps/nomenclature`

* :doc:`first_steps/file_formats`

* :doc:`first_steps/custom_building_blocks`

Available building blocks
-------------------------

Here you will find a list of the building blocks that are currently available in pyCOFBuilder.

* :doc:`building_blocks/L2`
* :doc:`building_blocks/T3`
* :doc:`building_blocks/S4`
* :doc:`building_blocks/D4`
* :doc:`building_blocks/H6`
* :doc:`building_blocks/Q`
* :doc:`building_blocks/Func`

Tutorials
---------

Here you will find a list of tutorials on how to use pyCOFBuilder to generate and manipulate COF structures.

* :doc:`tutorials/bb_creation`
* :doc:`tutorials/creating_mofs`


.. toctree::
   :caption: First steps
   :hidden:
   :maxdepth: 2

   first_steps/install
   first_steps/simple_use
   first_steps/nomenclature
   first_steps/file_formats
   first_steps/visualization
   
.. toctree::
   :caption: Available building blocks
   :hidden:
   :maxdepth: 1

   building_blocks/L2
   building_blocks/T3
   building_blocks/S4
   building_blocks/H6
   building_blocks/Q
   building_blocks/Func

.. toctree::
   :caption: Tutorials
   :hidden:
   :maxdepth: 1

   tutorials/bb_creation
   tutorials/creating_mofs

.. toctree::
   :maxdepth: 2
   :caption: Code contents:

   pycofbuilder

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
