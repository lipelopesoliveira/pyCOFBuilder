**********
Change Log
**********

v0.0.8.8 Beta (01/02/2025)
========================

New features
------------

- Refactor `IO_Tools` module to improve the code quality and readability.
- Add unit tests with pytest for the modules:
  - `IO_Tools`
  - `Tools`
  - `CJSON`
  - `Framework`
  - `BuildingBlock`
- Add the LON and LON-A nets.
- Add the FXT and FXT-A nets.
- Add the `geometry` module to handle the geometry identification of molecular structures.
- Now the building blocks are not save by default. The user can activate the save of the building blocks by setting the `save_bb` variable to `True` on the `Framework` class.

Bugfixes
--------

- Fix bond types not following the CCDC conventions.
- Fix a bug where the folder for building blocks was created even if the user did not want to save the building blocks.
- Fix a bug where the save_dir was not created if it did not exist.
- Fix the `out_path` variable on the `Framework` class to save the structure on the correct folder.


Known bugs
----------

- The `geometry` module is not fully implemented yet.
- The calculation of bonds is taking too long for large structures.
- The bonds calculation does not identify the bond types.
- The bond assignment is done after the structure is created, which can lead to errors in the structure creation. The bond assignment should be done before the structure creation, and the bonds betweeen building blocks should be assigned manually on the creation of the structure.
- The positioning of the building blocks on FXT and FXT-A nets is still not perfect.
- On FXT and FXT-A nets the cell parameters for R4 building blocks are not being calculated correctly.
- The BOR network is not working properly.
- The LON-A network is not implemented yet.
- Some S4 and R4 building blocks are not being positioned correctly on the structure.

v0.0.8.7 Beta (29/12/2024)
========================

New features
------------

- Added new connection groups:
  - CCH3O
  
- Added new building blocks:
  - L2: PRZN2

Bugfixes
--------

- Remove the `pwd`` library to ensure compatibility with Windows systems.
  
v0.0.8.6 Beta (03/09/2024)
========================

New features
------------

- Added new connection groups:
  - CCH3O
  
- Added new building blocks:
  - L2: PRZN

v0.0.8.5 Beta (16/04/2024)
========================

New features
------------

- Added new building blocks:
  - L2: DFDB
  - R4: TPDT

Bugfixes
--------

- Fixed the bug that prevent the generation of `SQL` and `SQL_A` nets with `R4` building block.


v0.0.8.4 Beta (16/04/2024)
========================

New features
------------

The new method Framework.make_cubic() generates a cubic super-cell with the structure, allowing the creation of membranes or perpendicular boxes easier. The result is not necessarily a cubic unit cell, but a cell with angles equal to 90 and cell parameters as close as possible to cubic.

- Added new building blocks:

  - L2
    - NAP2
  - S4
    - OTPR
    - TBPR
  - R4
    - ATTP
    - PRLN
    - TPLN
    - ETKB

Bugfixes
--------

- Fixed a bug in the creation of the `DIA` and `DIA_A` nets.
- Fix the position of the atoms on `CH2CN` and `CHO` connection groups.

v0.0.8.7 Beta (01/02/2025)
========================

Bugfixes
--------

- Remove the `pwd` library to ensure compatibility with Windows systems.


v0.0.8.3 Beta (16/04/2024)
========================

New features
------------

- Added new building blocks:
  - L2: DFDB
  - R4: TPDT

Bugfixes
--------

- Fixed the bug that prevent the generation of `SQL` and `SQL_A` nets with `R4` building block.


v0.0.6 Beta (02/03/2024)
========================

New features
------------

- A web-based documentation of pyCOFBuilder, as a result of #51
- Possibility to create 3D nets with `DIA` and `DIA_A` topology as a result of #54 
- Possibility to create 3D nets with `BOR` topology as a result of #54 
- Add new D4 organic cores:
  - ADAM
  - SBFE
  - TDAT
  - TKAT
  - TKPM
- Add new custom exceptions:
  - `BondLenghError` exception that is raised when the distance between two atoms on the structure are smaller than a distance thresshold. It is controlled by the `dist_threshold` variable on the `Framework` class (0.8 angstrom by default)
  - `BBConnectivityError` exception raised when the building block connectivity is not valid.
  - `ConnectionGroupError` exception raised when the connection group is not valid.
  - `MissingXError` exception raised when the custom building block is missing X atoms.
- The `CJSON` module now has the capability to read and write results from simulations. 
- Add the possibility to create MOF structures
- Add a new log system that can print on the screen or save on a file the log.

Bugfixes
--------

- It's now much easier to create and use custom building blocks.
- HXL-A and KDG are working properly now.

v0.0.2 Beta (17/06/2022)
========================

Added 
-----

- Add AA, AB1, AB2, AAl, AAt, ABC1 e ABC2 stakings for KDG net https://github.com/lipelopesoliveira/pyCOFBuilder/pull/23
- Add a new C6 HEXB buinding block derived from `hexaphenilbenzene <https://en.wikipedia.org/wiki/Hexaphenylbenzene>`__ https://github.com/lipelopesoliveira/pyCOFBuilder/pull/23 
- Code for creation of C6 building block https://github.com/lipelopesoliveira/pyCOFBuilder/pull/23
- AA, AB1, AB2, AAl, AAt, ABC1 e ABC2 stakings for HXL-A net https://github.com/lipelopesoliveira/pyCOFBuilder/pull/24
- Add a new C4 Buiding Block derived from 4,4',4'',4'''-(pyrene-1,3,6,8-tetrayl)tetrabenzene. https://github.com/lipelopesoliveira/pyCOFBuilder/pull/25
- Add AA, AB1x, AB1y, AB1xy, AB2, AAl, AAt, stakings for KGM and KGM-A net https://github.com/lipelopesoliveira/pyCOFBuilder/issues/18
- Add proper documentaion of the net methods 
- Add the class methods documentations

Know bugs
---------

- KGM and KGM-A nets do not generate the proper structure
- HXL-A and KDG stakings are not tested

v0.0.1 Alpha (09/06/2022)
=========================

Added
-----

- General structure of the code
- COF generation with HCB and HCB-A nets
- AA, AB1, AB2, AAl, AAt, ABC1 e ABC2 stakings for HCB and HCB-A nets
- Several types of organic cores, functional groups and conectors