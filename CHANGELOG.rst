**********
Change Log
**********


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