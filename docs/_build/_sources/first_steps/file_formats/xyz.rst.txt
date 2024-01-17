XYZ
===

The XYZ file format is a simple text-based format commonly used to represent molecular structures in three dimensions. 
It typically contains information about the atomic coordinates of a molecule. Each line in the file corresponds to an atom 
and includes the atom's element type, as well as its X, Y, and Z coordinates.

A typical XYZ file might look like this:

.. code-block:: bash

   Number of atoms
   Comment line (optional)
   Element1 x1 y1 z1
   Element2 x2 y2 z2
   ...
   ElementN xN yN zN


- "Number of atoms" indicates the total number of atoms in the molecular structure.
- The "Comment line" is optional and is often used to provide additional information about the file.
- Each subsequent line contains the element type (e.g., C for carbon, H for hydrogen) and the X, Y, and Z coordinates of an atom.

For example:

.. code-block:: bash

   5
   Lattice="5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.0 5.44" Properties=species:S:1:pos:R:3 Time=0.0
   O 0.00000 0.00000 0.00000
   H 0.75700 0.58600 0.00000
   H -0.75700 0.58600 0.00000


This XYZ file represents a water molecule with oxygen (O) at the origin and two hydrogen (H) atoms positioned accordingly.

Lattice is a Cartesian 3x3 matrix representation of the cell lattice vectors, with each vector stored as a column and the 9 
values listed in Fortran column-major order, i.e. in the form

.. code-block:: bash

   Lattice="R1x R1y R1z R2x R2y R2z R3x R3y R3z"

where R1x R1y R1z are the Cartesian x-, y- and z-components of the first lattice vector (a), R2x R2y R2z those of the second lattice 
vector (b) and R3x R3y R3z those of the third lattice vector (c).

The list of properties in the file is described by the Properties parameter, which should take the form of a 
series of colon separated triplets giving the name, format (R for real, I for integer) and number of columns of each property. 
For example:

.. code-block:: bash

   Properties="species:S:1:pos:R:3:vel:R:3:select:I:1"

indicates the first column represents atomic species, the next three columns represent atomic positions, the next three velocities,
and the last is an single integer called select. With this property definition, the line

.. code-block:: bash

   Si        4.08000000      4.08000000      1.36000000   0.00000000      0.00000000      0.00000000       1

would describe a silicon atom at position (4.08,4.08,1.36) with zero velocity and the select property set to 1.

The extended XYZ format is now also supported by the ``ase.io.read()`` and ``ase.io.write()`` functions in the Atomic Simulation Environment (ASE) toolkit, 
and by the Ovito visualisation tool (from v2.4 beta onwards).


It's crucial to emphasize that the XYZ file format is primarily designed for representing molecular structures rather than 
crystal structures. While it efficiently captures the three-dimensional coordinates of individual atoms within a molecule, 
it lacks the necessary information to describe the periodic arrangement of atoms in a crystal lattice. Crystal structures 
involve repeating units in three dimensions, and their representation requires additional information such as unit cell parameters,
symmetry operations, and space group details. Researchers and practitioners dealing with crystallographic data should explore
dedicated file formats like CIF (Crystallographic Information File) or other formats specifically tailored for crystal structures
to ensure accurate representation and analysis of periodic arrangements in a crystal lattice.