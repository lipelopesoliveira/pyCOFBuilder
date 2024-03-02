Protein Data Bank with Charges and Radii (PQR)
==============================================

The PQR (PQR) file format is a variant of the Protein Data Bank (PDB) format that includes additional 
information about the atomic charges and radii. PQR files are often used in molecular dynamics simulations, 
particularly with software like APBS (Adaptive Poisson-Boltzmann Solver), to include information about 
the electrostatic properties of biomolecular structures. Below is a simplified example of a PQR file:

.. code-block:: bash

   HEADER    SMALL MOLECULE CRYSTAL STRUCTURE
   CRYST1    5.0   5.0   5.0  90.00  90.00  90.00 P 1           1
   ATOM      1  O   WAT     1       0.000   0.000   0.000 -0.600  1.4
   ATOM      2  H1  WAT     1       0.757   0.586   0.000  0.300  1.2
   ATOM      3  H2  WAT     1      -0.757   0.586   0.000  0.300  1.2
   END

In this example:

- `ATOM` lines represent individual atoms in the structure.
- Columns contain information about atom number, atom type, residue name, residue number, X, Y, and Z coordinates, partial charge, and atomic radius.
- The partial charge (in this case, -0.600 for oxygen and 0.300 for hydrogen) is typically represented in units of electron charge (e).

The PQR file format is an extension of the PDB format, where the last two columns are added for charge and radius information. 
This allows users to include electrostatic information for molecular dynamics simulations or other calculations that involve 
the consideration of atomic charges and radii. The specific details of the PQR format may vary depending on the software or tool being used,
so it's essential to consult the documentation of the specific application for accurate information.


