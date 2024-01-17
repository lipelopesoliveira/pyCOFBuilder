Turbomole (coord)
=================

Turbomole is a quantum chemistry software suite that is widely used for electronic structure calculations. 
The Turbomole input file, often named `control`, plays a crucial role in defining the system, setting up the calculation,
and specifying various parameters. Below is a simplified example of a Turbomole input file for a molecular structure:

.. code-block:: bash

   $coord
      0.000000   0.000000   0.000000  O
      0.757000   0.586000   0.000000  H
     -0.757000   0.586000   0.000000  H
   $end

   $basis
      * basis set details for each element *
   $end

   $scf
      maxiter      200
      conver      1.0e-6
   $end


In this example:

- The `$coord` block specifies the atomic coordinates and element types. In this case, it represents a water molecule with oxygen (O) and hydrogen (H) atoms.
- The `$basis` block is typically used to define the basis set details for each element. In this simplified example, it is omitted, but in a real input file, you would specify the basis set for each atom.
- The `$scf` block contains parameters related to the self-consistent field (SCF) convergence, such as the maximum number of iterations (`maxiter`) and the convergence threshold (`conver`).

This example is quite basic, and a complete Turbomole input file would include additional sections specifying details like the type of calculation 
(e.g., geometry optimization, single-point energy), molecular symmetry, and various numerical parameters. The specific details and format may 
vary depending on the version of Turbomole and the type of calculation you are performing. Always consult the Turbomole documentation for accurate 
and up-to-date information on the input file format and parameters.
