Generating COF structures from custom Building Blocks
=====================================================

To create a COF structure from custom Building Blocks, you need the desired Building Block in the ``.xyz``, ``.gjf``, or ``.cjson`` format, as in the example below:

.. tabs::
   .. group-tab:: ChemicalJSON (.cjson)
      .. code-block:: json

         {
         "chemical json": 1,
         "name": "benzene",
         "formula": "Q2 C6 R22 R12",
         "atoms": {
            "elements": {
                  "type": [
                     "C",
                     "C",
                     "C",
                     "C",
                     "C",
                     "C",
                     "R1",
                     "R1",
                     "R2",
                     "R2",
                     "Q",
                     "Q"
                  ],
                  "number": [
                     6,
                     6,
                     6,
                     6,
                     6,
                     6,
                     0.0,
                     0.0,
                     0.0,
                     0.0,
                     0.0,
                     0.0
                  ]
            },
            "coords": {
                  "3d": [
                     1.21383527,
                     0.70080812,
                     0.0,
                     0.0,
                     1.40161624,
                     0.0,
                     -1.21383527,
                     0.70080812,
                     0.0,
                     -1.21383527,
                     -0.70080812,
                     0.0,
                     0.0,
                     -1.40161624,
                     0.0,
                     1.21383527,
                     -0.70080812,
                     0.0,
                     2.07986067,
                     1.20080812,
                     0.0,
                     -2.07986067,
                     -1.20080812,
                     0.0,
                     2.07986067,
                     -1.20080812,
                     0.0,
                     -2.07986067,
                     1.20080812,
                     0.0,
                     0.0,
                     2.40161624,
                     0.0,
                     0.0,
                     -2.40161624,
                     0.0
                  ]
            }
         },
         "properties": {
            "smiles": "[Q]C1=C([R2])C([R1])=C([Q])C([R2])=C1[R1]",
            "code": "BENZ",
            "xsmiles": "[*]C1=C([*])C([*])=C([*])C([*])=C1[*]",
            "xsmiles_label": "|$Q;;;R2;;R1;;Q;;R2;;R1$|"
         }
         }
   .. group-tab:: XYZ (.xyz)

      .. code-block:: bash

         12
         L2_BENZ building block
         C                  1.21383527    0.70080812    0.00000000
         C                  0.00000000    1.40161624    0.00000000
         C                 -1.21383527    0.70080812    0.00000000
         C                 -1.21383527   -0.70080812    0.00000000
         C                  0.00000000   -1.40161624    0.00000000
         C                  1.21383527   -0.70080812    0.00000000
         R1                 2.07986067    1.20080812    0.00000000
         R1                -2.07986067   -1.20080812    0.00000000
         R2                 2.07986067   -1.20080812    0.00000000
         R2                -2.07986067    1.20080812    0.00000000
         Q                  0.00000000    2.40161624    0.00000000
         Q                  0.00000000   -2.40161624    0.00000000

   .. group-tab:: GaussView (.gjf)

      .. code-block:: bash

         # hf/3-21g

         L2_BENZ building block

         0 1
         C                  1.21383527    0.70080812    0.00000000
         C                  0.00000000    1.40161624    0.00000000
         C                 -1.21383527    0.70080812    0.00000000
         C                 -1.21383527   -0.70080812    0.00000000
         C                  0.00000000   -1.40161624    0.00000000
         C                  1.21383527   -0.70080812    0.00000000
         R1                 2.07986067    1.20080812    0.00000000
         R1                -2.07986067   -1.20080812    0.00000000
         R2                 2.07986067   -1.20080812    0.00000000
         R2                -2.07986067    1.20080812    0.00000000
         Q                  0.00000000    2.40161624    0.00000000
         Q                  0.00000000   -2.40161624    0.00000000

Then, you will need to create a ``BuildingBlock`` object from the desired Building Block. The ``BuildingBlock`` class is used to create the COF structure. 
This method requires the path to the Building Block file and the name of the Building Block file and the type of connector. 
The connector atom is the atom that will be used to determine the correct atom to connect the Building Blocks. 
The method to create the ``BuildingBlock`` object is shown below:

.. code-block:: python

    import pycofbuilder as pcb

    building_block = pcb.BuildingBlock(
        name='BB_T3.xyz',
        out_dir='.',
        conector='NH2'
        )


Then, you can use the ``from_building_blocks`` method to generate the COF structure. 
The method requires the following parameters: the Building Blocks, the net type and the stacking/interpenetration degree.
The method to create the COF structure is shown below:

.. code-block:: python

    import pycofbuilder as pcb

    # Create the building blocks
    bb_t3 = pcb.BuildingBlock('BB_T3.xyz', out_dir='.', conector='NH2')
    bb_l2 = pcb.BuildingBlock('BB_L2.xyz', out_dir='.', conector='CHO')

    # Create the structure
    cof = pcb.Framework()
    cof.from_building_blocks(bb_t3, bb_l2, 'HCB_A', 'AA')