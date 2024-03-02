Generating COF structures from custom Building Blocks
=====================================================

To create a COF structure from custom Building Blocks, you need the desired Building Block in the ``.xyz``, ``.gjf``, or ``.cjson`` format, as in the example below:

.. tabs::
   .. group-tab:: XYZ (.xyz)

      .. code-block:: bash

        9
        T3 building block
        C          0.9828406      0.9992759      0.0000000
        C         -0.3739780      1.3508029      0.0000000
        C         -1.3568187      0.3515270      0.0000000
        C         -0.9828406     -0.9992759      0.0000000
        C          0.3739780     -1.3508029      0.0000000
        C          1.3568187     -0.3515270      0.0000000
        C         -1.9996089     -2.0330469      0.0000000
        H         -1.8107021     -3.1273182      0.0000000
        C          2.7604747     -0.7151886      0.0000000
        H          3.6136881     -0.0044549      0.0000000
        C         -0.7608658      2.7482355      0.0000000
        H         -1.8029860      3.1317731      0.0000000
        O          1.9855845      2.0187880      0.0000000
        H          1.7077038      3.0013254      0.0000000
        O         -2.7411139      0.7101726      0.0000000
        H         -3.4530759     -0.0217478      0.0000000
        O          0.7555294     -2.7289606      0.0000000
        H          1.7453721     -2.9795776      0.0000000
        X          0.0000000      3.7112465      0.0000000
        X         -3.2140338     -1.8556233      0.0000000
        X          3.2140338     -1.8556233      0.0000000

   .. group-tab:: GaussView (.gjf)

      .. code-block:: bash

        # hf/3-21g

        T3 building block

        0 1
        C          0.9828406      0.9992759      0.0000000
        C         -0.3739780      1.3508029      0.0000000
        C         -1.3568187      0.3515270      0.0000000
        C         -0.9828406     -0.9992759      0.0000000
        C          0.3739780     -1.3508029      0.0000000
        C          1.3568187     -0.3515270      0.0000000
        C         -1.9996089     -2.0330469      0.0000000
        H         -1.8107021     -3.1273182      0.0000000
        C          2.7604747     -0.7151886      0.0000000
        H          3.6136881     -0.0044549      0.0000000
        C         -0.7608658      2.7482355      0.0000000
        H         -1.8029860      3.1317731      0.0000000
        O          1.9855845      2.0187880      0.0000000
        H          1.7077038      3.0013254      0.0000000
        O         -2.7411139      0.7101726      0.0000000
        H         -3.4530759     -0.0217478      0.0000000
        O          0.7555294     -2.7289606      0.0000000
        H          1.7453721     -2.9795776      0.0000000
        X          0.0000000      3.7112465      0.0000000
        X         -3.2140338     -1.8556233      0.0000000
        X          3.2140338     -1.8556233      0.0000000
    
    .. group-tab:: ChemicalJSON (.cjson)
      .. code-block:: json

        {
        "chemical json": 1,
        "name": "T3_BENZ_CHO_OH",
        "formula": "C9O3H6X3",
        "atoms": {
            "elements": {
                "type": [
                    "C",
                    "C",
                    "C",
                    "C",
                    "C",
                    "C",
                    "C",
                    "X",
                    "H",
                    "C",
                    "X",
                    "H",
                    "C",
                    "X",
                    "H",
                    "O",
                    "H",
                    "O",
                    "H",
                    "O",
                    "H"
                ],
                "number": [
                    6,
                    6,
                    6,
                    6,
                    6,
                    6,
                    6,
                    0.0,
                    1,
                    6,
                    0.0,
                    1,
                    6,
                    0.0,
                    1,
                    8,
                    1,
                    8,
                    1,
                    8,
                    1
                ]
            },
            "coords": {
                "3d": [
                    0.9828406,
                    0.9992759,
                    0.0,
                    -0.373978,
                    1.3508029,
                    0.0,
                    -1.3568187,
                    0.351527,
                    0.0,
                    -0.9828406,
                    -0.9992759,
                    0.0,
                    0.373978,
                    -1.3508029,
                    0.0,
                    1.3568187,
                    -0.351527,
                    0.0,
                    -1.9996089,
                    -2.0330469,
                    0.0,
                    -3.2140338,
                    -1.8556233,
                    0.0,
                    -1.8107021,
                    -3.1273182,
                    0.0,
                    2.7604747,
                    -0.7151886,
                    0.0,
                    3.2140338,
                    -1.8556233,
                    0.0,
                    3.6136881,
                    -0.0044549,
                    0.0,
                    -0.7608658,
                    2.7482355,
                    0.0,
                    0.0,
                    3.7112465,
                    0.0,
                    -1.802986,
                    3.1317731,
                    0.0,
                    1.9855845,
                    2.018788,
                    0.0,
                    1.7077038,
                    3.0013254,
                    0.0,
                    -2.7411139,
                    0.7101726,
                    0.0,
                    -3.4530759,
                    -0.0217478,
                    0.0,
                    0.7555294,
                    -2.7289606,
                    0.0,
                    1.7453721,
                    -2.9795776,
                    0.0
                ]
            }
        },
        "partialCharges": {},
        "results": [],
        "properties": {
            "totalCharge": 0,
            "spinMultiplicity": 1
            }
        }

Then, you will need to create a ``BuildingBlock`` object from the desired Building Block. Don't forget to add the ``X`` atoms on the structure where the building blocks will be connected to each other.
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