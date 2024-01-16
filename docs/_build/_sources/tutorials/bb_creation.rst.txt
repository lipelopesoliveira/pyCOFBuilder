Creation of new building blocks
===============================

All the files required to generate a building block (organic core, connection group, and functional group) can be created in the same way:

- First you need to create the structure file using you preferred software (e.g. Avogadro, GaussView, etc.)
- Then, you need to convert it to a ChemicalJSON format (`.csjon`). Although this can be done it directly with the Avogadro software, we recommend to use the `ChemJSON` module avalilable within pyCOFBuilder. This module allows conversion from several molecular formats (*e.g.* `.xyz`, `.gjf`, etc) to `.cjson` in manner that is fully compatible with the pyCOFBuilder library.
- You can add some aditional properties such as the smile code, the 4-letter code, and the xsmiles label for the structure that is being created.  
- Save the `.cjson` file in the proper folder on the `pyCOFbuilder/src/pycofbuilder/data` folder.

.. role:: raw-html(raw)
   :format: html

.. important::
   Don't forget to add the special points (Q, X, or :raw-html:`R<sub>y</sub>`.) in the structure of the molecule you want to add as in the image below!
   
   

For organic cores by default the distance ``Q-C`` should be set to 0.5 angstroms and the distance ``Ry-C`` should be set to 1.0 angstrom.

For the connection groups and functional groups by default the distance from the connection points, ``X`` and ``R`` respectivelly, should be set in a way to get the final espected distance between the connection points of the building block and the linkers.

Below there is an example of how to create a new building block using the ``ChemicalJSON`` module.

.. code-block:: python
    
    from pycofbuilder.cjson import ChemJSON

    # Create an empty ChemJSON object
    new_BB = ChemJSON()

    # Read the file conatining the molecular structure
    new_BB.from_gjf(os.getcwd(), 'L2_BENZ.gjf')

    # Define the name of the molecule. This is just a label and it is not used in the creation of the building block
    new_BB.name = 'benzene'

    # Define the properties of the building block. Although this informations are not required to create the building block, it is 
    # recomended to add them in order to take full advantage of the pyCOFBuilder capabilities
    new_BB.properties = {
        "smiles": "[Q]C1=C([R2])C([R1])=C([Q])C([R2])=C1[R1]",
        "code": "BENZ",
        "xsmiles": "[*]C1=C([*])C([*])=C([*])C([*])=C1[*]",
        "xsmiles_label": "|$Q;;;R2;;R1;;Q;;R2;;R1$|",
    }

    # Save the ChemJSON object as a cjson file
    new_BB.write_cjson_file(os.getcwd(), 'BENZ.cjson')

