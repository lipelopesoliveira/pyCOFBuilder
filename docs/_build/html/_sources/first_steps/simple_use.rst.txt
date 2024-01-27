Generating COF structures
=========================

To create a specific COF, such as ``T3_BENZ_NH2_OH-L2_BENZ_CHO_H-HCB_A-AA``:

.. code-block:: python

    import pycofbuilder as pcb

    cof = pcb.Framework('T3_BENZ_CHO_OH-L2_BENZ_NH2_H-HCB_A-AA')
    cof.save(fmt='cif', supercell = [1, 1, 2], save_dir = '.')


You should see an output such as:

.. code-block:: python
    
    T3_BENZ_NH2_OH-L2_BENZ_CHO_H_H-HCB_A-AA                       hexagonal   P    P6/m # 175    12 sym. op.

A ``.cif`` file (the default save format is CIF, but it can be easily changed by setting other value on the ``fmt`` option) will be created in the ``out`` folder. 
The code will print out some information about the structure created. You can turn off this information output by setting the ``silente`` option to ``True``:

.. code-block:: python

    cof = pcb.Framework('T3_BENZ_CHO_OH-L2_BENZ_NH2_H-HCB_A-AA', silent = True)


Currently, it is possible to select the following formats:

- ``cif``
- ``xsf``
- ``pdb``
- ``cjson``
- ``vasp``
- ``turbomole``
- ``pqr``
- ``qe``
- ``gjf``
- ``xyz``
  
For more details see :doc:`file_formats`.


Besides, the variable ``cof`` now is a ``Framework`` object. This object has some attributes that can be accessed:

.. code-block:: python

    cof.name
    'T3_BENZ_NH2_OH-L2_BENZ_CHO_H-HCB_A-AA'
    cof.smiles
    '(N)C1=C(O)C((N))=C(O)C((N))=C1O.(C([H])=O)C1=C([H])C([H])=C((C([H])=O))C([H])=C1[H]'
    cof.lattice
    array([[ 22.49540055,   0.        ,   0.        ],
           [-11.24770028,  19.48158835,   0.        ],
           [  0.        ,   0.        ,   3.6       ]])
    cof.n_atoms
    72
    cof.space_group
    'P6/m'

