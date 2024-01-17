COFs and Building Blocks nomenclature
=====================================

COF nomenclature
----------------

In order to ensure greater reproducibility as well as quickly and easily access to relevant information from the COFs, I've developed a simple nomenclature to name this structures. 
Generally speaking, a COF can be described as building blocks connected by covalent bonds following a underlying net and presenting a specific stacking pattern/interpenetration degree. 
In this way, the COF name is represented as

``BuildingBlock1-BuildingBlock2-Net-Stacking``

pyCOFBuilder expects the COF name to be in this fixed format, with each of these characterists represented by "cards" separated by a dash (``-``). The "cards" are described as follows:

- ``BuildingBlock1``: The building block with the greater connectivity.
- ``BuildingBlock2``: The building block with the smaller connectivity.
- ``Net``: The net describing the reticular structure.
- ``Stacking``: The stacking (for 2D structures) or interpenetrating degree (for 3D structures)

Building blocks
---------------

To name the building blocks we also developed a set of rules. The building block can be described as

``Symmetry_Core_Connector_FunctionalGroupR1_FunctionalGroupR2_FunctionalGroupR3_...``

where:

- ``Symmetry``: The general symmetry of the building block. Also represents the connectivity of the building block. For 2D building blocks can be ``L2``, ``T3`` or ``S4``, and ``H6``.
- ``Core``: The 4 letters code referring to the building block core.
- ``Connector``: The type of functional group that will be used to assembly the COF structure. Ex.: ``NH2``, ``CHO``, ``CONHNH2``, etc.
- ``FunctionalGroupRN``: The Nth functional group in the structure. The number of Functional groups will change according to the availability of the core.

Note that every "card" for the building block name is separated by an underline (``_``) and every "card" for the COF name is separated by a dash (``-``). 
This makes it easy to split the COF name into useful information.

Symmetry
~~~~~~~~

The symmetry of the building block is represented by a code composed of a letter and a number. The letter represents the geometric 
figure of the building block and the number represents the connectivity of the building block. The table below shows the symbols,
connectivity numbers, and geometric figures used to represent the building blocks. *Please note that not all elements of the table are
currently implemented in pyCOFBuilder.*

.. list-table:: Symbols, connectivity numbers, and geometric figures used to represent the building blocks.
   :widths: 25 25 50
   :align: center
   :header-rows: 1

   * - Symbol
     - Connectivity
     - Geometric Figure
   * - L
     - 2
     - Line
   * - T
     - 3
     - Triangle
   * - S
     - 4
     - Square
   * - R
     - 4
     - Rectangle
   * - T
     - 4
     - Tetrahedron
   * - O
     - 6
     - Octahedron
   * - P
     - 6
     - Trigonal prism
   * - H
     - 6
     - Hexagon
   * - C
     - 8
     - Cube
   * - A
     - 8
     - Square antiprism
   * - E
     - 8
     - Octagon
   * - B
     - 12
     - Cuboctahedron
   * - I
     - 12
     - Icosahedron
   * - U
     - 12
     - Truncated tetrahedron
   * - X
     - 12
     - Hexagon prism


Organic Cores
~~~~~~~~~~~~~

The organic cores are denoted by a four-letter code designed to closely reflect the IUPAC name of the equivalent isolated molecule.
For instance, both the building blocks 1,3,5-triformylphloroglucinol and 1,4-diaminobenzene share the same organic core, 
representing a benzene molecule. Consequently, the code assigned to this organic core is ``BENZ``.

Although constructing a consistent representation for the majority of COFs found in the literature is relatively straightforward, 
certain types of covalent connections can pose a greater challenge due to their impact on the formation of the final structure. 
For example, COFs formed through boronic condensation between a boronic diacid, such as COF-1 `cote2005porous`_, 
feature only a linear building block, which directly constrains the way the string-based representation is constructed. 
However, these materials can be perceived as being formed by an ``HCB_A`` network, wherein the tritopic building block corresponds 
to the hexagon formed by O and B atoms (``BRXN``), while the linear ditopic block represents the benzene ring (``BENZ``) 
present in the original diboronic acid, as illustrated in the figure below. Consequently, this structure can be represented 
as ``T3_BRXN_BOH2-L2_BENZ_BOH2_H_H_HCB_A-AB``. 

The identical approach can be applied to COFs synthesized through the trimerization of dicyanobenzene, resulting in the formation 
of the triazine-based framework CTF-1. `kuhn2008porous`_ In this case, the structure can be constructed using the representation 
``T3_TRZN_CN-L2_BENZ_CN_H_H_HCB_A-AA``. 

.. image:: ../img/SBU.png
   :width: 700
   :align: center
   :alt: SBU approach to determine the organic core of a building block

Another challenge of this nomenclature approach arises when dealing with the representation of multicomponent structures. `huang2016multiple`_ 
These structures are composed of three or more building blocks that can occupy the same reticular site, preventing a straightforward decomposition 
into one of the existing topological networks. Thus, in the current version of pyCOFBuilder, such structures cannot be constructed.



Connection groups and Functional Groups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For connector groups and functional groups, which have fewer atoms and are chemically simpler, their composition is used for representation. 
For instance, the amine group, which can function as both a connector and a functional group, is represented by the code ``NH2``. 
The aldehyde group is denoted by the code ``CHO``, and so forth.

Nets
----

The nets are described a three letter code representing the underlying net adapted from the Reticular Chemistry Structure Resource (`RCSR`_).
The augmented nets are represented by the three letter code followed by a ``_A``. The nets are described as below

.. _`RCSR`: https://rcsr.anu.edu.au/

.. image:: ../img/2D_nets.png
   :width: 700
   :align: center
   :alt: pyCOFBuilder 2D nets


Stacking / Interpenetration
---------------------------

To represent the stacking pattern or interpenetrating class desired, the encoding string will depend on the selected net.
For 2D nets, available stacking patterns include ``AA``, ``AB1``, ``AB2``, ``AAl``, and ``AAt``, among others deppending on the topology. 
In the case of 3D nets, the encoding string is determined by the number of interpenetrating structures.

.. image:: ../img/staking_hex.png
   :width: 700
   :align: center
   :alt: Stacking patterns for hexagonal 2D nets


.. image:: ../img/staking_square.png
   :width: 700
   :align: center
   :alt: Stacking patterns for square 2D nets


For the interpenetration degree of 3D nets, the encoding string is determined by the number of interpenetrating structures.


.. rubric:: References

.. [cote2005porous] Cote, A.P.; Benin, A.I.; Ockwig, N.W.; O’Keeffe, M.; Matzger, A.J.; and Yaghi, O.M. “Porous, crystalline, covalent organic frameworks.” science, 2005. 310(5751):1166–1170

.. [kuhn2008porous] Kuhn, P.; Antonietti, M.; and Thomas, A. “Porous, covalent triazine-based frameworks prepared by ionothermal synthesis.” Angewandte Chemie International Edition, 2008. 47(18):3450–3453

.. [huang2016multiple] Huang, N.; Zhai, L.; Coupry, D.E.; Addicoat, M.A.; Okushita, K.; Nishimura, K.; Heine, T.; and Jiang, D. “Multiple-component covalent organic frameworks.” Nature communications, 2016. 7(1):12325

