Chemical JSON (CJSON)
=======================

Them Chemical JSON format is a adaptiation of the JSON file format to encode chemical information developed and mantained by the
OpenChemistry project on the `repository`_. The format is designed to be human readable and easy to parse by computers. The format is designed to be
extensible, allowing for the addition of new properties and features without breaking existing code. 

.. _`repository`: https://github.com/OpenChemistry/chemicaljson

Below is a simplified example of a ``cjson`` file representing a crystal structure of Rutile:

.. code-block:: python

   {
   "chemical json": 0,
   "name": "TiO2 rutile",
   "formula": "Ti 2 O 4",
   "unit cell": {
      "a": 2.95812,
      "b": 4.59373,
      "c": 4.59373,
      "alpha": 90.0,
      "beta":  90.0,
      "gamma": 90.0
   },
   "atoms": {
      "elements": {
         "number": [ 22, 22, 8, 8, 8, 8 ]
      },
      "coords": {
         "3d fractional": [ 0.00000, 0.00000, 0.00000,
                           0.50000, 0.50000, 0.50000,
                           0.00000, 0.30530, 0.30530,
                           0.00000, 0.69470, 0.69470,
                           0.50000, 0.19470, 0.80530,
                           0.50000, 0.80530, 0.19470 ]
      }
   }
   }