pyCOFBuilder
========================

**pyCOFBuilder** is a simple and powerfull python package to automatically assembly COF structures with specifics bulding blocks, topologies and functionalizations.

Learn more at [ToNano&Beyond](https://tonanoandbeyondblog.wordpress.com/)


## Installation

For detailed instructions see the [installation instructions](https://tonanoandbeyondblog.wordpress.com/).
If the requirements are already satisfied
```
pip install pycofbuilder
```

### Requirements
0. Python >= 3.7
1. pymatgen >= 2022.0.0


The Python dependencies are most easily satisfied using a conda
([anaconda](https://www.anaconda.com/distribution)/[miniconda](https://docs.conda.io/en/latest/miniconda.html))
installation by running

## Usage

```python
from pycofbuilder.building_block import Building_Block
from pycofbuilder.reticulum import Reticulum

BB_aldehyde = Building_Block(‘C3_BENZ_H_CHO’)
BB_amine = Building_Block(‘C2_HDZ_NH2’)

COF = Reticulum()
COF.create_hcb_a_structure(BB_aldehyde, BB_amine, stacking=‘AA’)

COF.save_cif()
COF.save_xyz(supercell=[1, 1, 2])
COF.save_qe()
COF.save_gjf()

```

For more exemples see _examples/_ and the [docs](https://github.com/lipelopesoliveira/pyCOFBuilder/examples.html)
for further examples.



## Citation

If you find **pyCOFBuilder** useful in your research please consider citing the paper:

F. L. Oliveira and P. M. Esteves,
*pyCOFBuilder: A Python Module for Automated Assembly of Covalent Organic Frameworks*

in preparation. [DOI](https://doi.org/)