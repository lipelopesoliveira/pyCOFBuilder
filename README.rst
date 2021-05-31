pyCOFBuilder
========================

**pyCOFBuilder** is a simple and powerfull python package to automatically assembly COF structures with specifics bulding blocks, topologies and functionalizations.

`Learn more <https://tonanoandbeyondblog.wordpress.com/>`_.


## Installation

For detailed instructions see the [installation instructions](https://tonanoandbeyondblog.wordpress.com/).
If the requirements are already satisfied
```
pip install cgbind
```

### Requirements
0. Python v. 3.7
1. pymatgen


The Python dependencies are most easily satisfied using a conda
([anaconda](https://www.anaconda.com/distribution)/[miniconda](https://docs.conda.io/en/latest/miniconda.html))
installation by running

## Usage

```python
from pyCOFBuilder import Building_Block, Reticulum

BB_aldehyde = Building_Block(‘C3_BENZ_H_CHO’)
BB_amine = Building_Block(‘C2_HDZ_NH2’)

COF = Reticulum()
COF.create_hcb_a_structure(BB_aldehyde, BB_amine, stacking=‘AA’)

COF.save_cif()
COF.save_xyz(supercell=[1, 1, 2])
COF.save_qe()
COF.save_gjf()

```


## Citation

If you find **pyCOFBuilder** useful in your research please consider citing the paper:

F. L. Oliveira and P. M. Esteves,
*pyCOFBuilder: A Python Module for Automated Assembly of Covalent Organic Frameworks*
in preparation. [DOI](https://doi.org/)