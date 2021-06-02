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

To create a specific COF, like C3_BENZ_CHO_H-C2_HDZ_NH2-HCB-A-AA:
```python
import pycofbuilder as COF

COF.build('C3_BENZ_CHO_H-C2_HDZ_NH2-HCB-A-AA')

```

A `.cif` file will be created in the `out` folder. 

For more exemples see _examples/_ and the [docs](https://github.com/lipelopesoliveira/pyCOFBuilder/examples.html)
for further examples.



## Citation

If you find **pyCOFBuilder** useful in your research please consider citing the paper:

F. L. Oliveira and P. M. Esteves,
*pyCOFBuilder: A Python Module for Automated Assembly of Covalent Organic Frameworks*

*in preparation.* [DOI](https://doi.org/)