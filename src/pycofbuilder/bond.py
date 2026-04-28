from dataclasses import dataclass

from pycofbuilder.atom import AtomicSite


@dataclass
class Bond:
    """
    Represents a chemical bond between two atoms.
    """

    atom_1: AtomicSite
    atom_2: AtomicSite
    bond_type: int
    bond_type_name: str
    bond_distance: float
