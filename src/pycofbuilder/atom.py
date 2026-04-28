from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar, Dict

import numpy as np
from numpy.typing import NDArray


def _load_periodic_table() -> Dict[str, Any]:
    """
    Loads the JSON file into memory.
    This function executes only once when the module is imported.
    """

    filename = Path(__file__).parent / "data" / "periodic_table.json"
    path = Path(filename)
    if not path.exists():
        print(f"Warning: '{filename}' not found. Atom data will be empty.")
        return {}

    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


# Load data into a "private" module-level variable.
# This ensures disk I/O happens only once.
_ATOM_DATABASE = _load_periodic_table()


@dataclass
class Atom:
    """
    Represents an atom. Data is automatically populated based on the symbol
    by querying the loaded JSON database.
    """

    # The only field required in the constructor
    symbol: str

    # Fields populated automatically (init=False)
    # The user does not provide these; they are fetched from the DB.
    full_name: str = field(init=False)
    atomic_number: int = field(init=False)
    atomic_mass: float = field(init=False)
    atomic_radius: float = field(init=False)
    covalent_radius: float = field(init=False)
    vdw_radius: float = field(init=False)
    polarizability: float = field(init=False)
    pauling_electronegativity: float = field(init=False)
    thermo_electronegativity: float = field(init=False)
    mulliken_electronegativity: float = field(init=False)
    sanderson_electronegativity: float = field(init=False)
    allen_electronegativity: float = field(init=False)
    ghosh_electronegativity: float = field(init=False)
    martynov_batsanov_electronegativity: float = field(init=False)

    # Optional field with a default value
    color: str = "#FFFFFF"

    def __post_init__(self):
        """
        Runs automatically after __init__.
        Fetches data from _ATOM_DATABASE and populates the fields.
        """
        # 1. Validation
        if self.symbol not in _ATOM_DATABASE:
            raise ValueError(
                f"Unknown atomic symbol or data not found: '{self.symbol}'"
            )

        data = _ATOM_DATABASE[self.symbol]

        # 2. Dynamic Population
        # Iterate over the dataclass fields to ensure we only grab defined attributes
        for field_name in self.__dataclass_fields__:
            if field_name in ["symbol", "color"]:
                continue  # Skip fields that are manually set or optional

            # Get value from JSON. If key is missing, defaults to None.
            value = data.get(field_name, None)

            # Assign the value to the object
            setattr(self, field_name, value)

        self.color = self._get_cpk_color()

    def _get_cpk_color(self) -> str:
        """Returns a standard CPK hex color based on the symbol.
        Colors taken from https://jmol.sourceforge.net/jscolors/
        """
        cpk_colors = {
            "H"	:"#FFFFFF",
            "He":"#D9FFFF",
            "Li":"#CC80FF",
            "Be":"#C2FF00",
            "B"	:"#FFB5B5",
            "C"	:"#909090",
            "N"	:"#3050F8",
            "O"	:"#FF0D0D",
            "F"	:"#90E050",
            "Ne":"#B3E3F5",
            "Na":"#AB5CF2",
            "Mg":"#8AFF00",
            "Al":"#BFA6A6",
            "Si":"#F0C8A0",
            "P"	:"#FF8000",
            "S"	:"#FFFF30",
            "Cl":"#1FF01F",
            "Ar":"#80D1E3",
            "K"	:"#8F40D4",
            "Ca":"#3DFF00",
            "Sc":"#E6E6E6",
            "Ti":"#BFC2C7",
            "V"	:"#A6A6AB",
            "Cr":"#8A99C7",
            "Mn":"#9C7AC7",
            "Fe":"#E06633",
            "Co":"#F090A0",
            "Ni":"#50D050",
            "Cu":"#C88033",
            "Zn":"#7D80B0",
            "Ga":"#C28F8F",
            "Ge":"#668F8F",
            "As":"#BD80E3",
            "Se":"#FFA100",
            "Br":"#A62929",
            "Kr":"#5CB8D1",
            "Rb":"#702EB0",
            "Sr":"#00FF00",
            "Y"	:"#94FFFF",
            "Zr":"#94E0E0",
            "Nb":"#73C2C9",
            "Mo":"#54B5B5",
            "Tc":"#3B9E9E",
            "Ru":"#248F8F",
            "Rh":"#0A7D8C",
            "Pd":"#006985",
            "Ag":"#C0C0C0",
            "Cd":"#FFD98F",
            "In":"#A67573",
            "Sn":"#668080",
            "Sb":"#9E63B5",
            "Te":"#D47A00",
            "I"	:"#940094",
            "Xe":"#429EB0",
            "Cs":"#57178F",
            "Ba":"#00C900",
            "La":"#70D4FF",
            "Ce":"#FFFFC7",
            "Pr":"#D9FFC7",
            "Nd":"#C7FFC7",
            "Pm":"#A3FFC7",
            "Sm":"#8FFFC7",
            "Eu":"#61FFC7",
            "Gd":"#45FFC7",
            "Tb":"#30FFC7",
            "Dy":"#1FFFC7",
            "Ho":"#00FF9C",
            "Er":"#00E675",
            "Tm":"#00D452",
            "Yb":"#00BF38",
            "Lu":"#00AB24",
            "Hf":"#4DC2FF",
            "Ta":"#4DA6FF",
            "W"	:"#2194D6",
            "Re":"#267DAB",
            "Os":"#266696",
            "Ir":"#175487",
            "Pt":"#D0D0E0",
            "Au":"#FFD123",
            "Hg":"#B8B8D0",
            "Tl":"#A6544D",
            "Pb":"#575961",
            "Bi":"#9E4FB5",
            "Po":"#AB5C00",
            "At":"#754F45",
            "Rn":"#428296",
            "Fr":"#420066",
            "Ra":"#007D00",
            "Ac":"#70ABFA",
            "Th":"#00BAFF",
            "Pa":"#00A1FF",
            "U"	:"#008FFF",
            "Np":"#0080FF",
            "Pu":"#006BFF",
            "Am":"#545CF2",
            "Cm":"#785CE3",
            "Bk":"#8A4FE3",
            "Cf":"#A136D4",
            "Es":"#B31FD4",
            "Fm":"#B31FBA",
            "Md":"#B30DA6",
            "No":"#BD0D87",
            "Lr":"#C70066",
            "Rf":"#CC0059",
            "Db":"#D1004F",
            "Sg":"#D90045",
            "Bh":"#E00038",
            "Hs":"#E6002E",
            "Mt":"#EB0026"
        }
        # Return specific color or a default 'Hot Pink' for generic/error
        return cpk_colors.get(self.symbol, "#FF1493")

    def __repr__(self):
        return f"Atom({self.symbol}: {self.full_name}, Z={self.atomic_number}, Mass={self.atomic_mass})"


@dataclass
class AtomicSite:
    """
    Represents an atom and its properties.
    """

    atom_type: str
    coordinates: NDArray[np.float64]
    label: str = ""
    properties: dict[str, Any] = field(default_factory=dict)
    index: int = -1
    atom: ClassVar[Atom] = field(init=False)

    def __post_init__(self):
        """Validation after initialization."""
        if not self.label:
            self.label = self.atom_type
        if not isinstance(self.coordinates, np.ndarray):
            self.coordinates = np.array(self.coordinates, dtype=float)

        setattr(self, "atom", Atom(self.atom_type))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, AtomicSite):
            return NotImplemented

        # Use np.allclose for float comparison safety
        coords_match = np.allclose(self.coordinates, other.coordinates, atol=1e-5)
        return (
            self.atom_type == other.atom_type
            and self.label == other.label
            and coords_match
        )

    def __repr__(self) -> str:
        x, y, z = self.coordinates
        return (
            f"Atom(type={self.atom_type}, label={self.label}, "
            f"pos=[{x:.4f}, {y:.4f}, {z:.4f}], index={self.index})"
        )

    def __add__(self, other) -> None:
        if isinstance(other, (float, int, np.ndarray)):
            self.coordinates = self.coordinates + other
        elif (
            isinstance(other, list)
            and len(other) == 3
            and all(isinstance(coord, (float, int)) for coord in other)
        ):
            self.coordinates = self.coordinates + np.array(other)
        else:
            raise TypeError(
                "Unsupported type for in-place addition with AtomicSite coordinates."
            )


    def __iadd__(self, other) -> None:
        if isinstance(other, (float, int, np.ndarray)):
            self.coordinates = self.coordinates + other
        elif (
            isinstance(other, list)
            and len(other) == 3
            and all(isinstance(coord, (float, int)) for coord in other)
        ):
            self.coordinates = self.coordinates + np.array(other)
        else:
            raise TypeError(
                "Unsupported type for in-place addition with AtomicSite coordinates."
            )
        
    def __setattr__(self, name, value) -> None:
        # Logic to intercept specific attribute changes
        if name == "atom_type":
            # If the type is changing, update the dependent 'atom' object
            # We use super() or object.__setattr__ to avoid recursion
            super().__setattr__("atom", Atom(value))
        
        elif name == "coordinates":
            # Check if the length of coordinates is 3
            if len(value) != 3:
                raise ValueError("Coordinates must be a 3-dimensional vector.")
            # Ensure coordinates are stored as a numpy array
            if not isinstance(value, np.ndarray):
                value = np.array(value, dtype=float)

        # Default behavior for all other attributes
        super().__setattr__(name, value)

    def copy(self) -> AtomicSite:
        """Create a deep copy of the AtomicSite."""
        from copy import deepcopy

        return deepcopy(self)
