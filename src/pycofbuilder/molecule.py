from __future__ import annotations

import json
import os
from abc import abstractmethod
from pathlib import Path
from typing import Iterator, Sequence

import numpy as np
from numpy.typing import NDArray

from pycofbuilder.atom import AtomicSite
from pycofbuilder.bond import Bond
from pycofbuilder.utils import read_mol_file, rotation_matrix_from_vectors


class Molecule:
    def __init__(
        self,
        atomTypes: Sequence[str] = [],
        cartesianCoords: Sequence[tuple[float, float, float]] = [],
        atomLabels: Sequence[str] | None = None,
        siteProperties: list[dict] | None = None,
        bonds: Sequence[tuple[int, int]] | None = None,
        bondTypes: Sequence[int] | None = None,
        properties: dict | None = None,
    ) -> None:

        self.sites = []

        for i in range(len(atomTypes)):
            self.sites.append(
                AtomicSite(
                    atom_type=atomTypes[i],
                    label=atomLabels[i] if atomLabels else atomTypes[i],
                    coordinates=np.array(cartesianCoords[i]),
                    properties=siteProperties[i] if siteProperties else {},
                    index=i,
                )
            )

        self.bonds = []
        if bonds and bondTypes:
            self.set_bonds(bonds, bondTypes)

        self._charge = properties.get("charge", 0) if properties else 0
        self._spin_multiplicity = (
            properties.get("spin_multiplicity", 1)
            if properties
            else self.guess_multiplicity()
        )

        self.properties = properties or {
            "charge": self.charge,
            "spin_multiplicity": self.spin_multiplicity,
        }

        self.check_partial_charges()

    def __len__(self) -> int:
        return len(self.sites)

    def __iter__(self) -> Iterator[AtomicSite]:
        return iter(self.sites)

    def __setitem__(self, idx: int, value: AtomicSite) -> None:
        self.sites[idx] = value

    def __getitem__(self, idx: int) -> AtomicSite:
        return self.sites[idx]

    def __delitem__(self, idx: int) -> None:
        del self.sites[idx]

    def __add__(self, to_add) -> None:
        if isinstance(to_add, (int, float, np.ndarray)):
            for site in self.sites:
                site.coordinates = site.coordinates + to_add

        elif (
            isinstance(to_add, list)
            and len(to_add) == 3
            and all(isinstance(coord, (int, float)) for coord in to_add)
        ):
            for site in self.sites:
                site.coordinates = site.coordinates + np.array(to_add)

        elif isinstance(to_add, Molecule):
            original_size = len(self.sites)
            for i, site in enumerate(to_add.sites):
                site.index = original_size + i

            self.sites.extend(to_add)
            self.bonds.extend(to_add.bonds)

    def __subtract__(self, to_subtract) -> None:
        if isinstance(to_subtract, (int, float, np.ndarray)):
            for site in self.sites:
                site.coordinates = site.coordinates - to_subtract
        elif (
            isinstance(to_subtract, list)
            and len(to_subtract) == 3
            and all(isinstance(coord, float) for coord in to_subtract)
        ):
            for site in self.sites:
                site.coordinates = site.coordinates - np.array(to_subtract)

    def __repr__(self) -> str:
        return "Molecule Summary\n" + "\n".join(map(repr, self))

    def __str__(self) -> str:
        outs = [
            f"Full Formula ({self.composition})",
            f"Reduced Formula: {self.composition}",
            f"Charge = {self._charge}, Spin Mult = {self.spin_multiplicity}",
            f"Sites ({len(self)}), Bonds ({len(self.bonds)})",
        ]
        for site in self:
            outs.append(
                f"{site.index:4} {site.atom_type:3} {' '.join([f'{site.coordinates[coord]:0.6f}'.rjust(12) for coord in range(3)])}"
            )

        outs += ["Bonds:"]

        for bond in self.bonds:
            outs.append(
                "{}_{} - {}_{}, {} ({}), {:0.4f} Å".format(
                    bond.atom_1.atom_type,
                    bond.atom_1.index,
                    bond.atom_2.atom_type,
                    bond.atom_2.index,
                    bond.bond_type_name,
                    bond.bond_type,
                    bond.bond_distance,
                )
            )
        return "\n".join(outs)

    def copy(self):
        """Create a deep copy of the molecule."""
        from copy import deepcopy

        return deepcopy(self)

    def check_partial_charges(self):
        # Check if the difference on charge and sum of partial charges is significant
        if abs(sum(self.partial_charges) - self.charge) > 1e-3:
            raise Warning(
                f"Total partial charge ({sum(self.partial_charges)}) does not match molecular charge ({self.charge})."
            )

    def remove_bonds(self, atom: AtomicSite) -> None:
        """Remove all bonds associated with a given atom."""
        self.bonds = [
            bond
            for bond in self.bonds
            if bond.atom_1 is not atom and bond.atom_2 is not atom
        ]

    def remove_atom(self, idx: int) -> None:
        self.remove_bonds(self.sites[idx])
        del self.sites[idx]
        self.reset_index()

    def reset_index(self) -> None:
        for i, atom in enumerate(self.sites):
            atom.index = i

    @property
    def charge(self) -> float:
        """Charge of molecule."""
        return self._charge

    @charge.setter
    def charge(self, value: float) -> None:
        """Set the charge of the molecule."""
        self._charge = value

    @property
    def spin_multiplicity(self) -> int:
        """Spin multiplicity of molecule."""
        return self._spin_multiplicity

    @spin_multiplicity.setter
    def spin_multiplicity(self, value: int) -> None:
        """Set the spin multiplicity of the molecule."""
        self._spin_multiplicity = value

    @property
    def nelectrons(self) -> float:
        """Number of electrons in the molecule."""
        n_electrons = 0.0
        for site in self:
            if site.atom is not None:
                n_electrons += site.atom.atomic_number
        n_electrons -= self.charge
        return n_electrons

    @abstractmethod
    def guess_multiplicity(self) -> int:
        """Guess the spin multiplicity based on the number of electrons."""
        n_electrons = self.nelectrons
        return 1 if n_electrons % 2 == 0 and n_electrons > 1 else 2

    @property
    def positions(self) -> NDArray:
        """Get the positions of all atoms in the molecule as a numpy array."""
        return np.array([site.coordinates for site in self])

    @positions.setter
    def positions(self, new_positions: NDArray) -> None:
        """Set the positions of all atoms in the molecule."""
        if new_positions.shape != (len(self), 3):
            raise ValueError("New positions must have shape (n_atoms, 3).")
        for i, site in enumerate(self):
            site.coordinates = new_positions[i]

    @property
    def labels(self) -> list[str]:
        """Labels of the atoms in the molecule."""
        return [site.label for site in self]

    @labels.setter
    def labels(self, new_labels: list[str]) -> None:
        """Set the labels of the atoms in the molecule."""
        if len(new_labels) != len(self):
            raise ValueError(
                "Length of new_labels must match number of atoms in the molecule."
            )
        for i, site in enumerate(self):
            site.label = new_labels[i]

    @property
    def atom_types(self) -> list[str]:
        """Atom types of the molecule."""
        return [site.atom_type for site in self]

    @atom_types.setter
    def atom_types(self, new_types: list[str]) -> None:
        """Set the atom types of the molecule."""
        if len(new_types) != len(self):
            raise ValueError(
                "Length of new_types must match number of atoms in the molecule."
            )
        for i, site in enumerate(self):
            site.atom_type = new_types[i]

    @property
    def composition(self):
        """Composition of the molecule."""
        from collections import Counter

        comp_dict = Counter([site.atom_type for site in self])
        # Convert to string
        comp = "".join([f"{el}{comp_dict[el]}" for el in sorted(comp_dict.keys())])
        return comp

    @property
    def partial_charges(self) -> list[float]:
        """List of partial charges for each atom in the molecule."""
        return [site.properties.get("partial_charge", 0.0) for site in self]

    @partial_charges.setter
    def partial_charges(self, charges: list[float]) -> None:
        """Set the partial charges for each atom in the molecule."""
        if len(charges) != len(self):
            raise ValueError(
                "Length of charges must match number of atoms in the molecule."
            )
        for i, site in enumerate(self):
            site.properties["partial_charge"] = charges[i]

        self.check_partial_charges()

    @property
    def n_atoms(self) -> int:
        """Number of atoms in the molecule."""
        return len(self)

    @property
    def n_bonds(self) -> int:
        """Number of bonds in the molecule."""
        return len(self.bonds)
    
    @property
    def center_of_mass(self) -> NDArray[np.float64]:
        """Calculate the center of mass of the molecule."""
        total_mass = 0.0
        weighted_positions = np.zeros(3)

        for site in self.sites:
            mass = site.atom.atomic_mass
            total_mass += mass
            weighted_positions += mass * np.array(site.coordinates)

        return weighted_positions / total_mass

    @property
    def geometrical_center(self) -> NDArray[np.float64]:
        """Calculate the geometric center of the molecule."""
        all_coords = np.array([site.coordinates for site in self.sites])
        return np.mean(all_coords, axis=0)

    @abstractmethod
    def get_distance(self, i: int, j: int) -> float:
        """Get the distance between i and j atoms in the molecule."""
        atom1 = self.sites[i]
        atom2 = self.sites[j]
        return float(
            np.linalg.norm(np.array(atom1.coordinates) - np.array(atom2.coordinates))
        )

    @abstractmethod
    def get_distances(self) -> np.ndarray:
        """Get the distance matrix for the molecule."""
        n_atoms = len(self)
        dist_matrix = np.zeros((n_atoms, n_atoms))
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                dist = self.get_distance(i, j)
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist
        return dist_matrix

    def get_closest_atom(self, atom_index: int) -> int:
        """Get the index of the closest atom to the given atom index."""

        distances = self.get_distances()[atom_index]

        # Get the index of the closest atom (excluding itself)
        closest_index = int(
            np.argmin(
                [
                    dist if idx != atom_index else np.inf
                    for idx, dist in enumerate(distances)
                ]
            )
        )
        return closest_index

    def set_bonds(
        self, bonds: Sequence[tuple[int, int]], bondTypes: Sequence[int]
    ) -> None:
        """Set the bonds for the molecule.

        Bond index should start at 1.
        """

        bond_dict = {
            "1": "single",
            "2": "double",
            "3": "triple",
            "4": "aromatic",
            "5": "ionic",
            "6": "hydrogen",
            "7": "metallic",
            "8": "coordinate",
        }

        for i, (index1, index2) in enumerate(bonds):
            bond_type_name = bond_dict.get(str(bondTypes[i]), "unknown")
            self.bonds.append(
                Bond(
                    atom_1=self.sites[index1 - 1],
                    atom_2=self.sites[index2 - 1],
                    bond_type=bondTypes[i],
                    bond_type_name=bond_type_name,
                    bond_distance=self.get_distance(index1 - 1, index2 - 1),
                )
            )

    def add_site_property(self, property_name: str, values: Sequence | list) -> None:
        """Add a property to a site. Note: This is the preferred method
        for adding magnetic moments, selective dynamics, and related
        site-specific properties to a structure/molecule object.

        Examples:
            structure.add_site_property("magmom", [1.0, 0.0])
            structure.add_site_property("selective_dynamics", [[True, True, True], [False, False, False]])

        Args:
            property_name (str): The name of the property to add.
            values (list): A sequence of values. Must be same length as
                number of sites.

        Raises:
            ValueError: if len(values) != number of sites.

        Returns:
            SiteCollection: self with site property added.
        """
        if len(values) != len(self):
            raise ValueError(
                f"{len(values)=} must equal sites in structure={len(self)}"
            )
        for site, val in zip(self, values, strict=True):
            site.properties[property_name] = val

    def get_site_property(self, property_name: str) -> list:
        """Get a site property.

        Args:
            property_name (str): The name of the property to get.

        Raises:
            KeyError: if property_name not found.

        Returns:
            list: List of property values for each site.
        """
        values = []
        for site in self:
            if property_name not in site.properties:
                raise KeyError(
                    f"Property {property_name} not found in site properties. Available properties: {list(site.properties.keys())}"
                )
            values.append(site.properties[property_name])
        return values

    def centralize(self, by_atom: str | None = None) -> None:
        """
        Centralize the molecule by subtracting the mean coordinates.
        If by_x is True, centralize by the mean position of the X atoms.

        Parameters
        ----------
        by_atom : str | None
            If provided, centralize by the mean position of the atoms of this type.
        """

        if by_atom:
            x_coords = [
                site.coordinates for site in self.sites if site.atom_type == by_atom
            ]
            mean_coords = np.mean(x_coords, axis=0)
        else:
            mean_coords = self.geometrical_center

        for site in self.sites:
            site.coordinates = site.coordinates - mean_coords

    def bond_atoms(self, index_1, index_2, bond_type: int = 1) -> None:
        """
        Create a bond between two atoms in the molecule.

        Parameters
        ----------
        index_1 : int
            Index of the first atom.
        index_2 : int
            Index of the second atom.
        bond_type : int, optional
            Type of the bond. Default is 1 (single bond).
        """

        bond_dict = {
            "1": "single",
            "2": "double",
            "3": "triple",
            "4": "aromatic",
            "5": "ionic",
            "6": "hydrogen",
            "7": "metallic",
            "8": "coordinate",
        }

        bond_type_name = bond_dict.get(str(bond_type), "unknown")

        self.bonds.append(
            Bond(
                atom_1=self.sites[index_1],
                atom_2=self.sites[index_2],
                bond_type=bond_type,
                bond_type_name=bond_type_name,
                bond_distance=self.get_distance(index_1, index_2),
            )
        )

    def replace_by_atom_group(
            self,
            atom_index: int,
            group: Molecule,
            type_to_replace: str = "R") -> None:
        """
        Repleace an atom in the molecule with a group of atoms from another molecule.

        The new molecule will be added such as the atom with the same atom type as the replaced occupies the same position.
        The bond between the replaced atom and the closest atom in the current molecule will be replaced by a bond between the
        closest atom and the new group.

        Parameters
        ----------
        atom_index : int
            Index of the atom to be replaced.
        group : Molecule
            Molecule object representing the group to be added.
        type_to_replace : str
            Atom type in the group molecule that will be aligned to the replaced atom. Default is "R".
        """

        closest_atom_index_molecule = self.get_closest_atom(atom_index)
        # Get the alignment vector in the current molecule
        alignment_vector = (
            self.sites[atom_index].coordinates
            - self.sites[closest_atom_index_molecule].coordinates
        )

        group = group.copy()

        group.centralize(by_atom=type_to_replace)

        # Get the atom of type_to_replace on the group
        r_atom_indices = [
            i for i, site in enumerate(group.sites) if site.atom_type == type_to_replace
        ]

        if not r_atom_indices:
            raise ValueError(
                f"The group molecule must contain an atom of type '{type_to_replace}' to align with the replaced atom."
            )

        if len(r_atom_indices) > 1:
            raise ValueError(
                f"The group molecule must contain only one atom of type '{type_to_replace}' to align with the replaced atom."
            )

        r_atom_index = r_atom_indices[0]

        closest_atom_index_group = group.get_closest_atom(r_atom_index)

        align_vector_group = (
            group.sites[closest_atom_index_group].coordinates
            - group.sites[r_atom_index].coordinates
        )

        Rot_m = rotation_matrix_from_vectors(align_vector_group, alignment_vector)

        # Rotate and translade the conector group to Q position in the strucutre
        group.positions = (
            np.dot(group.positions, np.transpose(Rot_m))
            + self.sites[atom_index].coordinates
        )

        # Add the group to the molecule
        for i, site in enumerate(group.sites):
            site.index = len(self.sites) + i

        self.sites.extend(group.sites)
        self.bonds.extend(group.bonds)

        # Bond the atoms together
        self.bond_atoms(
            closest_atom_index_molecule,
            closest_atom_index_group + len(self.sites) - len(group.sites),
            bond_type=1,
        )

        self.remove_atom(atom_index)

        self.remove_atom(r_atom_index + len(self.sites) - len(group.sites))

    def view(self) -> None:
        """
        Visualize the molecule using the ase viewrer.
        """

        from ase import Atoms
        from ase.visualize import view

        symbols = [site.atom_type for site in self.sites]
        positions = [site.coordinates for site in self.sites]

        # Replace 'R' atoms with 'X' for visualization
        for i, symbol in enumerate(symbols):
            if symbol in ["R", "Xe", "Q"]:
                symbols[i] = "X"

        ase_molecule = Atoms(symbols=symbols, positions=positions)

        # Set initial charges if charge property exists
        if "partial_charge" in self.sites[0].properties:
            charges = [site.properties["partial_charge"] for site in self.sites]
            ase_molecule.set_initial_charges(charges)

        view(ase_molecule)

    def save_mol(self, path, filename) -> None:
        """
        Saves the molecule to a .mol file (V2000 format).

        Parameters
        ----------
        path : str
            Path to save the .mol file.
        filename : str
            Name of the .mol file (without extension).
        """

        # Check if the filename has an extension
        if filename.endswith(".mol"):
            filename = filename[:-4]

        mol_txt = [f"{filename}\n\nCreated by pyCOFBuilder"]

        mol_txt += [
            f"{len(self.sites):>3}{len(self.bonds):>3}  0  0  0  0  0  0  0  0  0    0"
        ]

        for atom in self.sites:
            mol_txt += [
                "{:>10.4f}{:>10.4f}{:>10.4f} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0".format(
                    *atom.coordinates, atom.atom_type
                )
            ]

        for bond in self.bonds:
            mol_txt += [
                "{:>3}{:>3}{:>3}  0  0  0  0".format(
                    bond.atom_1.index + 1, bond.atom_2.index + 1, bond.bond_type
                )
            ]
        mol_txt += ["M  END"]

        with open(f"{path}\\{filename}.mol", "w") as f:
            f.write("\n".join(mol_txt))

    def save_xyz(self, path, filename) -> None:
        """
        Saves the molecule to a .xyz file on the file standard for poremake. 
        It includes partial charges and bond information.
        Parameters
        ----------
        path : str
            Path to save the .xyz file.
        filename : str
            Name of the .xyz file (without extension).
        """

        if filename.endswith(".xyz"):
            filename = filename[:-4]

        bond_dict = {1: "S", 2: "D", 3: "T", 4: "A"}

        X_index = [atom.index for atom in self.sites if atom.atom_type == "X"]

        xyz_txt = [f"{self.n_atoms}", "  ".join(map(str, X_index))]

        for atom in self.sites:
            xyz_txt += [
                "{:3s}    {:12.7f}  {:12.7f}  {:12.7f}  {:12.7f}".format(
                    atom.atom_type, *atom.coordinates, atom.properties.get("partial_charge", 0.0)
                )
            ]

        for bond in self.bonds:
            xyz_txt += [
                "{:3}  {:3}   {:3}".format(
                    bond.atom_1.index, bond.atom_2.index, bond_dict[bond.bond_type]
                )
            ]

        with open(os.path.join(path, filename + ".xyz"), "w") as f:
            f.write("\n".join(xyz_txt))

    def from_mol(self, file_name) -> None:
        """
        Load a molecule from a .mol file.
        """

        atomTypes, cartPos, partialCharges, bonds, bondTypes = read_mol_file(file_name)

        for i in range(len(atomTypes)):
            self.sites.append(
                AtomicSite(
                    atom_type=atomTypes[i],
                    label=atomTypes[i],
                    coordinates=np.array(cartPos[i]),
                    index=i,
                    properties={"partial_charge": partialCharges[i]},
                )
            )

        self.set_bonds(bonds, bondTypes)
