from dataclasses import dataclass
from ase.cell import Cell
import numpy as np
from pathlib import Path

from pycofbuilder.topology_tools import parse_cdg_file

@dataclass
class Node:
    """
    Represents the node of a graph from a reticular material
    """
    id: int
    connectivity: int
    x: float
    y: float
    z: float
    
    def __repr__(self) -> str:
        return f"Node(id={self.id}, connectivity={self.connectivity}, x={self.x}, y={self.y}, z={self.z})"
    def __str__(self) -> str:
        return f"Node {self.id} ({self.connectivity}): ({self.x}, {self.y}, {self.z})"

@dataclass
class EdgeCenter:
    """
    Represents the edge center of a graph from a reticular material
    """
    x: float
    y: float
    z: float
    
    def __repr__(self) -> str:
        return f"EdgeCenter(x={self.x}, y={self.y}, z={self.z})"
    def __str__(self) -> str:
        return f"EdgeCenter: ({self.x}, {self.y}, {self.z})"

class Topology:
    """
    Represents the topology of a reticular material
    """
    def __init__(self, name: str = "") -> None:
        self._name: str = name
        self._spacegroup: str = "P1"
        self._cell: Cell = Cell(np.eye(3))
        self._nodes: list[Node] = []
        self._edge_centers: list[EdgeCenter] = []
        
        
    def __repr__(self) -> str:
        return "Topology(name={}, spacegroup={}, nodes={}, edge_centers={})".format(
            self.name,
            self.spacegroup,
            len(self.nodes),
            len(self.edge_centers)
            )
    
    def __str__(self) -> str:
        topo_str = f"Topology {self.name} \nSpacegroup: {self.spacegroup}\n"

        topo_str += "Cell: {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}\n".format(*self.cell.cellpar())
        for line in self.cell.array:
            topo_str += "  {:<10.3f} {:<10.3f} {:<10.3f}\n".format(*line)

        topo_str += f"Nodes:\n"
        for node in self.nodes:
            topo_str += f"  {node}\n"

        topo_str += f"Edge Centers:\n"
        for edge_center in self.edge_centers:
            topo_str += f"  {edge_center}\n"

        return topo_str
    
    @property
    def name(self) -> str:
        return self._name
    @name.setter
    def name(self, value: str) -> None:
        self._name = value

    @property
    def spacegroup(self) -> str:
        return self._spacegroup
    @spacegroup.setter
    def spacegroup(self, value: str) -> None:
        self._spacegroup = value

    @property
    def cell(self) -> Cell:
        return self._cell
    @cell.setter
    def cell(self, value: list[float]) -> None:
        if len(value) == 3:
            self._cell = Cell(value)
        elif len(value) == 6:
             self._cell = Cell.fromcellpar(value)
        else:
            raise ValueError("Cell must be a list of 3 or 6 floats")

    @property
    def nodes(self) -> list[Node]:
        return self._nodes
    @nodes.setter
    def nodes(self, value: list[Node]) -> None:
        self._nodes = value

    @property
    def edge_centers(self) -> list[EdgeCenter]:
        return self._edge_centers
    @edge_centers.setter
    def edge_centers(self, value: list[EdgeCenter]) -> None:
        self._edge_centers = value

    def read_cgd_file(self, cgd_file: str | Path, expand_by_symmetry: bool = True) -> None:
        """
        Parses a CGD file to populate the topology's nodes and edge centers
        """
        parsed_file = parse_cdg_file(cgd_file, expand_by_symmetry=expand_by_symmetry)

        self.name = parsed_file['name']
        self.spacegroup = parsed_file['spacegroup']
        self.cell = parsed_file['cell']

        self.nodes = [Node(**node) for node in parsed_file['nodes']]
        self.edge_centers = [EdgeCenter(*edge_center) for edge_center in parsed_file['edge_centers']]

    @classmethod
    def from_cgd_file(cls, cgd_file: str | Path, expand_by_symmetry: bool = True) -> 'Topology':
        """
        Parses a CGD file to populate the topology's nodes and edge centers
        """
        parsed_file = parse_cdg_file(cgd_file, expand_by_symmetry=expand_by_symmetry)

        topology = cls()
        topology.name = parsed_file['name']
        topology.spacegroup = parsed_file['spacegroup']
        topology.cell = parsed_file['cell']

        topology.nodes = [Node(**node) for node in parsed_file['nodes']]
        topology.edge_centers = [EdgeCenter(*edge_center) for edge_center in parsed_file['edge_centers']]

        return topology

    def to_cgd(self, output_file: str | Path) -> None:
        """
        Exports the topology to a CGD file
        """

        gcd_file = [
            'CRYSTAL',
            f"NAME {self.name}",
            f"GROUP {self.spacegroup}",
            "CELL {:<10.5f} {:<10.5f} {:<10.5f} {:<10.5f} {:<10.5f} {:<10.5f}".format(*self.cell.cellpar()),
        ]

        for node in self.nodes:
            gcd_file.append("NODE {:2} {:2} {:<10.5f} {:<10.5f} {:<10.5f}".format(node.id, node.connectivity, node.x, node.y, node.z))
        for edge_center in self.edge_centers:
            gcd_file.append("# EDGE_CENTER {:<10.5f} {:<10.5f} {:<10.5f}".format(edge_center.x, edge_center.y, edge_center.z))

        gcd_file.append("END")

        with open(output_file, 'w') as f:
            f.write("\n".join(gcd_file))