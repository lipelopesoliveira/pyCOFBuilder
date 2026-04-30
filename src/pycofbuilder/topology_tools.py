import re
from pathlib import Path
import pymatgen.core as mg

def expand_topology_by_symmetry(parsed_data: dict) -> dict:
    """
    Expand the topology by applying the symmetry operations of the space group to the nodes and edge centers.

    This function uses pymatgen to create a structure from the parsed data, applying the symmetry operations of the specified space group,
    so it is important to ensure that the space group on the .cgd file is compatible with pymatgen's space group definitions.

    Parameters
    ----------
    parsed_data (dict): A dictionary containing the parsed information from the .cdg file

    Returns
    -------
    expanded_data (dict): A dictionary containing the expanded information, with the same format as the input but with the nodes and edge centers expanded by symmetry.
    """

    species = ['C'] * len(parsed_data["nodes"]) + ['H'] * len(parsed_data["edge_centers"])
    positions = [tuple([node["x"], node["y"], node["z"]]) for node in parsed_data["nodes"]] + parsed_data["edge_centers"]
    properties = {
        'connectivity': [node["connectivity"] for node in parsed_data["nodes"]] + [2 for _ in parsed_data["edge_centers"]],
        'id': [node["id"] for node in parsed_data["nodes"]] + [-i - 1 for i in range(len(parsed_data["edge_centers"]))],
        'kind': ['node'] * len(parsed_data["nodes"]) + ['edge_center'] * len(parsed_data["edge_centers"])
    }

    structure = mg.Structure.from_spacegroup(
            sg=parsed_data["spacegroup"],
            lattice=mg.Lattice.from_parameters(*parsed_data["cell"]),
            species=species,
            coords=positions,
            site_properties=properties,
        )
    
    expandded_data = {
        "name": parsed_data["name"],
        "spacegroup": parsed_data["spacegroup"],
        "cell": parsed_data["cell"],
        "nodes": [
            {
                'id': site.properties['id'],
                'connectivity': site.properties['connectivity'],
                'x': float(site.frac_coords[0]),
                'y': float(site.frac_coords[1]),
                'z': float(site.frac_coords[2])
            } for site in structure if site.properties['kind'] == 'node'],
        "edge_centers": [site.frac_coords.tolist() for site in structure if site.properties['kind'] == 'edge_center']
    }

    return expandded_data

def parse_cdg_file(file_path: str | Path, expand_by_symmetry: bool = True) -> dict:
    """
    Parse a .cdg file and extract the relevant information into a structured dictionary.

    Parameters
    ----------
    - file_path (str | Path): The path to the .cdg file.
    - expand_by_symmetry (bool): Whether to expand the structure by symmetry.

    Returns
    -------
    dict: A dictionary containing the parsed information.
    """
    # Read the file content
    with open(file_path, 'r') as file:
        content = file.read()

    parsed_data = parse_cdg_content(content)

    if expand_by_symmetry:
        parsed_data = expand_topology_by_symmetry(parsed_data)
        
    return parsed_data

def parse_val(val_str: str) -> float:
    """
    Helper function to convert strings (decimals or fractions) to floats.

    Parameters:
    val_str (str): The string to convert, which can be a decimal (e.g. "3.5") or a fraction (e.g. "1/2").

    Returns:
    float: The converted float value.
    """

    if '/' in val_str:
        numerator, denominator = val_str.split('/')
        return float(numerator) / float(denominator)
    return float(val_str)

def parse_cdg_content(content: str) -> dict:
    """
    Parse the content of a .cdg file and extract the relevant information into a structured dictionary.

    Parameters:
    content (str): The content of the .cdg file.

    Returns:
    dict: A dictionary containing the parsed information.
    """
    parsed_data = {
        "name": None,
        "spacegroup": None,
        "cell": [],
        "nodes": [],
        "edges": [],
        "edge_centers": []
    }

    # Extract NAME or ID (?:NAME|ID) means match "NAME" or "ID", but don't capture it as a group.
    name_match = re.search(r"(?:NAME|ID)\s+(\S+)", content)
    if name_match:
        parsed_data["name"] = name_match.group(1)

    # Extract space group information
    group_match = re.search(r"GROUP\s+(\S+)", content)
    if group_match:
        parsed_data["spacegroup"] = group_match.group(1)

    # Extract CELL parameters
    cell_pattern = r"CELL\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)"
    cell_match = re.search(cell_pattern, content)
    if cell_match:
        parsed_data["cell"] = [parse_val(val) for val in cell_match.groups()]

    # Extract NODE information
    node_pattern = r"NODE\s+(\d+)\s+(\d+)\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)"
    for match in re.finditer(node_pattern, content):
        parsed_data["nodes"].append({
            "id": int(match.group(1)),
            "connectivity": int(match.group(2)),
            "x": parse_val(match.group(3)),
            "y": parse_val(match.group(4)),
            "z": parse_val(match.group(5))
        })

    # Extract EDGE information
    edge_pattern = r"EDGE\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)"
    for match in re.finditer(edge_pattern, content):
        parsed_data["edges"].append({
            "start": (parse_val(match.group(1)), parse_val(match.group(2)), parse_val(match.group(3))),
            "end": (parse_val(match.group(4)), parse_val(match.group(5)), parse_val(match.group(6)))
        })

    # Extract EDGE_CENTER information
    edge_center_pattern = r"EDGE_CENTER\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)\s+([\d\.\-\/]+)"
    for match in re.finditer(edge_center_pattern, content):
        parsed_data["edge_centers"].append(
            (parse_val(match.group(1)), parse_val(match.group(2)), parse_val(match.group(3)))
        )

    # If there is no edge center information, we can calculate it from the edge information
    if not parsed_data["edge_centers"] and parsed_data["edges"]:
        for edge in parsed_data["edges"]:
            start = edge["start"]
            end = edge["end"]
            center = (
                (start[0] + end[0]) / 2,
                (start[1] + end[1]) / 2,
                (start[2] + end[2]) / 2
            )
            parsed_data["edge_centers"].append(center)

    return parsed_data
