from pathlib import Path
from pycofbuilder.molecule import Molecule

class BuildingBlock(Molecule):

    def __init__(
            self,
            database_path: Path,
            ) -> None:

        super().__init__()

        self.database_path: Path = database_path

        self.allowed_shapes: list[str] = [c.stem for c in (self.database_path / 'core').glob('*') if c.is_dir()]
        self.allowed_connectors: list[str] = [c.stem for c in (self.database_path / 'connGroups').glob('*')]
        self.allowed_funcGroups: list[str] = [c.stem for c in (self.database_path / 'funcGroups').glob('*')]

    def __repr__(self) -> str:
        return (
            "BuildingBlock(shape={}, core={}, connector={}, funcGroups={})".format(
                self.shape, self.core, self.connector, self.funcGroups
            )
        )
    
    def _check_name(self, name: str) -> None:
        if len(name.split('_')) < 3:
            raise ValueError(
                f"Name '{name}' is not in the correct format. " \
                 "Expected format: 'SHAPE_CORE_CONNECTOR_FUNCGROUP1_FUNCGROUP2_...'"
                 )
        
    def _check_shape(self) -> None:
        if self.shape not in self.allowed_shapes:
            raise ValueError(
                f"Shape '{self.shape}' is not allowed. Allowed shapes are: {self.allowed_shapes}"
                )
    
    def _check_core(self) -> None:
        pass

    def _check_connector(self) -> None:
        if self.connector not in self.allowed_connectors:
            raise ValueError(
                f"Connector '{self.connector}' is not allowed. Allowed connectors are: {self.allowed_connectors}"
                )
    
    def _check_funcGroups(self) -> None:
        # Get the list of functionalizable sites in the building block
        func_sites = list(set([s for s in self.atom_types if s.startswith("R")]))

        if len(self.funcGroups) > len(func_sites):
            raise ValueError(
                f"Number of functional groups ({len(self.funcGroups)}) exceeds the number of functionalizable sites ({len(func_sites)})."
                )
    
        for funcGroup in self.funcGroups:
            if funcGroup not in self.allowed_funcGroups:
                raise ValueError(
                    f"Functional group '{funcGroup}' is not allowed. Allowed functional groups are: {self.allowed_funcGroups}"
                    )

        self.funcGroups = self.funcGroups + ["H"] * (len(func_sites) - len(self.funcGroups))

    def add_conector(self) -> None:
        """
        Replace all the "Q" atoms in the building block with the connector molecule.
        The connector molecule is read from the database based on the connector name.
        """

        connector: Molecule = Molecule()
        connector.from_mol(self.database_path / 'connGroups' / f"{self.connector}.mol")

        while "Q" in self.atom_types:
            self.replace_by_atom_group(
                atom_index=self.atom_types.index("Q"),
                group=connector,
                type_to_replace="Q")

    def add_funcGroups(self) -> None:
        """
        Replace all the "RX" atoms in the building block with the functional group molecule.
        The functional group molecule is read from the database based on the functional group name.
        """
        for i, funcGroup_name in enumerate(self.funcGroups):
            funcGroup: Molecule = Molecule()
            funcGroup.from_mol(self.database_path / 'funcGroups' / f"{funcGroup_name}.mol")

            print(f"Adding functional group {funcGroup_name} as R{i+1}...")
            while f"R{i+1}" in self.atom_types:
                self.replace_by_atom_group(
                    atom_index=self.atom_types.index(f"R{i+1}"),
                    group=funcGroup,
                    type_to_replace="R")
                
        # Replace any remaining "RX" atoms with the first H functional group
        remaining_groups = list(set([s for s in self.atom_types if s.startswith("R")]))

        for group in remaining_groups:
            while group in self.atom_types:
                funcGroup: Molecule = Molecule()
                funcGroup.from_mol(self.database_path / 'funcGroups' / "H.mol")
                self.replace_by_atom_group(
                    atom_index=self.atom_types.index(group),
                    group=funcGroup,
                    type_to_replace="R")

    def from_name(self, name: str) -> None:
        # Extract the shape, core, connector, and functional groups from the name
        self._check_name(name)

        self.shape: str = name.split('_')[0]
        self._check_shape()

        self.core: str = name.split('_')[1]
        self._check_core()

        self.from_mol(self.database_path / 'core' / self.shape / f"{self.core}.mol")

        self.connector: str = name.split('_')[2]
        self._check_connector()
        self.add_conector()
        
        self.funcGroups: list[str] = name.split('_')[3:]
        self._check_funcGroups()
        self.add_funcGroups()

        # Rebuild the name based on the shape, core, connector, and functional groups
        self.name: str = "{shape}_{core}_{connector}_{funcGroups}".format(
            shape=self.shape,
            core=self.core,
            connector=self.connector,
            funcGroups='_'.join(self.funcGroups)
        )

    def from_file(self, file_path: Path) -> None:
        """
        Read the building block from a file.
        """
        self.name = file_path.stem

        # Check file extension
        if file_path.suffix == '.mol':
            self.from_mol(file_path)
        else:
            raise NotImplementedError(
                f"File format '{file_path.suffix}' is not supported. Supported formats are: .mol"
                )
