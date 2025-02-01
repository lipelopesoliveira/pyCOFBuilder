import importlib
import pytest

# List of dependencies (only the package names)
dependencies = [
    "numpy",
    "scipy",
    "pymatgen",
    "simplejson",
    "pandas",
    "tqdm",
    "gemmi",
    "ase",
]


@pytest.mark.parametrize("package", dependencies)
def test_package_installed(package):
    """
    Tests if the package can be imported, which indicates it is installed.
    """
    try:
        importlib.import_module(package)
    except ImportError as e:
        pytest.fail(f"Package '{package}' is not installed. Details: {e}")
