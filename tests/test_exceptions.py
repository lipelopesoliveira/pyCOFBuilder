from pycofbuilder.exceptions import (BondLenghError,
                                     BBConnectivityError,
                                     ConnectionGroupError,
                                     MissingXError)


def test_bond_length_error_message():
    # Create an instance of BondLenghError with specific parameters
    e = BondLenghError(atom_1="H", atom_2="O", dist=0.5, threshold=0.8)
    message = str(e)
    # Check that the message contains the expected parts
    assert "H" in message
    assert "O" in message
    assert "0.5" in message
    assert "0.8" in message
    # Optionally, check the complete expected message
    expected_message = "WARNING: Atoms H and O are closer than 0.5 A, 0.8"
    assert message == expected_message


def test_bb_connectivity_error_message():
    # Create an instance of BBConnectivityError with specific parameters
    e = BBConnectivityError(connectivity=3, found_connectivity=2)
    message = str(e)
    # Check that the message contains the expected parts
    assert "3" in message
    assert "2" in message
    expected_message = "ERROR: The building block connectivity should be 3 buy is 2"
    assert message == expected_message


def test_connection_group_error_message():
    # Create an instance of ConnectionGroupError with specific parameters
    e = ConnectionGroupError(conn_1="GroupA", conn_2="GroupB")
    message = str(e)
    # Check that the message contains the expected parts
    assert "GroupA" in message
    assert "GroupB" in message
    expected_message = "ERROR: The connection group GroupA not compatible with GroupB"
    assert message == expected_message


def test_missing_x_error_message():
    # Create an instance of MissingXError (without parameters)
    e = MissingXError()
    message = str(e)
    # Check that the message contains the expected warning
    expected_message = "No X points found in the structure!"
    assert message == expected_message
