"""Reader and writer for the AO-style ``calibrate`` binary solutions file.

This is a copy of the ``AOSolutions`` structure and the
``load_aosolutions_file`` / ``save_aosolutions_file`` functions from
``flint.calibrate.aocalibrate`` (https://github.com/flint-crew/flint), with the
plotting and flagging helpers removed. Flint uses this code to read the files
written by ``calibrate``, so testing our output with it is the compatibility
check that matters.
"""

from __future__ import annotations

import struct
from pathlib import Path
from typing import NamedTuple

import numpy as np

HEADER_FORMAT = "8s6I2d"
HEADER_INTRO = b"MWAOCAL\0"
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)


class AOSolutions(NamedTuple):
    """Structure to load an AO-style solutions file"""

    path: Path
    """Path of the solutions file loaded"""
    nsol: int
    """Number of time solutions"""
    nant: int
    """Number of antenna in the solution file"""
    nchan: int
    """Number of channels in the solution file"""
    npol: int
    """Number of polarisations in the file"""
    bandpass: np.ndarray
    """Complex data representing the antenna Jones. Shape is (nsol, nant, nchan, npol)"""

    @classmethod
    def load(cls, path: Path) -> AOSolutions:
        """Load in an AO-style solution file. See `load_aosolutions_file`."""
        return load_aosolutions_file(solutions_path=path)

    def save(self, output_path: Path) -> Path:
        """Save the instance of AOSolution to a standard aosolution binary file

        Args:
            output_path (Path): Location to write the file to

        Returns:
            Path: Location the file was written to
        """
        return save_aosolutions_file(aosolutions=self, output_path=output_path)


def save_aosolutions_file(aosolutions: AOSolutions, output_path: Path) -> Path:
    """Save a AOSolutions file to the ao-standard binary format.

    Args:
        aosolutions (AOSolutions): Instance of the solutions to save
        output_path (Path): Output path to write the files to

    Returns:
        Path: Path the file was written to
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(str(output_path), "wb") as out_file:
        out_file.write(
            struct.pack(
                HEADER_FORMAT,
                HEADER_INTRO,
                0,  # File type, only 0 mode available
                0,  # Structure type, 0 model available only
                aosolutions.nsol,
                aosolutions.nant,
                aosolutions.nchan,
                aosolutions.npol,
                0.0,  # time start
                0.0,  # time end
            )
        )
        aosolutions.bandpass.astype("<c16").tofile(out_file)

    return output_path


def load_aosolutions_file(solutions_path: Path) -> AOSolutions:
    """Load in an AO-style solutions file

    Args:
        solutions_path (Path): The path of the solutions file to load

    Returns:
        AOSolutions: Structure container the deserialized solutions file
    """
    assert solutions_path.exists() and solutions_path.is_file(), (
        f"{solutions_path!s} either does not exist or is not a file. "
    )

    with open(solutions_path) as in_file:
        _junk = np.fromfile(in_file, dtype="<i4", count=2)

        header = np.fromfile(in_file, dtype="<i4", count=10)
        file_type = header[0]
        assert file_type == 0, f"Expected file_type of 0, found {file_type}"

        structure_type = header[1]
        assert structure_type == 0, (
            f"Expected structure_type of 0, found {structure_type}"
        )

        nsol, nant, nchan, npol = header[2:6]
        sol_shape = (nsol, nant, nchan, npol)

        bandpass = np.fromfile(in_file, dtype="<c16", count=np.prod(sol_shape)).reshape(
            sol_shape
        )

        return AOSolutions(
            path=solutions_path,
            nsol=int(nsol),
            nant=int(nant),
            nchan=int(nchan),
            npol=int(npol),
            bandpass=bandpass,
        )
