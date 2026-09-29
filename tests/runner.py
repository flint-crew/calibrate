"""Locate and run the calibrate programs.

The command helpers build the same command lines that Flint builds in
``flint.calibrate.aocalibrate`` and ``flint.predict.addmodel``, so these tests
exercise ``calibrate`` the way the pipeline uses it.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from collections.abc import Iterable
from pathlib import Path
from typing import NamedTuple

from tests.msgen import MODEL_1934

PROGRAMS = ("calibrate", "applysolutions", "addmodel")


class Binaries(NamedTuple):
    """Paths to one build of the three programs"""

    calibrate: Path
    """The calibrate program"""
    applysolutions: Path
    """The applysolutions program"""
    addmodel: Path
    """The addmodel program"""


def _is_exe(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)


def find_binaries(search: Iterable[Path], use_path: bool = True) -> Binaries | None:
    """Find the programs in the first directory that holds all of them.

    Args:
        search (Iterable[Path]): Directories to look in, in order
        use_path (bool, optional): Fall back to ``$PATH``. Defaults to True.

    Returns:
        Binaries | None: The programs, or None when they can not be found
    """
    for directory in search:
        paths = [Path(directory) / p for p in PROGRAMS]
        if all(_is_exe(p) for p in paths):
            return Binaries(*paths)

    if use_path:
        found = [shutil.which(p) for p in PROGRAMS]
        if all(found):
            return Binaries(*[Path(p) for p in found])

    return None


def run(
    program: Path, args: list[str], check: bool = True
) -> subprocess.CompletedProcess:
    """Run one of the programs, capturing its output.

    Args:
        program (Path): Program to run
        args (list[str]): Its arguments
        check (bool, optional): Fail if the exit code is non-zero. Defaults to True.

    Returns:
        subprocess.CompletedProcess: The finished process
    """
    result = subprocess.run(
        [str(program), *args], capture_output=True, text=True, check=False
    )
    if check and result.returncode != 0:
        raise AssertionError(
            f"{program.name} {' '.join(args)} exited with {result.returncode}\n"
            f"stdout:\n{result.stdout[-4000:]}\nstderr:\n{result.stderr[-4000:]}"
        )
    return result


def calibrate_args(
    ms_path: Path,
    solutions_path: Path,
    model: Path | None = MODEL_1934,
    datacolumn: str = "DATA",
    minuv: float | None = 600,
    i: int | None = 100,
    extra: Iterable[str] = (),
) -> list[str]:
    """Arguments for ``calibrate`` as Flint's ``create_calibrate_cmd`` builds them.

    Flint's default is
    ``calibrate -datacolumn DATA -m <model> -minuv 600 -i 100 <ms> <solutions>``.

    Args:
        ms_path (Path): Measurement set to calibrate
        solutions_path (Path): Solutions file to write
        model (Path | None, optional): Sky model. None calibrates against MODEL_DATA. Defaults to MODEL_1934.
        datacolumn (str, optional): Column to calibrate. Defaults to "DATA".
        minuv (float | None, optional): Minimum baseline length in metres. Defaults to 600.
        i (int | None, optional): Maximum number of iterations. Defaults to 100.
        extra (Iterable[str], optional): Additional options, placed before the positional arguments. Defaults to ().

    Returns:
        list[str]: The arguments
    """
    args = ["-datacolumn", datacolumn]
    if model is not None:
        args += ["-m", str(model)]
    if minuv is not None:
        args += ["-minuv", str(minuv)]
    if i is not None:
        args += ["-i", str(i)]
    args += list(extra)
    return [*args, str(ms_path), str(solutions_path)]


def applysolutions_args(
    ms_path: Path,
    solutions_path: Path,
    datacolumn: str = "DATA",
    copy: bool = True,
    extra: Iterable[str] = (),
) -> list[str]:
    """Arguments for ``applysolutions`` as Flint's ``create_apply_solutions_cmd`` builds them.

    Args:
        ms_path (Path): Measurement set to correct
        solutions_path (Path): Solutions file to apply
        datacolumn (str, optional): Column to correct. Defaults to "DATA".
        copy (bool, optional): Write CORRECTED_DATA (``-copy``) instead of overwriting (``-nocopy``). Defaults to True.
        extra (Iterable[str], optional): Additional options. Defaults to ().

    Returns:
        list[str]: The arguments
    """
    return [
        "-datacolumn",
        datacolumn,
        "-copy" if copy else "-nocopy",
        *extra,
        str(ms_path),
        str(solutions_path),
    ]


def addmodel_args(
    ms_path: Path, model: Path = MODEL_1934, datacolumn: str = "DATA", mode: str = "a"
) -> list[str]:
    """Arguments for ``addmodel`` as Flint's ``add_model_options_to_command`` builds them.

    Args:
        ms_path (Path): Measurement set to modify
        model (Path, optional): Sky model. Defaults to MODEL_1934.
        datacolumn (str, optional): Column to modify. Defaults to "DATA".
        mode (str, optional): a(dd), s(ubtract), c(opy) or z(ero). Defaults to "a".

    Returns:
        list[str]: The arguments
    """
    return ["-datacolumn", datacolumn, "-m", mode, str(model), str(ms_path)]
