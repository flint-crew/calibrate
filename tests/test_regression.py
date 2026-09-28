"""Compare this build with a reference build, output for output.

Set ``BASELINE_BIN_DIR`` to a build of the reference version (normally
``main``) to run these tests. Every existing command line option must give
byte-identical solutions and identical corrected visibilities, which is how
we make sure a change can not break the pipeline.

The measurement sets here have unit weights, so the weighting fixes are not
expected to change anything. Autocorrelation rows are left out of the
applysolutions comparison because the reference build writes stale data into
them (see ``test_known_issues.py``).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from tests import msgen
from tests.runner import (
    Binaries,
    addmodel_args,
    applysolutions_args,
    calibrate_args,
    run,
)

CALIBRATE_CASES = {
    "flint_default": {},
    "model_data_column": {"model": None},
    "no_minuv": {"minuv": None},
    "default_iterations": {"i": None},
    "intervals": {"extra": ["-t", "1"]},
    "scan_selection": {"extra": ["-startscan", "1", "-endscan", "2"]},
    "scan_selection_intervals": {"extra": ["-startscan", "1", "-t", "1"]},
    "refmode_0": {"extra": ["-refmode", "0"]},
    "refmode_1": {"extra": ["-refmode", "1"]},
    "refmode_2": {"extra": ["-refmode", "2"]},
    "diag": {"extra": ["-diag"]},
    "scalar": {"extra": ["-scalar"]},
    "rotation": {"extra": ["-rotation"]},
    "maxuv": {"extra": ["-maxuv", "3000"]},
    "accuracy": {"extra": ["-a", "1e-3", "1e-5"]},
    "few_iterations": {"i": 5},
    "threads": {"extra": ["-j", "1"]},
    "quiet": {"extra": ["-quiet"]},
}


@pytest.fixture(scope="module")
def regression_ms(tmp_path_factory, template_ms, template_model) -> Path:
    """A noisy, unit weight MS with one scan per timestep, shared read-only"""
    out_dir = Path(tmp_path_factory.mktemp("regression"))
    ms_path = out_dir / "regression.ms"
    msgen.make_calibration_ms(
        template_ms=template_ms,
        output_ms=ms_path,
        model_data=template_model,
        noise_sigma=np.full(36, 1.0),
    )
    msgen.set_scans_per_timestep(ms_path)
    return ms_path


def _fresh_copy(ms_path: Path, directory: Path, name: str) -> Path:
    out = directory / name
    shutil.copytree(ms_path, out)
    return out


def _calibrate(binaries: Binaries, ms_path: Path, sol_path: Path, **kwargs) -> Path:
    run(binaries.calibrate, calibrate_args(ms_path, sol_path, **kwargs))
    return sol_path


def _cross_rows(ms_path: Path) -> np.ndarray:
    ant1, ant2 = msgen.antennas(ms_path)
    return ant1 != ant2


@pytest.mark.parametrize("case", CALIBRATE_CASES, ids=list(CALIBRATE_CASES))
def test_calibrate_solutions_identical(
    case, regression_ms, bins, baseline_bins, tmpdir
):
    """Solutions are byte-identical to the reference build"""
    kwargs = CALIBRATE_CASES[case]
    tmpdir = Path(tmpdir)
    new = _calibrate(bins, regression_ms, tmpdir / "new.bin", **kwargs)
    old = _calibrate(baseline_bins, regression_ms, tmpdir / "old.bin", **kwargs)

    assert new.read_bytes() == old.read_bytes()


def test_calibrate_plot_files_identical(regression_ms, bins, baseline_bins, tmpdir):
    """The -p phase and gain plot files are identical to the reference build"""
    tmpdir = Path(tmpdir)
    outputs = {}
    for label, binaries in (("new", bins), ("old", baseline_bins)):
        phases, gains = tmpdir / f"{label}_phases.txt", tmpdir / f"{label}_gains.txt"
        _calibrate(
            binaries,
            regression_ms,
            tmpdir / f"{label}.bin",
            extra=["-p", str(phases), str(gains)],
        )
        outputs[label] = (phases.read_bytes(), gains.read_bytes())

    assert outputs["new"] == outputs["old"]


APPLY_CASES = {
    "copy": {"copy": True},
    "nocopy": {"copy": False},
    "scan_selection": {"extra": ["-startscan", "1", "-endscan", "2"]},
    "presets": {"extra": ["-s", "1", "0", "0", "1"]},
}


@pytest.mark.parametrize("case", APPLY_CASES, ids=list(APPLY_CASES))
@pytest.mark.parametrize("intervals", [False, True], ids=["one_interval", "t1"])
def test_applysolutions_identical(
    case, intervals, regression_ms, bins, baseline_bins, tmpdir
):
    """Corrected cross-correlations are identical to the reference build.

    The solutions come from the reference calibrate, so this also checks that
    the new applysolutions reads solutions written by the old calibrate.
    """
    if intervals and case == "scan_selection":
        pytest.skip("Scan selection with several intervals is a known issue")
    kwargs = APPLY_CASES[case]
    tmpdir = Path(tmpdir)
    sol_path = _calibrate(
        baseline_bins,
        regression_ms,
        tmpdir / "sols.bin",
        extra=["-t", "1"] if intervals else [],
    )

    column = "CORRECTED_DATA" if kwargs.get("copy", True) else "DATA"
    outputs = {}
    for label, binaries in (("new", bins), ("old", baseline_bins)):
        ms_path = _fresh_copy(regression_ms, tmpdir, f"{label}.ms")
        run(binaries.applysolutions, applysolutions_args(ms_path, sol_path, **kwargs))
        outputs[label] = msgen.get_column(ms_path, column)

    cross = _cross_rows(regression_ms)
    if case == "scan_selection":
        scans = msgen.get_column(regression_ms, "SCAN_NUMBER")
        cross &= scans >= 1
    np.testing.assert_array_equal(outputs["new"][cross], outputs["old"][cross])


def test_new_solutions_apply_with_old_applysolutions(
    regression_ms, bins, baseline_bins, tmpdir
):
    """Solutions written by the new calibrate work with the old applysolutions"""
    tmpdir = Path(tmpdir)
    outputs = {}
    for label, binaries in (("new", bins), ("old", baseline_bins)):
        sol_path = _calibrate(binaries, regression_ms, tmpdir / f"{label}.bin")
        ms_path = _fresh_copy(regression_ms, tmpdir, f"{label}.ms")
        run(baseline_bins.applysolutions, applysolutions_args(ms_path, sol_path))
        outputs[label] = msgen.get_column(ms_path, "CORRECTED_DATA")

    cross = _cross_rows(regression_ms)
    np.testing.assert_array_equal(outputs["new"][cross], outputs["old"][cross])


@pytest.mark.parametrize("correlations", [[0, 3], [0]], ids=["2pol", "1pol"])
def test_applysolutions_fewer_polarisations_identical(
    correlations, regression_ms, bins, baseline_bins, tmpdir
):
    """Applying 4-pol solutions to 2-pol and 1-pol data is unchanged"""
    tmpdir = Path(tmpdir)
    sol_path = _calibrate(baseline_bins, regression_ms, tmpdir / "sols.bin")
    reduced = msgen.reshape_ms(
        regression_ms, tmpdir / "reduced.ms", correlations=correlations
    )

    outputs = {}
    for label, binaries in (("new", bins), ("old", baseline_bins)):
        ms_path = _fresh_copy(reduced, tmpdir, f"{label}.ms")
        run(binaries.applysolutions, applysolutions_args(ms_path, sol_path))
        outputs[label] = msgen.get_column(ms_path, "CORRECTED_DATA")

    cross = _cross_rows(reduced)
    np.testing.assert_array_equal(outputs["new"][cross], outputs["old"][cross])


@pytest.mark.parametrize("mode", ["a", "s", "c", "z"])
@pytest.mark.parametrize("datacolumn", ["DATA", "NEW_MODEL"])
def test_addmodel_identical(
    mode, datacolumn, regression_ms, bins, baseline_bins, tmpdir
):
    """addmodel output is identical to the reference build.

    Autocorrelations of a newly created column are left out, because the
    reference build never writes them (see ``test_known_issues.py``).
    """
    if datacolumn == "NEW_MODEL" and mode != "c":
        pytest.skip("A new column is only used with copy mode by Flint")
    tmpdir = Path(tmpdir)
    outputs = {}
    for label, binaries in (("new", bins), ("old", baseline_bins)):
        ms_path = _fresh_copy(regression_ms, tmpdir, f"{label}.ms")
        run(binaries.addmodel, addmodel_args(ms_path, datacolumn=datacolumn, mode=mode))
        outputs[label] = msgen.get_column(ms_path, datacolumn)

    rows = _cross_rows(regression_ms) if datacolumn == "NEW_MODEL" else slice(None)
    np.testing.assert_array_equal(outputs["new"][rows], outputs["old"][rows])
