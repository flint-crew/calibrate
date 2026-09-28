"""Rows stored with ANTENNA1 > ANTENNA2.

``VisibilityArray::ValuePtr`` files every baseline under (smaller antenna,
larger antenna) without conjugating the data or model, so a reversed row's
V_qp = V_pq^H is solved as if it were V_pq. Flint PR #1 ("Correct flips")
made applysolutions swap the antennas of reversed rows in the same way.

For a measurement set where every row is reversed, with diagonal gains and
an unpolarised model (the bandpass case), this is self consistent: calibrate
solves for conj(J), and applysolutions undoes it, so the corrected data is
right. It is not right for leakage terms, or for a measurement set that mixes
both orders. These tests record what works and what does not.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tests import msgen
from tests.aosolutions import AOSolutions
from tests.runner import applysolutions_args, calibrate_args, run

# See test_calibrate.py for how these were chosen
GAIN_TOLERANCE = 3e-3
MEDIAN_TOLERANCE = 1e-3
WORST_BASELINE_TOLERANCE = 1e-2


def _reverse(ms_path: Path, every_row: bool) -> None:
    ant1, ant2 = msgen.antennas(ms_path)
    rows = ant1 != ant2
    if not every_row:
        rows &= np.arange(len(ant1)) % 2 == 1
    msgen.reverse_baselines(ms_path, rows)


def _calibrate_and_apply(bins, ms_path: Path, tmpdir, **kwargs) -> np.ndarray:
    sol_path = Path(tmpdir) / "sols.bin"
    run(bins.calibrate, calibrate_args(ms_path, sol_path, **kwargs))
    run(bins.applysolutions, applysolutions_args(ms_path, sol_path))
    return AOSolutions.load(sol_path).bandpass[0]


def _corrected_error(ms_path: Path) -> np.ndarray:
    ant1, ant2 = msgen.antennas(ms_path)
    cross = ant1 != ant2
    corrected = msgen.get_column(ms_path, "CORRECTED_DATA")[cross]
    model = msgen.get_column(ms_path, "MODEL_DATA")[cross]
    return np.abs(corrected - model) / np.max(np.abs(model))


def _gain_error(solutions: np.ndarray, jones: np.ndarray, ms_path: Path) -> float:
    ant1, ant2 = msgen.antennas(ms_path)
    return float(np.max(msgen.gain_product_error(solutions, jones, ant1, ant2)))


@pytest.mark.parametrize(
    "kwargs",
    [
        {"extra": ["-diag"]},
        {},
        {"model": None, "extra": ["-diag"]},
    ],
    ids=["diag", "full_jones_solve", "model_data_column"],
)
def test_all_reversed_diagonal_gains_are_corrected(
    cal_ms_factory, bins, tmpdir, kwargs
):
    """Every row reversed, diagonal gains: the corrected data is right.

    The solutions are the inverse of conj(J), not of J: their phases have the
    opposite sign to the same data stored in the usual order.
    """
    ms_path, truth = cal_ms_factory(diagonal=True)
    _reverse(ms_path, every_row=True)
    solutions = _calibrate_and_apply(bins, ms_path, tmpdir, **kwargs)

    error = _corrected_error(ms_path)
    assert np.median(error) < MEDIAN_TOLERANCE
    assert np.max(error) < WORST_BASELINE_TOLERANCE
    assert _gain_error(solutions, np.conj(truth.jones), ms_path) < GAIN_TOLERANCE
    assert _gain_error(solutions, truth.jones, ms_path) > 0.5


@pytest.mark.xfail(
    strict=True,
    reason="Reversed rows are solved without conjugating them, so solutions "
    "for an all-reversed MS are conj(J)^-1 rather than J^-1",
)
def test_all_reversed_solutions_follow_physical_convention(
    cal_ms_factory, bins, tmpdir
):
    """The same sky and gains give the same solutions whatever the row order"""
    ms_path, truth = cal_ms_factory(diagonal=True)
    _reverse(ms_path, every_row=True)
    solutions = _calibrate_and_apply(bins, ms_path, tmpdir, extra=["-diag"])

    assert _gain_error(solutions, truth.jones, ms_path) < GAIN_TOLERANCE


@pytest.mark.xfail(
    strict=True,
    reason="With every row reversed the XY and YX terms are swapped, so "
    "leakage can not be solved",
)
def test_all_reversed_leakage_is_corrected(cal_ms_factory, bins, tmpdir):
    """Every row reversed, full-Jones gains with leakage"""
    ms_path, _ = cal_ms_factory()
    _reverse(ms_path, every_row=True)
    _calibrate_and_apply(bins, ms_path, tmpdir)

    error = _corrected_error(ms_path)
    assert np.median(error) < MEDIAN_TOLERANCE
    assert np.max(error) < WORST_BASELINE_TOLERANCE


@pytest.mark.xfail(
    strict=True,
    reason="Reversed rows claim the opposite phase difference to normal rows, "
    "so an MS with both orders can not be solved",
)
def test_mixed_order_is_corrected(cal_ms_factory, bins, tmpdir):
    """Half the rows reversed, diagonal gains"""
    ms_path, _ = cal_ms_factory(diagonal=True)
    _reverse(ms_path, every_row=False)
    solutions = _calibrate_and_apply(bins, ms_path, tmpdir, extra=["-diag"])

    assert np.all(np.isfinite(solutions))
    error = _corrected_error(ms_path)
    assert np.median(error) < MEDIAN_TOLERANCE
    assert np.max(error) < WORST_BASELINE_TOLERANCE
