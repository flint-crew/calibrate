"""Rows stored with ANTENNA1 > ANTENNA2.

A row stored as (p, q) with p > q holds V_pq = V_qp^H. calibrate turns such
rows into V_qp (data and model conjugate-transposed) before solving, and
applysolutions corrects every row as S_a1 V S_a2^H with its own antennas. So
the row order of a measurement set does not change the solutions.

Before this was fixed, reversed rows were solved as if they were V_qp, and
applysolutions (Flint PR #1) swapped their antennas to match. That was only
self consistent for an MS with every row reversed, diagonal gains and an
unpolarised model, and even then the solutions were conj(J)^-1. Leakage and
MSs with both orders were solved wrongly. Each test here fails on that code.
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


@pytest.mark.parametrize("every_row", [True, False], ids=["all_reversed", "mixed"])
@pytest.mark.parametrize(
    "diagonal, kwargs",
    [
        (True, {"extra": ["-diag"]}),
        (True, {}),
        (False, {}),
        (False, {"model": None}),
    ],
    ids=["diag", "diag_gains_full_solve", "leakage", "leakage_model_data_column"],
)
def test_reversed_rows_are_calibrated(
    cal_ms_factory, bins, tmpdir, every_row, diagonal, kwargs
):
    """Reversed rows recover the true gains, and are corrected, like normal rows"""
    ms_path, truth = cal_ms_factory(diagonal=diagonal)
    _reverse(ms_path, every_row=every_row)
    solutions = _calibrate_and_apply(bins, ms_path, tmpdir, **kwargs)

    assert np.all(np.isfinite(solutions))
    assert _gain_error(solutions, truth.jones, ms_path) < GAIN_TOLERANCE
    error = _corrected_error(ms_path)
    assert np.median(error) < MEDIAN_TOLERANCE
    assert np.max(error) < WORST_BASELINE_TOLERANCE


def test_row_order_does_not_change_solutions(cal_ms_factory, bins, tmpdir):
    """The same data in either row order gives the same solutions"""
    normal_ms, _ = cal_ms_factory(name="normal")
    reversed_ms, _ = cal_ms_factory(name="reversed")
    _reverse(reversed_ms, every_row=True)

    solutions = {}
    for label, ms_path in (("normal", normal_ms), ("reversed", reversed_ms)):
        sol_path = Path(tmpdir) / f"{label}.bin"
        run(bins.calibrate, calibrate_args(ms_path, sol_path))
        solutions[label] = AOSolutions.load(sol_path).bandpass

    np.testing.assert_allclose(
        solutions["reversed"], solutions["normal"], rtol=1e-5, atol=1e-6
    )


@pytest.mark.parametrize("every_row", [True, False], ids=["all_reversed", "mixed"])
def test_solutions_transfer_between_row_orders(cal_ms_factory, bins, tmpdir, every_row):
    """Solutions from a normally ordered MS correct a reversed copy of it"""
    normal_ms, _ = cal_ms_factory(name="normal")
    sol_path = Path(tmpdir) / "normal.bin"
    run(bins.calibrate, calibrate_args(normal_ms, sol_path))

    reversed_ms, _ = cal_ms_factory(name="reversed")
    _reverse(reversed_ms, every_row=every_row)
    run(bins.applysolutions, applysolutions_args(reversed_ms, sol_path))

    error = _corrected_error(reversed_ms)
    assert np.median(error) < MEDIAN_TOLERANCE
    assert np.max(error) < WORST_BASELINE_TOLERANCE
