"""Measurement sets whose rows have ANTENNA1 > ANTENNA2 (e.g. SKA-Low).

A row stored as (p, q) with p > q holds V_pq = V_qp^H. calibrate turns every
such row back into V_qp (data and model conjugate-transposed) before solving,
and applysolutions corrects every row as S_a1 V S_a2^H with its own antennas.

Solutions for a measurement set with every row reversed are written as the
complex conjugate of the physical solutions, conj(J^-1). That is what
calibrate has always written for such data, so solution files made by
earlier versions stay valid (``test_regression.py`` checks this against a
reference build). A measurement set that mixes both orders is rejected.

Before this was fixed, reversed rows were solved as if they were V_qp and
applysolutions swapped their antennas to match. Parallel hands (XX/YY) came
out right for diagonal gains and an unpolarised sky, but the cross hands
(XY/YX) and leakage were wrong. The tests marked "fails on main" below show
this against the old code.
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
NANT = 36


def reverse_all_rows(ms_path: Path) -> None:
    """Store every cross-correlation as (higher antenna, lower antenna)"""
    ant1, ant2 = msgen.antennas(ms_path)
    msgen.reverse_baselines(ms_path, ant1 != ant2)


def polarised_model(model: np.ndarray) -> np.ndarray:
    """Turn an unpolarised model into one with Q, U and V.

    Uses fractional Stokes q=0.1, u=0.05, v=0.02, so the brightness matrix
    [[I+Q, U+iV], [U-iV, I-Q]] has non-zero cross hands.
    """
    q, u, v = 0.1, 0.05, 0.02
    xx = model[..., 0]
    out = np.empty_like(model)
    out[..., 0] = (1 + q) * xx
    out[..., 1] = (u + 1j * v) * xx
    out[..., 2] = (u - 1j * v) * xx
    out[..., 3] = (1 - q) * xx
    return out


def _calibrate_and_apply(bins, ms_path: Path, tmpdir, **kwargs) -> np.ndarray:
    Path(tmpdir).mkdir(parents=True, exist_ok=True)
    sol_path = Path(tmpdir) / "sols.bin"
    run(bins.calibrate, calibrate_args(ms_path, sol_path, **kwargs))
    run(bins.applysolutions, applysolutions_args(ms_path, sol_path))
    return AOSolutions.load(sol_path).bandpass[0]


def _corrected_error(ms_path: Path, correlations=slice(None)) -> np.ndarray:
    ant1, ant2 = msgen.antennas(ms_path)
    cross = ant1 != ant2
    corrected = msgen.get_column(ms_path, "CORRECTED_DATA")[cross][..., correlations]
    model = msgen.get_column(ms_path, "MODEL_DATA")[cross]
    return np.abs(corrected - model[..., correlations]) / np.max(np.abs(model))


def _assert_corrected(ms_path: Path, correlations=slice(None)) -> None:
    error = _corrected_error(ms_path, correlations)
    assert np.median(error) < MEDIAN_TOLERANCE
    assert np.max(error) < WORST_BASELINE_TOLERANCE


def _gain_error(solutions: np.ndarray, jones: np.ndarray, ms_path: Path) -> float:
    ant1, ant2 = msgen.antennas(ms_path)
    return float(np.max(msgen.gain_product_error(solutions, jones, ant1, ant2)))


@pytest.mark.parametrize(
    "kwargs", [{"extra": ["-diag"]}, {}], ids=["diag", "full_jones_solve"]
)
def test_bandpass_parallel_hands_and_convention(cal_ms_factory, bins, tmpdir, kwargs):
    """Diagonal gains, unpolarised sky: corrected data is right, and the
    solutions are conj(J^-1), the convention calibrate has always used for
    reversed data. Passes on main as well."""
    ms_path, truth = cal_ms_factory(diagonal=True)
    reverse_all_rows(ms_path)
    solutions = _calibrate_and_apply(bins, ms_path, tmpdir, **kwargs)

    _assert_corrected(ms_path)
    assert _gain_error(solutions, np.conj(truth.jones), ms_path) < GAIN_TOLERANCE
    assert _gain_error(solutions, truth.jones, ms_path) > 0.5


def test_polarised_sky_matches_normal_order(
    cal_ms_factory, template_model, bins, tmpdir
):
    """Diagonal gains, polarised sky: a reversed MS gives the same solutions
    (conjugated) and the same corrected data, all four correlations, as the
    same data in normal order.

    This compares with normal order rather than the truth because the solver
    only slowly finds the X-Y phase from a weakly polarised sky, in either
    order. Fails on main.
    """
    polarised = polarised_model(template_model)
    normal_ms, _ = cal_ms_factory(name="normal", diagonal=True, model_data=polarised)
    reversed_ms, _ = cal_ms_factory(
        name="reversed", diagonal=True, model_data=polarised
    )
    reverse_all_rows(reversed_ms)

    solutions, corrected = {}, {}
    for label, ms_path in (("normal", normal_ms), ("reversed", reversed_ms)):
        solutions[label] = _calibrate_and_apply(
            bins, ms_path, Path(tmpdir) / label, model=None, extra=["-diag"]
        )
        corrected[label] = msgen.get_column(ms_path, "CORRECTED_DATA")

    np.testing.assert_allclose(
        solutions["reversed"], np.conj(solutions["normal"]), rtol=1e-5, atol=1e-6
    )
    ant1, ant2 = msgen.antennas(normal_ms)
    cross = ant1 != ant2
    scale = np.max(np.abs(corrected["normal"][cross]))
    np.testing.assert_allclose(
        corrected["reversed"][cross],
        msgen.conjugate_transpose(corrected["normal"][cross]),
        atol=1e-5 * scale,
    )


def test_applysolutions_cross_hands_with_known_solutions(
    cal_ms_factory, template_model, bins, tmpdir
):
    """applysolutions alone: exact solutions in the reversed convention
    correct all four correlations of a polarised sky.

    Fails on main (the XY/YX correction).
    """
    ms_path, truth = cal_ms_factory(
        diagonal=True, model_data=polarised_model(template_model)
    )
    reverse_all_rows(ms_path)
    inverse = np.linalg.inv(truth.jones).reshape(NANT, -1, 4)
    sol_path = AOSolutions(
        path=Path(tmpdir) / "exact.bin",
        nsol=1,
        nant=NANT,
        nchan=inverse.shape[1],
        npol=4,
        bandpass=np.conj(inverse)[None],
    ).save(Path(tmpdir) / "exact.bin")

    run(bins.applysolutions, applysolutions_args(ms_path, sol_path))
    _assert_corrected(ms_path, correlations=[1, 2])
    _assert_corrected(ms_path)


@pytest.mark.parametrize("model", ["file", "model_data"])
def test_leakage_is_solved(cal_ms_factory, bins, tmpdir, model):
    """Full-Jones gains with leakage are solved and corrected.

    Fails on main, which solved the XY/YX-swapped equations.
    """
    ms_path, truth = cal_ms_factory()
    reverse_all_rows(ms_path)
    kwargs = {} if model == "file" else {"model": None}
    solutions = _calibrate_and_apply(bins, ms_path, tmpdir, **kwargs)

    assert _gain_error(solutions, np.conj(truth.jones), ms_path) < GAIN_TOLERANCE
    _assert_corrected(ms_path)


def test_noisy_reversed_matches_normal_order(cal_ms_factory, bins, tmpdir):
    """With noise and leakage, a reversed MS is solved as well as the same data
    in normal order: its solutions are exactly the conjugates."""
    noise = np.full(NANT, 1.0)
    normal_ms, _ = cal_ms_factory(name="normal", noise_sigma=noise)
    reversed_ms, _ = cal_ms_factory(name="reversed", noise_sigma=noise)
    reverse_all_rows(reversed_ms)

    solutions = {}
    for label, ms_path in (("normal", normal_ms), ("reversed", reversed_ms)):
        sol_path = Path(tmpdir) / f"{label}.bin"
        run(bins.calibrate, calibrate_args(ms_path, sol_path))
        solutions[label] = AOSolutions.load(sol_path).bandpass

    np.testing.assert_allclose(
        solutions["reversed"], np.conj(solutions["normal"]), rtol=1e-5, atol=1e-6
    )


@pytest.mark.parametrize("program", ["calibrate", "applysolutions"])
def test_mixed_row_order_is_rejected(cal_ms_factory, bins, tmpdir, ao_sols, program):
    """An MS with both row orders is an error in both programs.

    Fails on main, which gave wrong solutions instead.
    """
    ms_path, _ = cal_ms_factory()
    ant1, ant2 = msgen.antennas(ms_path)
    msgen.reverse_baselines(ms_path, (ant1 != ant2) & (np.arange(len(ant1)) % 2 == 1))

    if program == "calibrate":
        args = calibrate_args(ms_path, Path(tmpdir) / "sols.bin")
    else:
        args = applysolutions_args(ms_path, ao_sols)
    result = run(getattr(bins, program), args, check=False)

    assert result.returncode != 0
    assert "mixes row orders" in result.stderr


@pytest.mark.xfail(
    strict=True,
    reason="Solution files follow the row order of the MS they were made from "
    "(conjugated for reversed MSs), so they do not transfer to an MS with the "
    "opposite order. Unchanged from main.",
)
def test_solutions_transfer_between_row_orders(cal_ms_factory, bins, tmpdir):
    """Solutions from a normally ordered MS applied to a reversed copy"""
    normal_ms, _ = cal_ms_factory(name="normal")
    sol_path = Path(tmpdir) / "normal.bin"
    run(bins.calibrate, calibrate_args(normal_ms, sol_path))

    reversed_ms, _ = cal_ms_factory(name="reversed")
    reverse_all_rows(reversed_ms)
    run(bins.applysolutions, applysolutions_args(reversed_ms, sol_path))

    _assert_corrected(reversed_ms)
