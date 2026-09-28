"""calibrate recovers known Jones matrices, and applysolutions undoes them"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tests import msgen
from tests.aosolutions import AOSolutions
from tests.runner import applysolutions_args, calibrate_args, run

# Median error to expect from noise free data. The solver stops once the mean
# squared change per iteration drops below 1e-6, so individual baselines can be
# off by up to about 0.5%.
NOISE_FREE_TOLERANCE = 1e-3
WORST_BASELINE_TOLERANCE = 1e-2


def assert_restores_model(output: np.ndarray, model: np.ndarray) -> None:
    """Corrected visibilities match the model to the solver's precision"""
    scale = np.max(np.abs(model))
    error = np.abs(output - model) / scale
    assert np.median(error) < NOISE_FREE_TOLERANCE
    assert np.max(error) < WORST_BASELINE_TOLERANCE


def _solve(bins, ms_path: Path, tmpdir, **kwargs) -> AOSolutions:
    sol_path = Path(tmpdir) / "sols.bin"
    run(bins.calibrate, calibrate_args(ms_path, sol_path, **kwargs))
    return AOSolutions.load(sol_path)


def test_flint_default_recovers_gains(cal_ms_factory, bins, tmpdir):
    """Flint's default calibrate command solves noise free full-Jones data"""
    ms_path, truth = cal_ms_factory()
    sols = _solve(bins, ms_path, tmpdir)
    ant1, ant2 = msgen.antennas(ms_path)

    error = msgen.gain_product_error(sols.bandpass[0], truth.jones, ant1, ant2)
    assert np.all(np.isfinite(error))
    assert np.max(error) < NOISE_FREE_TOLERANCE


def test_model_data_column_recovers_gains(cal_ms_factory, bins, tmpdir):
    """Without -m the MODEL_DATA column is used as the model"""
    ms_path, truth = cal_ms_factory()
    sols = _solve(bins, ms_path, tmpdir, model=None)
    ant1, ant2 = msgen.antennas(ms_path)

    error = msgen.gain_product_error(sols.bandpass[0], truth.jones, ant1, ant2)
    assert np.max(error) < NOISE_FREE_TOLERANCE


def test_diag_recovers_diagonal_gains(cal_ms_factory, bins, tmpdir):
    """-diag solves diagonal Jones matrices and leaves the leakage at zero"""
    ms_path, truth = cal_ms_factory(diagonal=True)
    sols = _solve(bins, ms_path, tmpdir, extra=["-diag"])
    ant1, ant2 = msgen.antennas(ms_path)

    error = msgen.gain_product_error(sols.bandpass[0], truth.jones, ant1, ant2)
    assert np.max(error) < NOISE_FREE_TOLERANCE
    assert np.all(sols.bandpass[0][..., [1, 2]] == 0)


def test_noisy_data_recovers_gains(cal_ms_factory, bins, tmpdir):
    """Uniform noise gives solutions close to the truth"""
    ms_path, truth = cal_ms_factory(noise_sigma=np.full(36, 1.0))
    sols = _solve(bins, ms_path, tmpdir)
    ant1, ant2 = msgen.antennas(ms_path)

    error = msgen.gain_product_error(sols.bandpass[0], truth.jones, ant1, ant2)
    assert np.median(error) < 0.02


def test_solution_intervals(cal_ms_factory, bins, tmpdir):
    """-t 1 gives one solution interval per timestep"""
    ms_path, truth = cal_ms_factory()
    sols = _solve(bins, ms_path, tmpdir, extra=["-t", "1"])
    ant1, ant2 = msgen.antennas(ms_path)

    assert sols.nsol == 3
    for interval in range(sols.nsol):
        error = msgen.gain_product_error(
            sols.bandpass[interval], truth.jones, ant1, ant2
        )
        assert np.max(error) < NOISE_FREE_TOLERANCE


@pytest.mark.parametrize(
    "variant", ["weight_spectrum", "sigma_spectrum", "empty_weight_spectrum"]
)
def test_unit_weight_columns_give_identical_solutions(
    cal_ms_factory, bins, tmpdir, variant
):
    """With unit weights the choice of weight column does not matter"""
    ms_path, _ = cal_ms_factory(name="reference", noise_sigma=np.full(36, 1.0))
    reference = Path(tmpdir) / "reference.bin"
    run(bins.calibrate, calibrate_args(ms_path, reference))

    variant_ms, _ = cal_ms_factory(name=variant, noise_sigma=np.full(36, 1.0))
    msgen.set_weight_variant(variant_ms, variant)
    variant_sols = Path(tmpdir) / f"{variant}.bin"
    run(bins.calibrate, calibrate_args(variant_ms, variant_sols))

    assert variant_sols.read_bytes() == reference.read_bytes()


@pytest.mark.parametrize("copy", [True, False])
def test_applysolutions_restores_model(cal_ms_factory, bins, tmpdir, copy):
    """Applying the solutions turns DATA back into MODEL_DATA"""
    ms_path, _ = cal_ms_factory()
    sol_path = Path(tmpdir) / "sols.bin"
    run(bins.calibrate, calibrate_args(ms_path, sol_path))
    run(bins.applysolutions, applysolutions_args(ms_path, sol_path, copy=copy))

    output = msgen.get_column(ms_path, "CORRECTED_DATA" if copy else "DATA")
    model = msgen.get_column(ms_path, "MODEL_DATA")
    ant1, ant2 = msgen.antennas(ms_path)
    cross = ant1 != ant2

    assert_restores_model(output[cross], model[cross])
    if copy:
        with pytest.raises(AssertionError):
            np.testing.assert_allclose(
                msgen.get_column(ms_path, "DATA")[cross], model[cross], rtol=1e-2
            )


@pytest.mark.parametrize("correlations", [[0, 3], [0]])
def test_applysolutions_to_fewer_polarisations(
    cal_ms_factory, bins, tmpdir, correlations
):
    """Full-Jones solutions apply to 2-pol (XX,YY) and 1-pol (XX) data.

    Only the diagonal of each solution is used, so diagonal gains are
    recovered exactly.
    """
    ms_path, _ = cal_ms_factory(diagonal=True)
    sol_path = Path(tmpdir) / "sols.bin"
    run(bins.calibrate, calibrate_args(ms_path, sol_path, extra=["-diag"]))

    reduced = msgen.reshape_ms(
        ms_path, Path(tmpdir) / "reduced.ms", correlations=correlations
    )
    run(bins.applysolutions, applysolutions_args(reduced, sol_path))

    output = msgen.get_column(reduced, "CORRECTED_DATA")
    model = msgen.get_column(reduced, "MODEL_DATA")
    ant1, ant2 = msgen.antennas(reduced)
    cross = ant1 != ant2

    assert_restores_model(output[cross], model[cross])
