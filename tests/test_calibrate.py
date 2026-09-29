"""calibrate recovers known Jones matrices, and applysolutions undoes them"""

from __future__ import annotations

import re
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
# Largest per channel gain_product_error to expect from noise free data
# (about 1.5e-3); solutions in the wrong convention give errors of order 1
GAIN_TOLERANCE = 3e-3


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
    assert np.max(error) < GAIN_TOLERANCE


def test_model_data_column_recovers_gains(cal_ms_factory, bins, tmpdir):
    """Without -m the MODEL_DATA column is used as the model"""
    ms_path, truth = cal_ms_factory()
    sols = _solve(bins, ms_path, tmpdir, model=None)
    ant1, ant2 = msgen.antennas(ms_path)

    error = msgen.gain_product_error(sols.bandpass[0], truth.jones, ant1, ant2)
    assert np.max(error) < GAIN_TOLERANCE


def test_diag_recovers_diagonal_gains(cal_ms_factory, bins, tmpdir):
    """-diag solves diagonal Jones matrices and leaves the leakage at zero"""
    ms_path, truth = cal_ms_factory(diagonal=True)
    sols = _solve(bins, ms_path, tmpdir, extra=["-diag"])
    ant1, ant2 = msgen.antennas(ms_path)

    error = msgen.gain_product_error(sols.bandpass[0], truth.jones, ant1, ant2)
    assert np.max(error) < GAIN_TOLERANCE
    assert np.all(sols.bandpass[0][..., [1, 2]] == 0)


def test_noisy_data_recovers_gains(cal_ms_factory, bins, tmpdir):
    """Uniform noise gives solutions close to the truth"""
    ms_path, truth = cal_ms_factory(noise_sigma=np.full(36, 1.0))
    sols = _solve(bins, ms_path, tmpdir)
    ant1, ant2 = msgen.antennas(ms_path)

    error = msgen.gain_product_error(sols.bandpass[0], truth.jones, ant1, ant2)
    assert np.median(error) < 0.03


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
        assert np.max(error) < GAIN_TOLERANCE


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


# Weights


def heterogeneous_noise() -> np.ndarray:
    """Four very noisy antennas among otherwise clean ones"""
    sigma = np.full(36, 1.0)
    sigma[:4] = 30.0
    return sigma


def clean_antenna_error(ms_path: Path, sol_path: Path, truth: msgen.Truth) -> float:
    """Median gain error over baselines between the clean antennas"""
    sols = AOSolutions.load(sol_path)
    ant1, ant2 = msgen.antennas(ms_path)
    error = msgen.gain_product_error(
        sols.bandpass[0], truth.jones, ant1, ant2, exclude=set(range(4))
    )
    return float(np.median(error))


def test_noisy_antennas_are_downweighted(cal_ms_factory, bins, tmpdir):
    """With WEIGHT = 1/SIGMA^2 the noisy antennas barely affect the others.

    With uniform noise of sigma=2 the median error is about 0.03, so correctly
    weighted data with sigma=1 on the clean antennas must do at least as well.
    Using SIGMA itself as the weight gave an error of about 0.77.
    """
    ms_path, truth = cal_ms_factory(noise_sigma=heterogeneous_noise())
    sol_path = Path(tmpdir) / "sols.bin"
    run(bins.calibrate, calibrate_args(ms_path, sol_path))

    assert clean_antenna_error(ms_path, sol_path, truth) < 0.05


@pytest.mark.parametrize(
    "variant", ["weight_spectrum", "sigma_spectrum", "empty_weight_spectrum"]
)
def test_weight_columns_agree(cal_ms_factory, bins, tmpdir, variant):
    """Every weight column layout describing the same noise gives the same answer"""
    ms_path, _ = cal_ms_factory(name="reference", noise_sigma=heterogeneous_noise())
    reference = Path(tmpdir) / "reference.bin"
    run(bins.calibrate, calibrate_args(ms_path, reference))

    variant_ms, _ = cal_ms_factory(name=variant, noise_sigma=heterogeneous_noise())
    msgen.set_weight_variant(variant_ms, variant)
    variant_sols = Path(tmpdir) / f"{variant}.bin"
    result = run(bins.calibrate, calibrate_args(variant_ms, variant_sols))

    expected_column = {
        "weight_spectrum": "WEIGHT_SPECTRUM",
        "sigma_spectrum": "WEIGHT",
        "empty_weight_spectrum": "WEIGHT",
    }[variant]
    assert f"Using {expected_column} column for weights" in result.stdout
    np.testing.assert_allclose(
        AOSolutions.load(variant_sols).bandpass,
        AOSolutions.load(reference).bandpass,
        rtol=1e-5,
        atol=1e-6,
    )


def test_weightcolumn_selects_column(cal_ms_factory, bins, tmpdir):
    """-weightcolumn overrides the automatic choice of weight column.

    WEIGHT_SPECTRUM is made uniform, so it ignores the noisy antennas, while
    WEIGHT and SIGMA describe the noise. WEIGHT_SPECTRUM is chosen by default;
    forcing WEIGHT or SIGMA gives the better, noise weighted solutions.
    """
    ms_path, truth = cal_ms_factory(noise_sigma=heterogeneous_noise())
    flags = msgen.get_column(ms_path, "FLAG")
    msgen.add_array_column(
        ms_path, "WEIGHT_SPECTRUM", np.ones(flags.shape, dtype=np.float32)
    )

    errors = {}
    for column in (None, "WEIGHT", "SIGMA", "WEIGHT_SPECTRUM"):
        sol_path = Path(tmpdir) / f"{column}.bin"
        extra = ["-weightcolumn", column] if column else []
        result = run(bins.calibrate, calibrate_args(ms_path, sol_path, extra=extra))
        assert f"Using {column or 'WEIGHT_SPECTRUM'} column" in result.stdout
        errors[column] = clean_antenna_error(ms_path, sol_path, truth)

    assert errors["WEIGHT"] < 0.05
    assert errors["SIGMA"] == pytest.approx(errors["WEIGHT"], rel=1e-3)
    assert errors[None] == errors["WEIGHT_SPECTRUM"]
    assert errors[None] > 2 * errors["WEIGHT"]


@pytest.mark.parametrize("column", ["NOT_A_COLUMN", "SIGMA_SPECTRUM"])
def test_weightcolumn_rejects_bad_column(cal_ms_factory, bins, tmpdir, column):
    """An unsupported or missing weight column is an error"""
    ms_path, _ = cal_ms_factory()
    result = run(
        bins.calibrate,
        calibrate_args(
            ms_path, Path(tmpdir) / "sols.bin", extra=["-weightcolumn", column]
        ),
        check=False,
    )
    assert result.returncode != 0
    assert column in result.stderr


def test_retry_reproduces_first_attempt(cal_ms_factory, bins, tmpdir):
    """A retry starts from the same state as the first attempt.

    calibrate retries a channel that uses all its iterations, starting again
    from unity. Nothing else changes between the attempts, so the retry must
    reach exactly the same precision. The weights used to be applied a second
    time on the retry.
    """
    ms_path, _ = cal_ms_factory(noise_sigma=heterogeneous_noise())
    result = run(
        bins.calibrate,
        # One thread, so the log lines of different channels do not interleave
        calibrate_args(ms_path, Path(tmpdir) / "sols.bin", i=5, extra=["-j", "1"]),
    )

    first = dict(
        re.findall(r"Recalculating channel (\d+) \(accuracy=(\S+)\)\.", result.stdout)
    )
    final = dict(
        re.findall(
            r"finished calibrating channel (\d+) in \d+ iterations, precision=(\S+)\.$",
            result.stdout,
            flags=re.MULTILINE,
        )
    )
    assert len(first) > 0, "No channel was retried, the test needs fewer iterations"
    mismatched = {ch: (acc, final[ch]) for ch, acc in first.items() if final[ch] != acc}
    assert mismatched == {}


def test_single_channel_ms_gives_finite_solutions(cal_ms_factory, bins, tmpdir):
    """A one channel MS calibrated against a model file gives usable solutions"""
    ms_path, _ = cal_ms_factory()
    single = msgen.reshape_ms(
        ms_path, Path(tmpdir) / "single.ms", channels=slice(100, 101)
    )
    sol_path = Path(tmpdir) / "sols.bin"
    run(bins.calibrate, calibrate_args(single, sol_path))

    sols = AOSolutions.load(sol_path)
    assert sols.nchan == 1
    assert np.all(np.isfinite(sols.bandpass))


def test_autocorrelations_are_corrected(cal_ms_factory, bins, tmpdir):
    """Autocorrelation rows get the same S_a V S_a^H correction as other rows"""
    ms_path, _ = cal_ms_factory()
    sol_path = Path(tmpdir) / "sols.bin"
    run(bins.calibrate, calibrate_args(ms_path, sol_path))
    run(bins.applysolutions, applysolutions_args(ms_path, sol_path))

    sols = AOSolutions.load(sol_path)
    ant1, ant2 = msgen.antennas(ms_path)
    auto = ant1 == ant2
    data = msgen.get_column(ms_path, "DATA")[auto]
    corrected = msgen.get_column(ms_path, "CORRECTED_DATA")[auto]

    expected = msgen.apply_solutions(data, sols.bandpass[0], ant1[auto], ant2[auto])
    np.testing.assert_allclose(corrected, expected, rtol=1e-4, atol=1e-4)


def test_scan_selection_uses_matching_interval(cal_ms_factory, bins, tmpdir):
    """With -startscan the solution intervals line up with calibrate's.

    ``calibrate -startscan 1 -t 1`` on a three scan MS gives two intervals,
    for scans 1 and 2. applysolutions with the same scan selection must use
    interval 0 for scan 1 and interval 1 for scan 2.
    """
    ms_path, _ = cal_ms_factory()
    scans = msgen.set_scans_per_timestep(ms_path)
    before = msgen.get_column(ms_path, "DATA")

    # Interval 0 is the identity, interval 1 doubles every gain (x4 in power)
    bandpass = np.zeros((2, 36, 288, 4), dtype=np.complex128)
    bandpass[..., 0] = bandpass[..., 3] = 1.0
    bandpass[1] *= 2.0
    sol_path = AOSolutions(
        path=Path(tmpdir) / "sols.bin",
        nsol=2,
        nant=36,
        nchan=288,
        npol=4,
        bandpass=bandpass,
    ).save(Path(tmpdir) / "sols.bin")

    run(
        bins.applysolutions,
        applysolutions_args(ms_path, sol_path, copy=False, extra=["-startscan", "1"]),
    )
    after = msgen.get_column(ms_path, "DATA")
    ant1, ant2 = msgen.antennas(ms_path)
    cross = ant1 != ant2

    for scan, factor in ((0, 1.0), (1, 1.0), (2, 4.0)):
        rows = cross & (scans == scan)
        np.testing.assert_allclose(after[rows], factor * before[rows], rtol=1e-5)


# Channel blocks (-ch), memory (-absmem) and timestep range (-interval)


def block_bandpass(block_values: np.ndarray, nant: int = 36) -> np.ndarray:
    """Scalar solutions per channel block, shape (1, nant, nblock, 4)"""
    nblock = len(block_values)
    bandpass = np.zeros((1, nant, nblock, 4), dtype=np.complex128)
    bandpass[..., 0] = block_values
    bandpass[..., 3] = block_values
    return bandpass


def channel_blocks(nchan: int, nblock: int) -> np.ndarray:
    """Block index of every channel, laid out as calibrate does"""
    blocks = np.zeros(nchan, dtype=int)
    for cb in range(nblock):
        blocks[cb * nchan // nblock : (cb + 1) * nchan // nblock] = cb
    return blocks


def test_ch_1_matches_default(cal_ms_factory, bins, tmpdir):
    """-ch 1 is the default: one solution per channel"""
    ms_path, _ = cal_ms_factory(noise_sigma=np.full(36, 1.0))
    default = Path(tmpdir) / "default.bin"
    ch1 = Path(tmpdir) / "ch1.bin"
    run(bins.calibrate, calibrate_args(ms_path, default))
    run(bins.calibrate, calibrate_args(ms_path, ch1, extra=["-ch", "1"]))

    assert ch1.read_bytes() == default.read_bytes()


def test_channel_blocks_recover_gains(cal_ms_factory, bins, tmpdir):
    """-ch 4 solves one Jones matrix per 4 channels and applysolutions uses it"""
    ms_path, truth = cal_ms_factory(channel_block=4)
    sol_path = Path(tmpdir) / "sols.bin"
    run(bins.calibrate, calibrate_args(ms_path, sol_path, extra=["-ch", "4"]))

    sols = AOSolutions.load(sol_path)
    assert (sols.nsol, sols.nant, sols.nchan, sols.npol) == (1, 36, 72, 4)
    ant1, ant2 = msgen.antennas(ms_path)
    per_channel = np.repeat(sols.bandpass[0], 4, axis=1)
    error = msgen.gain_product_error(per_channel, truth.jones, ant1, ant2)
    assert np.max(error) < GAIN_TOLERANCE

    run(bins.applysolutions, applysolutions_args(ms_path, sol_path))
    output = msgen.get_column(ms_path, "CORRECTED_DATA")
    model = msgen.get_column(ms_path, "MODEL_DATA")
    cross = ant1 != ant2
    assert_restores_model(output[cross], model[cross])


def test_uneven_channel_blocks(cal_ms_factory, bins, tmpdir):
    """-ch 5 on 288 channels gives 57 blocks of 5 or 6 channels"""
    ms_path, _ = cal_ms_factory(channel_block=288)
    sol_path = Path(tmpdir) / "sols.bin"
    run(bins.calibrate, calibrate_args(ms_path, sol_path, extra=["-ch", "5"]))
    assert AOSolutions.load(sol_path).nchan == 288 // 5

    run(bins.applysolutions, applysolutions_args(ms_path, sol_path))
    ant1, ant2 = msgen.antennas(ms_path)
    cross = ant1 != ant2
    assert_restores_model(
        msgen.get_column(ms_path, "CORRECTED_DATA")[cross],
        msgen.get_column(ms_path, "MODEL_DATA")[cross],
    )


@pytest.mark.parametrize("nblock", [72, 57, 1])
def test_applysolutions_maps_channel_blocks(cal_ms_factory, bins, tmpdir, nblock):
    """Each channel gets the solution of the block it falls in"""
    ms_path, _ = cal_ms_factory()
    before = msgen.get_column(ms_path, "DATA")
    gains = 1.0 + np.arange(nblock) / nblock
    sol_path = AOSolutions(
        path=Path(tmpdir) / "sols.bin",
        nsol=1,
        nant=36,
        nchan=nblock,
        npol=4,
        bandpass=block_bandpass(gains),
    ).save(Path(tmpdir) / "sols.bin")

    run(bins.applysolutions, applysolutions_args(ms_path, sol_path))
    after = msgen.get_column(ms_path, "CORRECTED_DATA")

    factor = gains[channel_blocks(288, nblock)] ** 2
    np.testing.assert_allclose(after, before * factor[None, :, None], rtol=1e-5)


def test_applysolutions_rejects_too_many_channels(cal_ms_factory, bins, tmpdir):
    """A solutions file with more channels than the MS is an error"""
    ms_path, _ = cal_ms_factory()
    sol_path = AOSolutions(
        path=Path(tmpdir) / "sols.bin",
        nsol=1,
        nant=36,
        nchan=289,
        npol=4,
        bandpass=block_bandpass(np.ones(289)),
    ).save(Path(tmpdir) / "sols.bin")

    result = run(
        bins.applysolutions, applysolutions_args(ms_path, sol_path), check=False
    )
    assert result.returncode != 0


@pytest.mark.parametrize("extra", [[], ["-ch", "4"]], ids=["channels", "blocks"])
def test_absmem_passes_give_same_solutions(cal_ms_factory, bins, tmpdir, extra):
    """A small -absmem splits the band into passes without changing the result.

    This also checks that each pass predicts its own channel range.
    """
    ms_path, _ = cal_ms_factory(noise_sigma=np.full(36, 1.0))
    one_pass = Path(tmpdir) / "one_pass.bin"
    passes = Path(tmpdir) / "passes.bin"
    run(bins.calibrate, calibrate_args(ms_path, one_pass, extra=extra))
    result = run(
        bins.calibrate,
        calibrate_args(ms_path, passes, extra=[*extra, "-absmem", "0.002"]),
    )

    pass_count = int(re.search(r"\((\d+) passes\)", result.stdout).group(1))
    assert pass_count > 1
    assert passes.read_bytes() == one_pass.read_bytes()


def test_interval_uses_selected_timesteps(cal_ms_factory, bins, tmpdir):
    """-interval 0 2 ignores the corrupted third timestep"""
    ms_path, truth = cal_ms_factory()
    time = msgen.get_column(ms_path, "TIME")
    data = msgen.get_column(ms_path, "DATA")
    last = time == np.max(time)
    rng = np.random.default_rng(7)
    data[last] = 10 * (rng.normal(size=data[last].shape) + 1j)
    msgen.put_column(ms_path, "DATA", data)
    ant1, ant2 = msgen.antennas(ms_path)

    errors = {}
    for label, extra in (("all", []), ("interval", ["-interval", "0", "2"])):
        sol_path = Path(tmpdir) / f"{label}.bin"
        run(bins.calibrate, calibrate_args(ms_path, sol_path, extra=extra))
        sols = AOSolutions.load(sol_path)
        errors[label] = np.max(
            msgen.gain_product_error(sols.bandpass[0], truth.jones, ant1, ant2)
        )

    assert errors["interval"] < GAIN_TOLERANCE
    assert errors["all"] > 10 * GAIN_TOLERANCE


def test_interval_with_solution_intervals(cal_ms_factory, bins, tmpdir):
    """-interval 1 3 -t 1 gives one solution per selected timestep"""
    ms_path, truth = cal_ms_factory()
    sol_path = Path(tmpdir) / "sols.bin"
    run(
        bins.calibrate,
        calibrate_args(ms_path, sol_path, extra=["-interval", "1", "3", "-t", "1"]),
    )

    sols = AOSolutions.load(sol_path)
    assert sols.nsol == 2
    ant1, ant2 = msgen.antennas(ms_path)
    for interval in range(2):
        error = msgen.gain_product_error(
            sols.bandpass[interval], truth.jones, ant1, ant2
        )
        assert np.max(error) < GAIN_TOLERANCE


def test_applysolutions_interval(cal_ms_factory, bins, tmpdir):
    """applysolutions -interval 1 3 maps solution intervals like calibrate"""
    ms_path, _ = cal_ms_factory()
    before = msgen.get_column(ms_path, "DATA")
    time = msgen.get_column(ms_path, "TIME")
    _, timestep = np.unique(time, return_inverse=True)

    bandpass = np.concatenate(
        [block_bandpass(np.ones(288)), 2 * block_bandpass(np.ones(288))]
    )
    sol_path = AOSolutions(
        path=Path(tmpdir) / "sols.bin",
        nsol=2,
        nant=36,
        nchan=288,
        npol=4,
        bandpass=bandpass,
    ).save(Path(tmpdir) / "sols.bin")

    run(
        bins.applysolutions,
        applysolutions_args(
            ms_path, sol_path, copy=False, extra=["-interval", "1", "3"]
        ),
    )
    after = msgen.get_column(ms_path, "DATA")
    for step, factor in ((0, 1.0), (1, 1.0), (2, 4.0)):
        rows = timestep == step
        np.testing.assert_allclose(after[rows], factor * before[rows], rtol=1e-5)
