"""Bugs found while comparing this fork with upstream mwa-reduce.

Each test here describes the correct behaviour and is marked as an expected
failure against the current code. The commit that fixes a bug removes its
``xfail`` marker and moves the test next to the related tests, so a fix that
does not work fails the suite (``strict=True``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tests import msgen
from tests.aosolutions import AOSolutions
from tests.runner import addmodel_args, applysolutions_args, calibrate_args, run

NANT = 36


@pytest.mark.xfail(
    strict=True,
    reason="The model frequency is 0/0 when the MS has one channel (2a)",
)
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


@pytest.mark.xfail(
    strict=True,
    reason="applysolutions writes the previous row into autocorrelations (2b)",
)
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


@pytest.mark.xfail(
    strict=True,
    reason="A short read of the solutions payload is not detected (2d)",
)
def test_truncated_solutions_file_is_rejected(cal_ms_factory, bins, tmpdir, ao_sols):
    """applysolutions refuses a solutions file that ends early"""
    ms_path, _ = cal_ms_factory()
    truncated = Path(tmpdir) / "truncated.bin"
    truncated.write_bytes(ao_sols.read_bytes()[:300000])

    result = run(
        bins.applysolutions, applysolutions_args(ms_path, truncated), check=False
    )
    assert result.returncode != 0


@pytest.mark.xfail(
    strict=True,
    reason="The solution file error does not say the file could not be read (2d)",
)
def test_missing_solutions_file_reports_the_file(cal_ms_factory, bins, tmpdir):
    """A missing solutions file gives an error naming the problem"""
    ms_path, _ = cal_ms_factory()
    missing = Path(tmpdir) / "missing.bin"

    result = run(
        bins.applysolutions, applysolutions_args(ms_path, missing), check=False
    )
    assert result.returncode != 0
    assert "missing.bin" in result.stderr


@pytest.mark.xfail(
    strict=True,
    reason="applysolutions counts timesteps outside -startscan/-endscan (2e)",
)
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
    bandpass = np.zeros((2, NANT, 288, 4), dtype=np.complex128)
    bandpass[..., 0] = bandpass[..., 3] = 1.0
    bandpass[1] *= 2.0
    sol_path = AOSolutions(
        path=Path(tmpdir) / "sols.bin",
        nsol=2,
        nant=NANT,
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


@pytest.mark.xfail(
    strict=True,
    reason="addmodel does not write autocorrelations into a new column (2h)",
)
def test_addmodel_new_column_autocorrelations_are_zero(cal_ms_factory, bins):
    """A column created by addmodel has zero autocorrelations, not garbage"""
    ms_path, _ = cal_ms_factory()
    run(bins.addmodel, addmodel_args(ms_path, datacolumn="NEW_MODEL", mode="c"))

    predicted = msgen.get_column(ms_path, "NEW_MODEL")
    ant1, ant2 = msgen.antennas(ms_path)
    np.testing.assert_array_equal(predicted[ant1 == ant2], 0)


@pytest.mark.xfail(
    strict=True,
    reason="Known limitation, not fixed: rows with ANTENNA1 > ANTENNA2 are "
    "used without conjugating the data (see visibilityarray.h)",
)
def test_reversed_baselines_recover_gains(cal_ms_factory, bins, tmpdir):
    """Rows stored as (a2, a1) with V^H are calibrated like (a1, a2) rows"""
    ms_path, truth = cal_ms_factory()
    ant1, _ = msgen.antennas(ms_path)
    msgen.reverse_baselines(ms_path, np.arange(len(ant1)) % 2 == 1)

    sol_path = Path(tmpdir) / "sols.bin"
    run(bins.calibrate, calibrate_args(ms_path, sol_path, model=None))

    sols = AOSolutions.load(sol_path)
    ant1, ant2 = msgen.antennas(ms_path)
    error = msgen.gain_product_error(sols.bandpass[0], truth.jones, ant1, ant2)
    assert np.max(error) < 1e-3
