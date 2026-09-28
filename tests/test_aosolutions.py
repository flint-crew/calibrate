"""The binary solutions file, read and written the way Flint does it"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from tests import msgen
from tests.aosolutions import HEADER_INTRO, HEADER_SIZE, AOSolutions
from tests.runner import applysolutions_args, calibrate_args, run


def test_load_flint_solutions(ao_sols):
    """A solutions file made by calibrate for Flint loads with Flint's reader"""
    sols = AOSolutions.load(ao_sols)

    assert (sols.nsol, sols.nant, sols.nchan, sols.npol) == (1, 36, 288, 4)
    assert sols.bandpass.shape == (1, 36, 288, 4)
    assert sols.bandpass.dtype == np.complex128


def test_known_bad_solutions_contain_nans(ao_sols_known_bad):
    """The known bad Flint file keeps its NaN solutions"""
    sols = AOSolutions.load(ao_sols_known_bad)

    assert np.isnan(sols.bandpass).any()
    assert np.isfinite(sols.bandpass).any()


def test_roundtrip_is_byte_identical(ao_sols, tmpdir):
    """Writing a loaded file back out gives the same bytes"""
    sols = AOSolutions.load(ao_sols)
    out = sols.save(Path(tmpdir) / "roundtrip.bin")

    assert out.read_bytes() == ao_sols.read_bytes()


def test_calibrate_output_matches_flint_format(cal_ms_factory, bins, tmpdir):
    """calibrate writes a header and payload that Flint's reader understands"""
    ms_path, _ = cal_ms_factory()
    sol_path = Path(tmpdir) / "sols.bin"
    run(bins.calibrate, calibrate_args(ms_path, sol_path))

    raw = sol_path.read_bytes()
    assert raw[:8] == HEADER_INTRO
    sols = AOSolutions.load(sol_path)
    assert (sols.nsol, sols.nant, sols.nchan, sols.npol) == (1, 36, 288, 4)
    assert len(raw) == HEADER_SIZE + sols.bandpass.size * 16


def test_applysolutions_applies_flint_solutions(cal_ms_factory, bins, ao_sols):
    """applysolutions applies a real Flint solutions file as S_1 V S_2^H"""
    ms_path, _ = cal_ms_factory()
    run(bins.applysolutions, applysolutions_args(ms_path, ao_sols))

    sols = AOSolutions.load(ao_sols)
    data = msgen.get_column(ms_path, "DATA")
    corrected = msgen.get_column(ms_path, "CORRECTED_DATA")
    ant1, ant2 = msgen.antennas(ms_path)
    cross = ant1 != ant2

    expected = msgen.apply_solutions(
        data[cross], sols.bandpass[0], ant1[cross], ant2[cross]
    )
    finite = np.isfinite(expected)
    assert finite.any()
    np.testing.assert_allclose(
        corrected[cross][finite], expected[finite], rtol=1e-4, atol=1e-4
    )


def test_applysolutions_applies_known_bad_solutions(
    cal_ms_factory, bins, ao_sols_known_bad
):
    """NaN solutions give NaN data rather than an error"""
    ms_path, _ = cal_ms_factory()
    run(bins.applysolutions, applysolutions_args(ms_path, ao_sols_known_bad))

    corrected = msgen.get_column(ms_path, "CORRECTED_DATA")
    ant1, ant2 = msgen.antennas(ms_path)
    cross = ant1 != ant2
    assert np.isnan(corrected[cross]).any()
    assert np.isfinite(corrected[cross]).any()


def test_truncated_solutions_file_is_rejected(cal_ms_factory, bins, tmpdir, ao_sols):
    """applysolutions refuses a solutions file that ends early"""
    ms_path, _ = cal_ms_factory()
    truncated = Path(tmpdir) / "truncated.bin"
    truncated.write_bytes(ao_sols.read_bytes()[:300000])

    result = run(
        bins.applysolutions, applysolutions_args(ms_path, truncated), check=False
    )
    assert result.returncode != 0


def test_missing_solutions_file_reports_the_file(cal_ms_factory, bins, tmpdir):
    """A missing solutions file gives an error naming the problem"""
    ms_path, _ = cal_ms_factory()
    missing = Path(tmpdir) / "missing.bin"

    result = run(
        bins.applysolutions, applysolutions_args(ms_path, missing), check=False
    )
    assert result.returncode != 0
    assert "missing.bin" in result.stderr
