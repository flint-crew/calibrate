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
from tests.runner import calibrate_args, run


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
