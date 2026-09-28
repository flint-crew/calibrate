"""addmodel predicts the sky model into a column the way Flint uses it"""

from __future__ import annotations

import numpy as np
import pytest

from tests import msgen
from tests.runner import addmodel_args, run


@pytest.mark.parametrize("mode", ["a", "s", "c", "z"])
def test_addmodel_modes(cal_ms_factory, bins, mode):
    """Each mode combines DATA with the predicted model as documented"""
    ms_path, _ = cal_ms_factory()
    before = msgen.get_column(ms_path, "DATA")
    model = msgen.get_column(ms_path, "MODEL_DATA")
    ant1, ant2 = msgen.antennas(ms_path)
    cross = ant1 != ant2

    run(bins.addmodel, addmodel_args(ms_path, mode=mode))
    after = msgen.get_column(ms_path, "DATA")

    expected = {
        "a": before + model,
        "s": before - model,
        "c": model,
        "z": np.zeros_like(model),
    }[mode]
    scale = np.max(np.abs(model[cross]))
    np.testing.assert_allclose(after[cross], expected[cross], atol=1e-5 * scale)
    # addmodel never predicts autocorrelations, so an existing column keeps them
    np.testing.assert_array_equal(after[~cross], before[~cross])


def test_addmodel_new_column(cal_ms_factory, bins):
    """A column that does not exist yet is created (Flint's MODEL_DATA use)"""
    ms_path, _ = cal_ms_factory()
    model = msgen.get_column(ms_path, "MODEL_DATA")
    ant1, ant2 = msgen.antennas(ms_path)
    cross = ant1 != ant2

    run(bins.addmodel, addmodel_args(ms_path, datacolumn="NEW_MODEL", mode="c"))
    predicted = msgen.get_column(ms_path, "NEW_MODEL")

    assert predicted.shape == model.shape
    np.testing.assert_array_equal(predicted[cross], model[cross])


def test_addmodel_new_column_autocorrelations_are_zero(cal_ms_factory, bins):
    """A column created by addmodel has zero autocorrelations, not garbage"""
    ms_path, _ = cal_ms_factory()
    run(bins.addmodel, addmodel_args(ms_path, datacolumn="NEW_MODEL", mode="c"))

    predicted = msgen.get_column(ms_path, "NEW_MODEL")
    ant1, ant2 = msgen.antennas(ms_path)
    np.testing.assert_array_equal(predicted[ant1 == ant2], 0)
