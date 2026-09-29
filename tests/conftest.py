from __future__ import annotations

import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from tests import msgen
from tests.runner import Binaries, find_binaries, run


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


@pytest.fixture(scope="session")
def bins() -> Binaries:
    """The calibrate, applysolutions and addmodel programs under test.

    Taken from ``$CALIBRATE_BIN_DIR``, then ``<repo>/build``, then ``$PATH``.
    """
    env_dir = os.environ.get("CALIBRATE_BIN_DIR")
    search = [Path(env_dir)] if env_dir else []
    search.append(Path(__file__).parent.parent / "build")
    binaries = find_binaries(search, use_path=True)
    if binaries is None:
        pytest.skip("calibrate binaries not found, set CALIBRATE_BIN_DIR")
    return binaries


@pytest.fixture(scope="session")
def baseline_bins() -> Binaries:
    """A reference build to compare against, from ``$BASELINE_BIN_DIR``."""
    env_dir = os.environ.get("BASELINE_BIN_DIR")
    if not env_dir:
        pytest.skip("BASELINE_BIN_DIR not set")
    binaries = find_binaries([Path(env_dir)], use_path=False)
    if binaries is None:
        pytest.fail(f"No calibrate binaries found in BASELINE_BIN_DIR={env_dir}")
    return binaries


@pytest.fixture
def ms_example(tmpdir) -> Path:
    """The Flint example measurement set, extracted unmodified"""
    return msgen.unpack_example_ms(Path(tmpdir))


@pytest.fixture
def ao_sols(tmpdir) -> Path:
    """A real solutions file from a 1934-638 bandpass (from Flint)"""
    src = msgen.DATA_DIR / "SB39433.B1934-638.beam0.calibrate.bin"
    out = Path(tmpdir) / src.name
    shutil.copyfile(src, out)
    return out


@pytest.fixture
def ao_sols_known_bad(tmpdir) -> Path:
    """A real solutions file with many NaN solutions (from Flint)"""
    src = msgen.DATA_DIR / "SB38969.B1934-638.beam35.aocalibrate.bin"
    out = Path(tmpdir) / src.name
    shutil.copyfile(src, out)
    return out


@pytest.fixture(scope="session")
def template_ms(tmp_path_factory, bins) -> Path:
    """The example measurement set with UVWs and a predicted MODEL_DATA.

    This is shared by the whole session, so tests must copy it before
    changing anything (``cal_ms_factory`` does this).
    """
    out_dir = Path(tmp_path_factory.mktemp("template"))
    ms_path = msgen.prepare_template_ms(out_dir)
    run(
        bins.addmodel,
        ["-datacolumn", "MODEL_DATA", "-m", "c", str(msgen.MODEL_1934), str(ms_path)],
    )
    return ms_path


@pytest.fixture(scope="session")
def template_model(template_ms) -> np.ndarray:
    """MODEL_DATA of the template measurement set"""
    return msgen.get_column(template_ms, "MODEL_DATA")


@pytest.fixture
def cal_ms_factory(tmpdir, template_ms, template_model):
    """Factory making calibration measurement sets with known Jones matrices.

    Keyword arguments are passed to ``msgen.make_calibration_ms``;
    ``model_data`` defaults to the template's MODEL_DATA. Returns the
    path of the new measurement set and its ``Truth``.
    """
    counter = iter(range(1000))

    def _make(name: str | None = None, **kwargs) -> tuple[Path, msgen.Truth]:
        name = name or f"cal{next(counter)}"
        output_ms = Path(tmpdir) / f"{name}.ms"
        kwargs.setdefault("model_data", template_model)
        truth = msgen.make_calibration_ms(
            template_ms=template_ms, output_ms=output_ms, **kwargs
        )
        return output_ms, truth

    return _make
