"""Build measurement sets with known answers for testing ``calibrate``.

Everything starts from the small ASKAP measurement set that Flint uses in its
own tests (``SB39400.RACS_0635-31.beam0.small.ms``: 36 antennas, 288
channels, 3 timesteps, 4 linear polarisations, autocorrelations included).
Its visibilities are all NaN and flagged, so it is only used as a template:
the helpers here write a model, corrupt it with known Jones matrices and set
the flag and weight columns.

All random quantities are drawn from seeded generators, so every measurement
set produced here is identical between runs.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import NamedTuple

import numpy as np
from casacore.tables import makearrcoldesc, maketabdesc, table

DATA_DIR = Path(__file__).parent / "data"
EXAMPLE_MS_ZIP = DATA_DIR / "SB39400.RACS_0635-31.beam0.small.ms.zip"
EXAMPLE_MS_NAME = "SB39400.RACS_0635-31.beam0.small.ms"
MODEL_1934 = DATA_DIR / "1934-638.calibrate.txt"


_VALUE_TYPES = {
    np.complex64: "complex",
    np.complex128: "dcomplex",
    np.float32: "float",
    np.float64: "double",
    np.bool_: "boolean",
}


class Truth(NamedTuple):
    """The known answer baked into a synthetic measurement set"""

    jones: np.ndarray
    """Per antenna and channel Jones matrices, shape (nant, nchan, 2, 2)"""
    noise_sigma: np.ndarray | None
    """Per antenna noise standard deviation, or None when noise free"""


def unpack_example_ms(output_dir: Path) -> Path:
    """Extract the Flint example measurement set as is.

    Args:
        output_dir (Path): Directory to extract into

    Returns:
        Path: Path to the extracted measurement set
    """
    shutil.unpack_archive(EXAMPLE_MS_ZIP, output_dir)
    return output_dir / EXAMPLE_MS_NAME


def prepare_template_ms(output_dir: Path) -> Path:
    """Extract the example measurement set and give it usable UVWs.

    Args:
        output_dir (Path): Directory to extract into

    Returns:
        Path: Path to the template measurement set
    """
    ms_path = unpack_example_ms(output_dir)
    compute_uvw(ms_path)
    return ms_path


def copy_ms(ms_path: Path, output_path: Path) -> Path:
    """Deep copy a measurement set to a new location.

    Args:
        ms_path (Path): Measurement set to copy
        output_path (Path): Destination

    Returns:
        Path: The destination
    """
    shutil.copytree(ms_path, output_path)
    return output_path


def random_jones(
    nant: int,
    nchan: int,
    seed: int = 1234,
    amp_scatter: float = 0.2,
    leakage: float = 0.05,
    diagonal: bool = False,
    channel_block: int = 1,
) -> np.ndarray:
    """Draw random antenna Jones matrices.

    Args:
        nant (int): Number of antennas
        nchan (int): Number of channels
        seed (int, optional): Seed for the generator. Defaults to 1234.
        amp_scatter (float, optional): Fractional scatter of the gain amplitudes. Defaults to 0.2.
        leakage (float, optional): Amplitude of the off-diagonal terms. Defaults to 0.05.
        diagonal (bool, optional): If True the off-diagonal terms are zero. Defaults to False.
        channel_block (int, optional): Keep the Jones constant over blocks of this many channels. Defaults to 1.

    Returns:
        np.ndarray: Jones matrices of shape (nant, nchan, 2, 2)
    """
    rng = np.random.default_rng(seed)
    nblock = (nchan + channel_block - 1) // channel_block
    shape = (nant, nblock)

    jones = np.zeros((nant, nblock, 2, 2), dtype=np.complex128)
    for p in (0, 1):
        amp = 1.0 + amp_scatter * rng.uniform(-1, 1, size=shape)
        phase = rng.uniform(-np.pi, np.pi, size=shape)
        jones[..., p, p] = amp * np.exp(1j * phase)
    if not diagonal:
        for p, q in ((0, 1), (1, 0)):
            jones[..., p, q] = leakage * (
                rng.normal(size=shape) + 1j * rng.normal(size=shape)
            )

    return np.repeat(jones, channel_block, axis=1)[:, :nchan]


def corrupt(
    vis: np.ndarray, jones: np.ndarray, ant1: np.ndarray, ant2: np.ndarray
) -> np.ndarray:
    """Compute J_1 V J_2^H for every row and channel.

    Args:
        vis (np.ndarray): Visibilities, shape (nrow, nchan, 4)
        jones (np.ndarray): Jones matrices, shape (nant, nchan, 2, 2)
        ant1 (np.ndarray): First antenna of each row
        ant2 (np.ndarray): Second antenna of each row

    Returns:
        np.ndarray: Corrupted visibilities, shape (nrow, nchan, 4)
    """
    nrow, nchan, _ = vis.shape
    v = vis.reshape(nrow, nchan, 2, 2)
    j1 = jones[ant1]
    j2h = np.conj(np.swapaxes(jones[ant2], -1, -2))
    return (j1 @ v @ j2h).reshape(nrow, nchan, 4)


def apply_solutions(
    vis: np.ndarray,
    solutions: np.ndarray,
    ant1: np.ndarray,
    ant2: np.ndarray,
) -> np.ndarray:
    """Apply solutions from a calibrate file the way ``applysolutions`` does.

    The file stores the inverse of each antenna Jones, so this is simply
    ``corrupt`` with the stored matrices.

    Args:
        vis (np.ndarray): Visibilities, shape (nrow, nchan, 4)
        solutions (np.ndarray): One interval of solutions, shape (nant, nchan, 4)
        ant1 (np.ndarray): First antenna of each row
        ant2 (np.ndarray): Second antenna of each row

    Returns:
        np.ndarray: Corrected visibilities
    """
    nant, nchan, _ = solutions.shape
    return corrupt(vis, solutions.reshape(nant, nchan, 2, 2), ant1, ant2)


def add_array_column(
    ms_path: Path,
    name: str,
    value: np.ndarray,
    define_cells: bool = True,
) -> None:
    """Add (or replace) an array column and fill it.

    Args:
        ms_path (Path): Measurement set to modify
        name (str): Column name
        value (np.ndarray): Values for every row, shape (nrow, ...)
        define_cells (bool, optional): If False the column is declared, without
            a fixed shape, but no cell is written, as some writers do. Defaults to True.
    """
    with table(str(ms_path), readonly=False, ack=False) as tab:
        if name in tab.colnames():
            tab.removecols(name)
        shape_kwargs = {"shape": list(value.shape[1:]), "options": 4}  # FixedShape
        coldesc = makearrcoldesc(
            name,
            0,
            valuetype=_VALUE_TYPES[value.dtype.type],
            ndim=value.ndim - 1,
            **(shape_kwargs if define_cells else {}),
        )
        dminfo = {
            "TYPE": "TiledShapeStMan",
            "NAME": f"Tiled_{name}",
            "SPEC": {"DEFAULTTILESHAPE": [*value.shape[:0:-1], 64]},
        }
        tab.addcols(maketabdesc(coldesc), dminfo)
        if define_cells:
            tab.putcol(name, value)


def remove_columns(ms_path: Path, names: list[str]) -> None:
    """Remove columns from a measurement set, ignoring any that are missing.

    Args:
        ms_path (Path): Measurement set to modify
        names (list[str]): Columns to remove
    """
    with table(str(ms_path), readonly=False, ack=False) as tab:
        existing = [n for n in names if n in tab.colnames()]
        if existing:
            tab.removecols(existing)


def get_column(ms_path: Path, name: str) -> np.ndarray:
    """Read a whole column.

    Args:
        ms_path (Path): Measurement set to read
        name (str): Column name

    Returns:
        np.ndarray: The column values
    """
    with table(str(ms_path), ack=False) as tab:
        return tab.getcol(name)


def put_column(ms_path: Path, name: str, value: np.ndarray) -> None:
    """Write a whole existing column.

    Args:
        ms_path (Path): Measurement set to write
        name (str): Column name
        value (np.ndarray): Values for every row
    """
    with table(str(ms_path), readonly=False, ack=False) as tab:
        tab.putcol(name, value)


def antennas(ms_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read the ANTENNA1 and ANTENNA2 columns.

    Args:
        ms_path (Path): Measurement set to read

    Returns:
        tuple[np.ndarray, np.ndarray]: The two antenna columns
    """
    return get_column(ms_path, "ANTENNA1"), get_column(ms_path, "ANTENNA2")


def antenna_count(ms_path: Path) -> int:
    """Number of rows in the ANTENNA table.

    Args:
        ms_path (Path): Measurement set to read

    Returns:
        int: Number of antennas
    """
    with table(str(ms_path / "ANTENNA"), ack=False) as tab:
        return tab.nrows()


MAIN_ARRAY_COLUMNS = (
    "DATA",
    "MODEL_DATA",
    "CORRECTED_DATA",
    "FLAG",
    "WEIGHT_SPECTRUM",
    "SIGMA_SPECTRUM",
)
MAIN_POL_COLUMNS = ("WEIGHT", "SIGMA")


def reshape_ms(
    ms_path: Path,
    output_path: Path,
    channels: slice | None = None,
    correlations: list[int] | None = None,
) -> Path:
    """Copy a measurement set keeping a subset of channels and correlations.

    Args:
        ms_path (Path): Measurement set to copy
        output_path (Path): Destination
        channels (slice | None, optional): Channels to keep. Defaults to all.
        correlations (list[int] | None, optional): Indices of the correlations to keep. Defaults to all.

    Returns:
        Path: The new measurement set
    """
    channels = channels if channels is not None else slice(None)
    correlations = correlations if correlations is not None else [0, 1, 2, 3]
    copy_ms(ms_path, output_path)

    with table(str(output_path), ack=False) as tab:
        present = [
            c for c in (*MAIN_ARRAY_COLUMNS, *MAIN_POL_COLUMNS) if c in tab.colnames()
        ]
        values = {c: tab.getcol(c) for c in present}
        keywords = {c: tab.getcolkeywords(c) for c in (*present, "FLAG_CATEGORY")}
    # FLAG_CATEGORY is required by casacore but unused: recreate it empty so
    # its declared shape can not clash with the new data shape
    remove_columns(output_path, ["FLAG_CATEGORY"])
    with table(str(output_path), readonly=False, ack=False) as tab:
        tab.addcols(
            maketabdesc(makearrcoldesc("FLAG_CATEGORY", False, ndim=3)),
            {"TYPE": "StandardStMan", "NAME": "FlagCategory"},
        )
    for name in present:
        value = values[name]
        if name in MAIN_POL_COLUMNS:
            value = value[:, correlations]
        else:
            value = value[:, channels][:, :, correlations]
        add_array_column(output_path, name, np.ascontiguousarray(value))

    with table(str(output_path), readonly=False, ack=False) as tab:
        for name, kws in keywords.items():
            if kws:
                tab.putcolkeywords(name, kws)

    with table(str(output_path / "SPECTRAL_WINDOW"), readonly=False, ack=False) as tab:
        freq = tab.getcell("CHAN_FREQ", 0)[channels]
        for name in ("CHAN_FREQ", "CHAN_WIDTH", "EFFECTIVE_BW", "RESOLUTION"):
            tab.putcell(name, 0, tab.getcell(name, 0)[channels])
        tab.putcell("NUM_CHAN", 0, len(freq))
        tab.putcell("REF_FREQUENCY", 0, float(freq[0]))
        tab.putcell("TOTAL_BANDWIDTH", 0, float(np.sum(tab.getcell("CHAN_WIDTH", 0))))

    with table(str(output_path / "POLARIZATION"), readonly=False, ack=False) as tab:
        tab.putcell("CORR_TYPE", 0, tab.getcell("CORR_TYPE", 0)[correlations])
        tab.putcell("CORR_PRODUCT", 0, tab.getcell("CORR_PRODUCT", 0)[correlations])
        tab.putcell("NUM_CORR", 0, len(correlations))

    return output_path


def compute_uvw(ms_path: Path) -> np.ndarray:
    """Fill the UVW column from the antenna positions and phase centre.

    The example measurement set has its UVW column zeroed. This computes
    approximate UVWs using the Earth rotation angle as the sidereal time and
    no precession, nutation or aberration. They are not astrometrically
    exact, but ``addmodel`` and ``calibrate`` both predict from this column, so
    the test data stays self consistent.

    Args:
        ms_path (Path): Measurement set to modify

    Returns:
        np.ndarray: The UVW values written, shape (nrow, 3)
    """
    with table(str(ms_path / "ANTENNA"), ack=False) as tab:
        positions = tab.getcol("POSITION")
    with table(str(ms_path / "FIELD"), ack=False) as tab:
        ra, dec = tab.getcol("PHASE_DIR")[0, 0]

    ant1, ant2 = antennas(ms_path)
    time = get_column(ms_path, "TIME")
    baseline = positions[ant2] - positions[ant1]

    days = time / 86400.0 - 51544.5
    era = 2.0 * np.pi * (0.7790572732640 + 1.00273781191135448 * days)
    hour_angle = era - ra

    sin_h, cos_h = np.sin(hour_angle), np.cos(hour_angle)
    sin_d, cos_d = np.sin(dec), np.cos(dec)
    lx, ly, lz = baseline.T
    uvw = np.stack(
        [
            sin_h * lx + cos_h * ly,
            -sin_d * cos_h * lx + sin_d * sin_h * ly + cos_d * lz,
            cos_d * cos_h * lx - cos_d * sin_h * ly + sin_d * lz,
        ],
        axis=1,
    )
    put_column(ms_path, "UVW", uvw)
    return uvw


def set_scans_per_timestep(ms_path: Path) -> np.ndarray:
    """Give each timestep its own SCAN_NUMBER (0, 1, 2, ...).

    Args:
        ms_path (Path): Measurement set to modify

    Returns:
        np.ndarray: The new scan numbers
    """
    time = get_column(ms_path, "TIME")
    _, scan = np.unique(time, return_inverse=True)
    scan = scan.astype(np.int32)
    put_column(ms_path, "SCAN_NUMBER", scan)
    return scan


def make_calibration_ms(
    template_ms: Path,
    output_ms: Path,
    model_data: np.ndarray,
    seed: int = 1234,
    diagonal: bool = False,
    channel_block: int = 1,
    noise_sigma: np.ndarray | None = None,
    noise_seed: int = 42,
) -> Truth:
    """Create a measurement set whose DATA is a corrupted copy of a model.

    ``DATA = J_1 MODEL J_2^H (+ noise)``. MODEL_DATA holds the uncorrupted
    model, all flags are cleared and WEIGHT / SIGMA are set consistently with
    the noise (unity when noise free).

    Args:
        template_ms (Path): Measurement set to copy (normally the output of ``unpack_example_ms``)
        output_ms (Path): Path of the measurement set to create
        model_data (np.ndarray): Model visibilities, shape (nrow, nchan, 4)
        seed (int, optional): Seed for the Jones matrices. Defaults to 1234.
        diagonal (bool, optional): Draw diagonal Jones matrices only. Defaults to False.
        channel_block (int, optional): Keep the Jones constant over this many channels. Defaults to 1.
        noise_sigma (np.ndarray | None, optional): Per antenna noise standard deviation. Defaults to None.
        noise_seed (int, optional): Seed for the noise. Defaults to 42.

    Returns:
        Truth: The Jones matrices and noise used
    """
    copy_ms(template_ms, output_ms)
    ant1, ant2 = antennas(output_ms)
    nrow, nchan, npol = model_data.shape

    # addmodel does not predict autocorrelations, so give them a simple
    # unpolarised value to keep every row finite
    model_data = np.array(model_data, dtype=np.complex128)
    model_data[ant1 == ant2] = np.array([10.0, 0.0, 0.0, 10.0])
    jones = random_jones(
        nant=antenna_count(output_ms),
        nchan=nchan,
        seed=seed,
        diagonal=diagonal,
        channel_block=channel_block,
    )
    data = corrupt(model_data, jones, ant1, ant2)

    sigma_row = np.ones((nrow, npol), dtype=np.float32)
    if noise_sigma is not None:
        rng = np.random.default_rng(noise_seed)
        # Noise on a baseline is the quadrature sum of its two antennas
        baseline_sigma = np.sqrt(
            0.5 * (noise_sigma[ant1] ** 2 + noise_sigma[ant2] ** 2)
        )
        noise = rng.normal(size=data.shape) + 1j * rng.normal(size=data.shape)
        data = data + noise * baseline_sigma[:, None, None] / np.sqrt(2.0)
        sigma_row = np.repeat(baseline_sigma[:, None], npol, axis=1).astype(np.float32)

    put_column(output_ms, "DATA", data.astype(np.complex64))
    add_array_column(output_ms, "MODEL_DATA", model_data.astype(np.complex64))
    put_column(output_ms, "FLAG", np.zeros(data.shape, dtype=bool))
    put_column(output_ms, "FLAG_ROW", np.zeros(nrow, dtype=bool))
    put_column(output_ms, "SIGMA", sigma_row)
    put_column(output_ms, "WEIGHT", (1.0 / sigma_row**2).astype(np.float32))

    return Truth(jones=jones, noise_sigma=noise_sigma)


def conjugate_transpose(vis: np.ndarray) -> np.ndarray:
    """Swap the correlation order XY <-> YX and conjugate (V -> V^H).

    Args:
        vis (np.ndarray): Visibilities, shape (..., 4)

    Returns:
        np.ndarray: V^H for every visibility
    """
    return np.conj(vis[..., [0, 2, 1, 3]])


def reverse_baselines(ms_path: Path, rows: np.ndarray) -> None:
    """Store some rows with ANTENNA1 > ANTENNA2, keeping them physically identical.

    Swapping the antennas of a row means its visibilities become V^H and its
    UVW changes sign.

    Args:
        ms_path (Path): Measurement set to modify
        rows (np.ndarray): Boolean mask of the rows to reverse
    """
    ant1, ant2 = antennas(ms_path)
    new1 = np.where(rows, ant2, ant1)
    new2 = np.where(rows, ant1, ant2)
    put_column(ms_path, "ANTENNA1", new1)
    put_column(ms_path, "ANTENNA2", new2)
    uvw = get_column(ms_path, "UVW")
    uvw[rows] *= -1
    put_column(ms_path, "UVW", uvw)
    for name in ("DATA", "MODEL_DATA"):
        vis = get_column(ms_path, name)
        vis[rows] = conjugate_transpose(vis[rows])
        put_column(ms_path, name, vis)


def set_weight_variant(ms_path: Path, variant: str) -> None:
    """Rewrite the weight columns of a measurement set.

    The values are derived from the existing SIGMA column, so every variant
    describes the same noise.

    Args:
        ms_path (Path): Measurement set to modify
        variant (str): One of ``weight`` (WEIGHT and SIGMA only, the ASKAP
            default), ``weight_spectrum``, ``sigma_spectrum`` or
            ``empty_weight_spectrum`` (the column exists but has no cells).
            WEIGHT and SIGMA are required MS columns, so they always stay.
    """
    sigma = get_column(ms_path, "SIGMA")
    nchan = get_column(ms_path, "FLAG").shape[1]
    sigma_spectrum = np.repeat(sigma[:, None, :], nchan, axis=1)

    if variant == "weight":
        return
    elif variant == "weight_spectrum":
        add_array_column(
            ms_path, "WEIGHT_SPECTRUM", (1.0 / sigma_spectrum**2).astype(np.float32)
        )
    elif variant == "sigma_spectrum":
        add_array_column(ms_path, "SIGMA_SPECTRUM", sigma_spectrum.astype(np.float32))
    elif variant == "empty_weight_spectrum":
        add_array_column(
            ms_path,
            "WEIGHT_SPECTRUM",
            np.ones(sigma_spectrum.shape, dtype=np.float32),
            define_cells=False,
        )
    else:
        raise ValueError(f"Unknown weight variant {variant}")


def gain_product_error(
    solutions: np.ndarray,
    jones: np.ndarray,
    ant1: np.ndarray,
    ant2: np.ndarray,
    exclude: set[int] | None = None,
) -> np.ndarray:
    """Per channel error of the recovered solutions.

    A full-Jones solution against an unpolarised model is only defined up to a
    common unitary matrix, so antenna Jones matrices cannot be compared
    directly. Instead this checks the gauge invariant product
    ``S_1 J_1 J_2^H S_2^H``, which is the identity for perfect solutions.

    Args:
        solutions (np.ndarray): One interval of solutions, shape (nant, nchan, 4)
        jones (np.ndarray): True Jones matrices, shape (nant, nchan, 2, 2)
        ant1 (np.ndarray): First antenna of each baseline to check
        ant2 (np.ndarray): Second antenna of each baseline to check
        exclude (set[int] | None, optional): Antennas to leave out. Defaults to None.

    Returns:
        np.ndarray: Median absolute deviation from identity per channel, shape (nchan,)
    """
    keep = ant1 != ant2
    if exclude:
        keep &= ~np.isin(ant1, list(exclude)) & ~np.isin(ant2, list(exclude))
    ant1 = ant1[keep]
    ant2 = ant2[keep]

    nant, nchan, _ = solutions.shape
    identity = np.broadcast_to(np.eye(2), (len(ant1), nchan, 2, 2)).reshape(
        len(ant1), nchan, 4
    )
    corrected = apply_solutions(
        corrupt(identity.astype(np.complex128), jones, ant1, ant2),
        solutions,
        ant1,
        ant2,
    )
    return np.median(np.abs(corrected - identity), axis=(0, 2))
