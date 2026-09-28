# Tests

The test suite runs the built `calibrate`, `applysolutions` and `addmodel`
programs on small measurement sets whose correct answer is known. The helper
code and command lines follow Flint (`flint.calibrate.aocalibrate` and
`flint.predict.addmodel`), so the tests exercise the programs the way the
pipeline does.

## Requirements

- A build of this repository (see the main `README.md`)
- Python 3.11 or newer with `python-casacore>=3.6`, `numpy>=2` and `pytest`

```bash
pip install "python-casacore>=3.6" "numpy>=2" pytest
```

## Running

From the repository root:

```bash
# Uses the programs in ./build, or $PATH if they are not there
pytest

# Or point at a specific build directory
CALIBRATE_BIN_DIR=/path/to/build pytest
```

`ctest` in the build directory runs the same suite against that build.

### Comparing with a reference build

`test_regression.py` compares every output byte for byte with another build.
It is skipped unless `BASELINE_BIN_DIR` points at the reference programs:

```bash
git worktree add ../calibrate-main main
cmake -S ../calibrate-main -B ../build-main && make -C ../build-main
BASELINE_BIN_DIR=../build-main CALIBRATE_BIN_DIR=build pytest
```

Run it before merging a change: existing options must keep giving identical
solutions and corrected data.

## What is tested

| File | Checks |
| --- | --- |
| `test_aosolutions.py` | The `.bin` solutions format, read with Flint's reader; applying real Flint solution files |
| `test_calibrate.py` | Known Jones matrices are recovered (full-Jones, `-diag`, `-t`, `-ch`, `-interval`, `-absmem`, noisy data); weighting and `-weightcolumn`; `applysolutions` restores the model for 4, 2 and 1 polarisations, including autocorrelations |
| `test_addmodel.py` | The add, subtract, copy and zero modes, and creating a new column |
| `test_regression.py` | Same output as `BASELINE_BIN_DIR` for all existing options (byte-identical with `REGRESSION_EXACT=1`, see the module docstring); for fully reversed MSs, solutions and XX/YY match and `.bin` files work across versions |
| `test_reversed_baselines.py` | MSs with every row ANTENNA1 > ANTENNA2 (e.g. SKA-Low): solutions in main's conjugate convention, all four correlations corrected, leakage solved; mixed row orders are rejected |

## Test data

`tests/data` holds files copied from Flint's test data, see
`tests/data/README.md`. `tests/msgen.py` turns the example measurement set into
calibration test cases: it computes UVWs, predicts the 1934-638 model with
`addmodel`, corrupts it with seeded random Jones matrices and optionally adds
noise with matching WEIGHT and SIGMA columns.
