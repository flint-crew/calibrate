# Calibrate

Forked from [André Offringa's MWA data reduction code](https://github.com/ICRAR/mwa-reduce).

If any of this code is used, please cite [Offringa et al. (2016)](https://doi.org/10.1093/mnras/stw310).

Colloquially referred to as `aocalibrate` in Flint.

Note that [MWA Hyperdrive](https://github.com/MWATelescope/mwa_hyperdrive) is a more modern, and well-supported implementation of these algorithems.

## Provenance

This code was copied from [ICRAR/mwa-reduce](https://github.com/ICRAR/mwa-reduce), written by André Offringa with contributions from Marcin Sokolowski, Natasha Hurley-Walker, Paul Hancock, Tim Galvin and others. The first commit here (`146fc39`, Emil Lenc, 2021) is that copy, without its history. The upstream history up to [`90b202d`](https://github.com/ICRAR/mwa-reduce/commit/90b202d2d6cc5ee25e7f18dd1fecc8c6bfebe3d3) (March 2021, the last upstream commit before the copy) has since been merged into this repository's history, so `git log 90b202d` shows it.

The copied files do not all come from the same upstream version: `calibrator.cpp`, `calibrationmethod.cpp` and `solutionapplier.h` match upstream from around 2015 (`83f8275`), `predicter` and `mspredicter` from October 2018, and `model/` from later upstream code. Later upstream commits have not been merged; fixes from them are ported individually.

## Installation

### Prerequisites

- **C++ compiler** (g++ recommended)
- **CASA (Common Astronomy Software Applications)** libraries
- **GSL (GNU Scientific Library)**
- **Boost libraries**

### Quick Build

```bash
# Create build directory
mkdir -p build
cd build

# Configure and build
cmake ../
make

# Install (optional)
sudo make install
```

### Custom Casacore Location

If casacore is not in standard system locations:

```bash
cmake ../ -DCMAKE_PREFIX_PATH=/path/to/casacore
```

**Example for conda:**
```bash
cmake ../ -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/casacore
```

### Clean Rebuild

```bash
cd build
rm -rf *
cmake ../
make
```

### Troubleshooting

- **CASA libraries not found**: Use `-DCMAKE_PREFIX_PATH=/path/to/casacore`
- **GSL/Boost not found**: Install via package manager (`apt-get`, `brew`, etc.)
- **Build errors**: Ensure casacore version 3.6+ is installed

## Usage

These tools perform "MitchCal":

> a direction-independent full-polarization self-calibration. This is performed with the mitchcal tool, which is the authors’ custom implementation of the algorithm described by Mitchell et al. (2008)

In the current form of this repo, only a single time-step, frequency dependent solution is computed. This is suitable for a bandpass calibration.

The CLI hooks are:

```bash
calibrate
# Usage: calibrate [-p <phases.txt> <gains.txt>] [-refmod <0|1|2> [-minuv <min uvw dist in m>] [-maxuv <min uvw dist in m>] [-startscan <scan>] [-endscan <scan>] [-a <min-accuracy> <stop-accuracy>] [-i <niter>] [-j <threads>] [-m <model>] [-scalar] [-diag] [-rhs <rhs solutions>] [-rotation] [-t timesteps] [-datacolumn <name>] [-quiet] <measurementset.ms> <solutions.bin>

# This will calculate "static" phase offsets for all stations. It produces approximate least-squares solutions.
# refmode=0 process all baselines; =1 only include baselines to reference antenna; =2 exclude baselines to reference antenna.
```

```bash
applysolutions
# Usage: applysolutions [-datacolumn <name>] [-gflag <solutions-flag-file.txt>] [-startscan <scan>] [-endscan <scan>] [-copy/-nocopy] [-s xx xy yx yy] <ms> <gains-bin-file>
# Will apply the found solution matrices.
# Options:
#   -copy/-nocopy Don't(/do) alter the original DATA column but store the corrected data in the CORRECTED_DATA (this is std CASA behaviour)
#     default: -copy
```

```bash
addmodel
# Usage: addmodel [-usemodelcol] [-datacolumn <COLUMN>] [-m <a|s|c|z>] [-n <σ>] <model> <ms>
# Modify visibilities using a model. If -usemodelcol is specified then the MODEL_DATA column is used as the source model otherwise the specified component model file is used. The modification to use is defined with the mode switch(-m) where a=add model to visibilities (default), s=subtract model from visibilities, c=copy model to visibilities, z=zero visibilities.
```