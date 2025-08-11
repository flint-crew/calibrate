# Build Instructions for Calibrate Project

## Prerequisites

- **C++ compiler** (g++ recommended)
- **CASA (Common Astronomy Software Applications)** libraries
- **GSL (GNU Scientific Library)**
- **Boost libraries** (thread, system)

## Quick Build

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

## Custom Casacore Location

If casacore is not in standard system locations:

```bash
cmake ../ -DCMAKE_PREFIX_PATH=/path/to/casacore
```

**Example for conda:**
```bash
cmake ../ -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/casacore
```

## Clean Rebuild

```bash
cd build
rm -rf *
cmake ../
make
```

## Troubleshooting

- **CASA libraries not found**: Use `-DCMAKE_PREFIX_PATH=/path/to/casacore`
- **GSL/Boost not found**: Install via package manager (`apt-get`, `brew`, etc.)
- **Build errors**: Ensure casacore version 3.6+ is installed
