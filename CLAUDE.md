# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`hifigps` is a Python package for HI filament stacking simulation analysis (arxiv:2411.03988). It provides data processing, plotting, and calculation utilities for 21cm astronomy research.

## Development Commands

```bash
# Install package with test dependencies
pip install .[test]

# Run tests
pytest tests/

# Run tests with coverage
pytest tests/ --cov=hifigps --cov-report=xml

# Run a single test file
pytest tests/test_scripts.py

# Build package
pip install .
```

## Architecture

### Package Structure

```
src/hifigps/
├── __init__.py           # Package entry, exports all modules
├── scripts/              # CLI commands (hifigps-convolve, hifigps-stack, hifigps-find-fuzzy)
├── calculation.py         # Astronomical calculations (beamsize, sensitivity, redshift)
├── constant.py            # 21cm astronomy constants (rest_frequency, c_light, etc.)
├── data.py                # HDF5 I/O operations
├── plot.py               # Main plotting (heatmap, histogram, line, arcs)
├── plot_custom.py        # Extended/custom plotting functionality
├── stack.py              # Stacking procedure helpers
├── estimate*.py          # Signal estimation modules
├── halo*.py              # Halo component fitting/subtraction
├── bins.py                # Binning helper functions
└── utils.py              # Fitting utilities
```

### CLI Commands

All CLI scripts use `argparse`. Key commands:

- **hifigps-convolve**: Convolve map with FAST main beam. Positional: `map_file out_path nworker`
- **hifigps-find-fuzzy**: Find inner fuzzy particles in TNG simulations. Positional: `base snap`
- **hifigps-stack**: Galaxy pairwise stacking. Uses named args (see `--help`)

### Core Data Flow

1. **Stack workflow**: raw data → `data.py` I/O → `bins.py` processing → `stack.py` stacking → `plot.py` visualization
2. **Estimation workflow**: stacked results → `estimate*.py` signal estimation → error analysis
3. **Halo subtraction**: `halo_opt.py` (numerical) or `crafts_stack` (analytical)

## Key Patterns

- HDF5 data access via `hifigps.data.read_h5()` / `save_h5()`
- Masked arrays used throughout (`np.ma.masked_array`)
- Global `mask_map` variable used in multiprocessing pool for `hifigps-stack`
- Scripts are registered as console entry points in `pyproject.toml`
