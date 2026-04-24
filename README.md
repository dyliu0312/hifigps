![](example/img_tngcube_stack_fit_res.png)
# HIFiGPS (HI Filament Detection with Galaxy Pairwise Stacking)

![arXiv](https://img.shields.io/badge/arXiv-2411.03988-orange)
![Python](https://img.shields.io/badge/python-3.8+-blue)
![Tests](https://github.com/dyliu0312/hifigps/actions/workflows/test.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow)](https://opensource.org/licenses/MIT)

## Overview

Python toolkit for **HI filament stacking simulation** (arXiv:2411.03988) — includes data processing, plotting, and astronomical calculations.

## Modules

| Module | Description |
|--------|-------------|
| `data.py` | HDF5 file I/O (save, read) |
| `plot.py` | Figures (heatmap, histogram, line, arcs) |
| `constant.py` | 21cm astronomy constants |
| `calculation.py` | Beamsize, sensitivity, redshift calculations |
| `stack.py` | Stacking procedure helpers |
| `halo.py` / `halo_new.py` | Halo component fitting & subtraction |
| `utils.py` | Fitting utilities (coordinates, mask) |
| `estimate.py` / `estimate_fixwidth.py` | Signal level estimation |
| `bins.py` | Binning helper functions |

## CLI Commands

| Command | Description |
|---------|-------------|
| `hifigps-convolve` | Beam convolution (FAST main beam) |
| `hifigps-stack` | Galaxy pairwise stacking |
| `hifigps-find-fuzzy` | Find inner fuzzy particles for filament-only map construction |

## Installation

```sh
git clone https://github.com/dyliu0312/hifigps.git
cd hifigps
pip install .
```

or directly:
```sh
pip install git+https://github.com/dyliu0312/hifigps.git
```

## Dependencies

**Required:** numpy, matplotlib, h5py, astropy, scipy

**Optional:** tqdm (for `hifigps-stack`), [illustris_python](https://github.com/illustristng/illustris_python.git) (for `hifigps-find-fuzzy`)

Tested on Python ≥ 3.8.

## Usage

```py
from hifigps.calculation import freq2z, u
freq2z(1.3*u.GHz)
```

## Examples

- **Filament signal estimation:** [example](example) folder
- **Functionality tutorials:** [tutorial](tutorial) folder
- **Running `hifigps-stack`:** See [pair_stack.sh](slurm/pair_stack.sh) for parameters
