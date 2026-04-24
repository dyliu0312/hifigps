"""
hifigps - A collection of python tools used in my HI filament stacking simulation work, arxiv:2411.03988.
"""

__version__ = "0.2.0"
__author__ = "dyliu <dyliu0312@gmail.com>"


from . import (
    bins,
    calculation,
    constant,
    data,
    estimate,
    estimate_fixwidth,
    halo,
    halo_opt,
    plot,
    plot_custom,
    stack,
    utils,
)

__all__ = [
    "bins",
    "calculation",
    "constant",
    "data",
    "estimate",
    "estimate_fixwidth",
    "halo",
    "halo_opt",
    "plot",
    "plot_custom",
    "stack",
    "utils",
]
