import os
import subprocess
import numpy as np
import pytest
from hifigps.data import save_h5, read_h5

def test_cli_convolve(tmp_path):
    # Prepare dummy data
    data_path = tmp_path / "test_data.h5"
    out_dir = tmp_path / "output"
    os.makedirs(out_dir, exist_ok=True)

    test_map = np.random.randn(2, 100, 100).astype(np.float32)
    save_h5(str(data_path), ["T"], [test_map])

    # Run CLI command
    cmd = ["hifigps-convolve", str(data_path), str(out_dir) + "/", "2"]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0
    # The output filename should be: out_dir / test_data_convolvedFASTbeam.hdf5
    expected_out = out_dir / "test_data_convolvedFASTbeam.hdf5"
    assert os.path.exists(expected_out)

    # Verify result
    con_data = read_h5(str(expected_out), "T")
    assert con_data.shape == test_map.shape

def test_cli_stack(tmp_path):
    # Prepare dummy map and pair catalog
    map_path = tmp_path / "map.h5"
    pair_path = tmp_path / "pairs.h5"
    out_dir = tmp_path / "stack_out"
    os.makedirs(out_dir, exist_ok=True)

    # 1. Save dummy map with bins
    # Shape: (freq, x, y)
    map_data = np.random.randn(10, 20, 20).astype(np.float32)
    f_bins = np.linspace(1.0, 1.4, 11)
    x_bins = np.linspace(-5, 5, 21)
    y_bins = np.linspace(-5, 5, 21)
    mask = np.zeros_like(map_data, dtype=bool)

    save_h5(str(map_path),
            ["T", "mask", "f_bin_edge", "x_bin_edge", "y_bin_edge"],
            [map_data, mask, f_bins, x_bins, y_bins])

    # 2. Save dummy pair catalog
    # pos: (N, 6) -> [ra1, dec1, freq1, ra2, dec2, freq2]
    # We use 2 pairs
    pos = np.zeros((2, 6))
    pos[:, 0] = 0.0 # ra1
    pos[:, 1] = 0.0 # dec1
    pos[:, 2] = 1.2 # freq1
    pos[:, 3] = 1.0 # ra2
    pos[:, 4] = 1.0 # dec2
    pos[:, 5] = 1.2 # freq2
    is_ra = np.array([True, False])

    save_h5(str(pair_path), ["is_ra", "pos"], [is_ra, pos])

    # Run CLI command with command-line arguments
    cmd = [
        "hifigps-stack",
        "--map-base", str(tmp_path),
        "--map-prefix", "map.h5",
        "--paircat-base", str(tmp_path),
        "--paircat-prefix", "pairs.h5",
        "--out-base", str(out_dir),
        "--out-prefix", "stack_res.h5",
        "--nfs", "2",
        "--ssize", "1",
        "--nworker", "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"stderr: {result.stderr}"
    expected_out = out_dir / "stack_res.h5"
    assert os.path.exists(expected_out)
