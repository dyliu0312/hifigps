"""
Simulate the beam smearing effect.
"""

import gc
import multiprocessing as mp
import sys
import time

from astropy.convolution import Gaussian2DKernel, convolve_fft
from hifigps.calculation import get_beam_npix
from hifigps.data import delete_files, get_filename, is_exist, np, read_h5, save_h5


def convolve(data, kernel, **kwargs) -> np.ndarray:
    """
    convolve data with kernel using `convolve_fft` with periodic boundary 'warp'.
    """
    return convolve_fft(data, kernel, boundary="warp", **kwargs)


def beam_kernel_fast(z=0.1):
    """
    get the kernel of FAST main beam, assuming a ideal Gaussian beam model.
    """
    return Gaussian2DKernel(get_beam_npix(z))


def convolve_save(data, kernel, path, key="T", dtype=np.float32, **kwargs) -> None:
    """
    get and save the convolution result.

    Args:
        data: the input array.
        kernel: the convolution kernel.
        path: the output file path.
        key: the save key of output file, default if 'T'.
        dtype: the dtype of output file, default is `np.float32` or f4.
        **kwargs: other keyword args for `convolve_fft`.
    """
    con_data = convolve(data, kernel, **kwargs).astype(dtype)
    save_h5(
        path,
        [
            key,
        ],
        [
            con_data,
        ],
    )


def reconstruct_map(files, output, key="T") -> None:
    """
    load the seprately convolved results, reconstruct the final convolved map.

    Args:
        files: the list of path of the seprately convolved results.
        output: the output path of final map.
        key: the save key of the final map.
    """
    reconstruct = np.stack([read_h5(fn, key) for fn in files])
    save_h5(
        output,
        [key],
        [
            reconstruct,
        ],
    )


if __name__ == "__main__":
    args = sys.argv

    if len(args) != 4:
        print("Usage: python convolve.py [map_file_path] [out_path] [nworker]\n")
        print(
            "Example: python convolve.py '/home/dyliu/data/test.hdf5' '/home/dyliu/output/' 4 "
        )
        sys.exit()

    # --- default args ---
    Z = 0.1  # redshift.
    KEY = "T"  # the key of MAP data.
    DTYPE = np.float32  # the dtype of the MAP data, you may modify this for higher preccision, but **beware of memory usage**.
    TEMP_PREFIX = "convolved_temp_out_"  # the prefix of temprary output file.
    OUT_SUFFIX = "_convolvedFASTbeam.hdf5"  # the suffix of final output file.
    WAIT_TIME = 60  # seconds, waiting for the multiprocessing result saving processes to be finished, to do further MAP reconstruction.

    # --- input args ---
    map_path, out_path, nworker = args[1:]
    print(f"load data from {map_path}")
    print(f"save result at {out_path}")
    print(f"set {nworker} workers")

    # --- pre process ---
    print("---processing start----")
    nworker = int(nworker)
    beam_kernel = beam_kernel_fast(Z)
    filename = get_filename(map_path)
    out_name = out_path + filename + OUT_SUFFIX
    if is_exist(out_name):
        print(f"This map is already processed, the output file is {out_name}!")
        sys.exit()

    # --- load the map into the memory, **beware of the MEMORY USAGE**---
    t0 = time.time()
    map_data = read_h5(map_path, KEY)
    print(f"Loading the MAP data successfully with {time.time() - t0} seconds.")

    # --- multiprocess ---
    t1 = time.time()
    temp_out_file = []
    with mp.Pool(processes=nworker) as pool:
        for i, map_slice in enumerate(map_data):
            temp_out = out_path + filename + "_" + TEMP_PREFIX + f"{i}.hdf5"
            temp_out_file.append(temp_out)
            pool.apply(convolve_save, (map_slice, beam_kernel, temp_out, KEY, DTYPE))
            # pool.apply_async(convolve_save, (map_slice, beam_kernel, temp_out, KEY, DTYPE)) # it fail to create file

    print(f"Convolution process finished with {time.time() - t1} seconds.")

    # --- wait the save processes fully finished---
    while True:
        if all(
            [is_exist(fn) for fn in temp_out_file]
        ):  # check exist of the temporary output files
            break
        time.sleep(WAIT_TIME)
        print("Waiting for the temporary output files to be saved")

    # --- del MAP data and collect space---
    del map_data
    gc.collect()

    # --- reconstruct the sky map, **beware of this process require DOUBLE memory space of MAP data**---
    t2 = time.time()
    reconstruct_map(temp_out_file, out_name)
    print(f"Save convolved results finished with {time.time() - t2} seconds.")

    # --- del temp output ---
    delete_files(temp_out_file)
    print("Delete temporary output files finished!")

    print("---processing finished---")
