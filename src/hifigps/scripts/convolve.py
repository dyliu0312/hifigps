"""
Simulate the beam smearing effect.
"""

import argparse
import gc
import multiprocessing as mp
import time

from astropy.convolution import Gaussian2DKernel, convolve_fft
from hifigps.calculation import get_beam_npix
from hifigps.data import delete_files, get_filename, is_exist, np, read_h5, save_h5


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convolve a map with FAST main beam.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Example:
  hifigps-convolve /data/test.hdf5 /output/ 4
""",
    )
    parser.add_argument("map_file", help="Path to the input map file (.hdf5)")
    parser.add_argument("out_path", help="Path to the output directory")
    parser.add_argument("nworker", type=int, help="Number of worker processes")
    parser.add_argument(
        "-z", "--redshift", type=float, default=0.1, help="Redshift (default: 0.1)"
    )
    parser.add_argument(
        "-k",
        "--key",
        default="T",
        help="Key of map data in HDF5 file (default: T)",
    )
    parser.add_argument(
        "-d",
        "--dtype",
        default="float32",
        choices=["float32", "float64"],
        help="Data type of output (default: float32)",
    )
    return parser.parse_args()


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


def main():
    args = parse_args()

    Z = args.redshift
    KEY = args.key
    DTYPE = np.float32 if args.dtype == "float32" else np.float64
    TEMP_PREFIX = "convolved_temp_out_"
    OUT_SUFFIX = "_convolvedFASTbeam.hdf5"
    WAIT_TIME = 60

    nworker = args.nworker
    map_path = args.map_file
    out_path = args.out_path

    print(f"load data from {map_path}")
    print(f"save result at {out_path}")
    print(f"set {nworker} workers")

    print("---processing start----")
    beam_kernel = beam_kernel_fast(Z)
    filename = get_filename(map_path)
    out_name = out_path + filename + OUT_SUFFIX
    if is_exist(out_name):
        print(f"This map is already processed, the output file is {out_name}!")
        return

    t0 = time.time()
    map_data = read_h5(map_path, KEY)
    print(f"Loading the MAP data successfully with {time.time() - t0} seconds.")

    t1 = time.time()
    temp_out_file = []
    with mp.Pool(processes=nworker) as pool:
        for i, map_slice in enumerate(map_data):
            temp_out = out_path + filename + "_" + TEMP_PREFIX + f"{i}.hdf5"
            temp_out_file.append(temp_out)
            pool.apply(convolve_save, (map_slice, beam_kernel, temp_out, KEY, DTYPE))

    print(f"Convolution process finished with {time.time() - t1} seconds.")

    while True:
        if all([is_exist(fn) for fn in temp_out_file]):
            break
        time.sleep(WAIT_TIME)
        print("Waiting for the temporary output files to be saved")

    del map_data
    gc.collect()

    t2 = time.time()
    reconstruct_map(temp_out_file, out_name)
    print(f"Save convolved results finished with {time.time() - t2} seconds.")

    delete_files(temp_out_file)
    print("Delete temporary output files finished!")

    print("---processing finished---")


if __name__ == "__main__":
    main()
