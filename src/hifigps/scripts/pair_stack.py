"""
main script to finish galaxy pairwise stacking process.
"""

import argparse
import gc
import multiprocessing as mp
import os
import sys
from typing import List, Sequence  # pyright: ignore[reportDeprecated]

import h5py as h5
import numpy as np
from hifigps.bins import (  # type: ignore
    get_ids_edge,
    set_resbins,
)
from hifigps.data import (  # type: ignore
    is_exist,
    read_h5,
    save_h5,
    split_data_generator,
)
from hifigps.stack import (  # type: ignore
    cut_freq,
    hist_data_3d,
)
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Galaxy pairwise stacking script.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Example:
  hifigps-stack --map-base /data/ --map-prefix map \\
                --paircat-base /data/ --paircat-prefix paircat \\
                --out-base /output/ --out-prefix stack \\
                --nfs 10 --ssize 1000 --nworker 24
""",
    )

    # Required arguments
    parser.add_argument("--map-base", required=True, help="Base path to input map")
    parser.add_argument("--map-prefix", required=True, help="Prefix of input map file")
    parser.add_argument(
        "--paircat-base", required=True, help="Base path to pair catalog"
    )
    parser.add_argument(
        "--paircat-prefix", required=True, help="Prefix of pair catalog file"
    )
    parser.add_argument("--out-base", required=True, help="Base path to output")
    parser.add_argument("--out-prefix", required=True, help="Prefix of output file")
    parser.add_argument(
        "--nfs", type=int, required=True, help="Number of frequency slices to stack"
    )
    parser.add_argument(
        "--ssize", type=int, required=True, help="Split size for processing"
    )

    # Optional arguments
    parser.add_argument(
        "--map-masked",
        type=lambda x: x.lower() == "true",
        default=True,
        help="Whether input map is masked (default: True)",
    )
    parser.add_argument(
        "--map-keys",
        default="T,mask,f_bin_edge,x_bin_edge,y_bin_edge",
        help="Comma-separated keys for map data (default: T,mask,f_bin_edge,x_bin_edge,y_bin_edge)",
    )
    parser.add_argument(
        "--paircat-keys",
        default="is_ra,pos",
        help="Comma-separated keys for pair catalog (default: is_ra,pos)",
    )
    parser.add_argument(
        "--out-keys",
        default="Signal,Mask",
        help="Comma-separated keys for output (default: Signal,Mask)",
    )
    parser.add_argument(
        "--nworker",
        type=int,
        default=None,
        help="Number of workers (default: CPU count)",
    )
    parser.add_argument(
        "--random-flip",
        type=lambda x: x.lower() == "true",
        default=True,
        help="Random flip signal (default: True)",
    )
    parser.add_argument(
        "--halfwidth", type=float, default=3.0, help="Stack result map half-width (default: 3.0)"
    )
    parser.add_argument(
        "--npix-x", type=int, default=120, help="Stack result map X pixels (default: 120)"
    )
    parser.add_argument(
        "--npix-y", type=int, default=120, help="Stack result map Y pixels (default: 120)"
    )
    parser.add_argument(
        "--skip-exist",
        type=lambda x: x.lower() == "true",
        default=False,
        help="Skip existing splits (default: False)",
    )
    parser.add_argument(
        "--compression",
        default="gzip",
        choices=["gzip", "lzf", "none"],
        help="Compression method (default: gzip)",
    )

    return parser.parse_args()


def join_path(base: str, prefix: str, extension: str = ".h5") -> str:
    # Ensure extension starts with a dot
    if not extension.startswith("."):
        extension = "." + extension
    path = os.path.join(base, prefix)
    if not path.endswith(extension):
        return path + extension
    else:
        return path


def get_ipm(pair_catalog, is_pro_ra, map_bins, nfreqslice, hist_bins):
    """
    get indivual pair's map.

    params:
        paircat (list): the paircat of two sources.
        pro_ra (bool): whether to project on ra axis.
        mapbins (tuple): the map bins.
        nfs (int): the number of frequency slices to extract.
        histbins (tuple): the histogram bins.

    return: the indivual pair's map.
    """

    # get map index
    ira1, idec1, ifreq1 = get_ids_edge(pair_catalog[0], map_bins)
    ira2, idec2, ifreq2 = get_ids_edge(pair_catalog[1], map_bins)

    # dimensions
    nfreq, nra, ndec = (
        mask_map.shape
    )  # ignore the error of using before assinment, this is necessray for multiprocessing

    # Extracting individual pair cut
    if nfreqslice == 0:
        signal = np.ma.zeros([nra, ndec])
        flag = np.zeros([nra, ndec], dtype=bool)

        if ifreq1 == ifreq2:
            signal = mask_map[ifreq1]
        else:
            if is_pro_ra:
                ifids = ifreq1 + float(ifreq2 - ifreq1) / float(ira2 - ira1) * (
                    np.arange(nra) - ira1
                )
            else:
                ifids = ifreq1 + float(ifreq2 - ifreq1) / float(idec2 - idec1) * (
                    np.arange(ndec) - idec1
                )

            freq_indices = np.round(ifids).astype(int)
            cut_freq_indices, flag = cut_freq(
                freq_indices, flag, is_pro_ra, 0, nfreq - 1
            )
            signal = (
                mask_map[cut_freq_indices, range(nra)]
                if is_pro_ra
                else mask_map[cut_freq_indices, :, range(nra)]
            )

    else:
        total_fs = 2 * nfreqslice + 1
        signal = np.ma.zeros([total_fs, nra, ndec])
        flag = np.zeros([total_fs, nra, ndec], dtype=bool)

        if ifreq1 == ifreq2:
            freq_indices = np.arange(ifreq1 - nfreqslice, ifreq1 + nfreqslice + 1)
            vallid = (freq_indices >= 0) & (freq_indices <= nfreq - 1)
            if vallid.sum() != total_fs:
                np.clip(freq_indices, 0, nfreq - 1, out=freq_indices)
                flag[~vallid] = True

            signal = mask_map[freq_indices]

        else:
            if is_pro_ra:
                ifids = ifreq1 + float(ifreq2 - ifreq1) / float(ira2 - ira1) * (
                    np.arange(nra) - ira1
                )
            else:
                ifids = ifreq1 + float(ifreq2 - ifreq1) / float(idec2 - idec1) * (
                    np.arange(ndec) - idec1
                )

            freq_indices = np.round(ifids).astype(int)
            for i in range(total_fs):
                cut_freq_indices, flag[i] = cut_freq(
                    freq_indices - nfreqslice + i, flag[i], is_pro_ra, 0, nfreq - 1
                )
                signal[i] = (
                    mask_map[cut_freq_indices, range(nra)]
                    if is_pro_ra
                    else mask_map[cut_freq_indices, :, range(nra)]
                )

    # add up the flag mask
    signal.mask += flag

    # hist data
    p1 = [ira1, idec1]
    p2 = [ira2, idec2]

    s = hist_data_3d(signal, p1, p2, hist_bins)

    return s


def stack_mp(
    nworker, pair_catalog, is_pro_ra, map_bins, nfreqslice, hist_bins, random_flip=True
):
    """
    stacking in multiprocessing async pool
    """
    # multiprocess
    total_pairs = len(pair_catalog)

    p = mp.Pool(nworker)
    res = [
        p.apply_async(
            get_ipm, (pair_catalog[i], is_pro_ra[i], map_bins, nfreqslice, hist_bins)
        )
        for i in range(total_pairs)
    ]

    p.close()
    p.join()

    # get the result
    if random_flip:
        if nfreqslice == 0:
            signal = [
                np.flip(i.get(), axis=0) if np.random.choice([True, False]) else i.get()
                for i in res
            ]  # upside down
            signal = [
                np.flip(s, axis=1) if np.random.choice([True, False]) else s
                for s in signal
            ]  # left right
        else:
            signal = [
                np.flip(i.get(), axis=1) if np.random.choice([True, False]) else i.get()
                for i in res
            ]  # upside down
            signal = [
                np.flip(s, axis=2) if np.random.choice([True, False]) else s
                for s in signal
            ]  # left right
    else:
        signal = [i.get() for i in res]

    # average those result.
    s = np.ma.array(signal).mean(axis=0)
    return s


def stack_run(
    output: str,
    pair_catalog: np.ndarray,
    is_pro_ra: np.ndarray,
    map_bins: Sequence[np.ndarray],
    hist_bins: Sequence[np.ndarray],
    nfreqslice: int,
    split_size: int = 500,
    nworker: int = 24,
    random_flip: bool = True,
    savekeys: List[str] = ["Signal", "Mask"],
    compression: str = "gzip",
    skip_exist: bool = True,
):
    """
    split data to stack and save seprately.

    Parameters
    ----------
    output : str
        output file name
    pair_catalog : np.ndarray
        pair catalog
    is_pro_ra : np.ndarray
        boolean array, indicating if the pair stacking is gona projection along ra (i.e. ra longer than dec)
    map_bins : list or tuple
        list of center bins of the map (ra, dec, freq)
    hist_bins : list or tuple
        list of center bins of the histogram (x, y)
    nfreqslice : int
        number of frequency slice to extract during stack
    split_size : int, optional
        number of pairs to process in one split, by default 500
    nworker : int, optional
        number of workers to use, by default 24
    random_flip : bool, optional
        if random flip the signal of pairs, by default True
    savekeys : list, optional
        list of keys of output dataset to save, by default ["Signal", "Mask"]
    compression : str, optional
        compression method of h5 file, by default "gzip"
    skip_exist : bool, optional
        if skip the split that already exist in the output file, by default True

    Returns
    -------
    None
    """

    # splite data for multiprocessing

    pbar = tqdm(
        enumerate(split_data_generator(split_size, pair_catalog, is_pro_ra)),
        total=len(pair_catalog) // split_size + 1,
        desc="Processing splits",
    )

    if skip_exist:
        print("--- skip_exist enabled, checking existing results ---")
        if is_exist(output):
            with h5.File(output, "r") as f:
                exist_keys = [i for i in f.keys()]
            if exist_keys == []:
                skip_exist = False
                print("--- No existing result found, disable skipping ---")
            else:
                print("--- Found existing result in output file, enable skipping ---")
        else:
            skip_exist = False
            print("--- Output file not found, disable skipping ---")

    for i, (ipair, ipra) in pbar:
        pbar.set_postfix({"split": f"{i}"})
        groupname = str(split_size) + "_" + str(i)
        if skip_exist and groupname in exist_keys:  # pyright: ignore[reportPossiblyUnboundVariable]
            continue

        s = stack_mp(nworker, ipair, ipra, map_bins, nfreqslice, hist_bins, random_flip)

        # save result
        save_h5(output, savekeys, [s.data, s.mask], groupname, compression=compression)


def main():

    args = parse_args()

    INPUT_MAP_BASE = args.map_base
    INPUT_MAP_PREFIX = args.map_prefix
    INPUT_PAIECAT_BASE = args.paircat_base
    INPUT_PAIRCAT_PREFIX = args.paircat_prefix
    OUTPUT_STACK_BASE = args.out_base
    OUTPUT_STACK_PREFIX = args.out_prefix
    OUTPUT_STACK_DATA_KEYS = args.out_keys.split(",")
    NFS = args.nfs
    SSIZE = args.ssize
    INPUT_MAP_MASKED = args.map_masked
    INPUT_MAP_KEYS = args.map_keys.split(",")
    INPUT_PAIECAT_KEYS = args.paircat_keys.split(",")
    nworker = args.nworker if args.nworker else mp.cpu_count()
    RANDOM_FLIP = args.random_flip
    HALFWIDTH = args.halfwidth
    NPIX_X = args.npix_x
    NPIX_Y = args.npix_y
    SKIP_EXIST = args.skip_exist
    COMPRESSION = args.compression if args.compression != "none" else None

    print("---- Starting Stacking Script ----")
    print("--- Initializing ---")
    # --- Print loaded configuration ---
    print("Processing with the following configuration:")
    print("--------------")
    print(f" Map: {join_path(INPUT_MAP_BASE, INPUT_MAP_PREFIX)}")
    print(f" Map Masked: {INPUT_MAP_MASKED}")
    print(f" Pair Catalog: {join_path(INPUT_PAIECAT_BASE, INPUT_PAIRCAT_PREFIX)}")
    print(f" Output Stack: {join_path(OUTPUT_STACK_BASE, OUTPUT_STACK_PREFIX)}")
    print(f" Output Stack Data Keys: {OUTPUT_STACK_DATA_KEYS}")
    print(f" NFS (frequency slices to stack): {NFS}")
    print(f" NWORKER (number of workers): {nworker}")
    print(f" SSIZE (catalog split size for processing): {SSIZE}")
    print(f" RANDOM_FLIP (randomly flip individual pair map): {RANDOM_FLIP}")
    print(f" COMPRESSION (compression method): {COMPRESSION}")
    print(f" SKIP_EXIST (skip existing stacks): {SKIP_EXIST}")
    print(f" HALFWIDTH (stack result map half-width): {HALFWIDTH}")
    print(f" NPIX_X (stack result map X pixels): {NPIX_X}")
    print(f" NPIX_Y (stack result map Y pixels): {NPIX_Y}")
    print(
        f" Stacked result map size: {2 * HALFWIDTH}x{2 * HALFWIDTH} with {NPIX_X}x{NPIX_Y} pixels"
    )
    print("--------------")

    # --- Construct tile-specific file paths ---

    map_path = join_path(INPUT_MAP_BASE, INPUT_MAP_PREFIX)
    paircat_path = join_path(INPUT_PAIECAT_BASE, INPUT_PAIRCAT_PREFIX)
    output_path = join_path(OUTPUT_STACK_BASE, OUTPUT_STACK_PREFIX)

    # --- Validate Input Files ---
    if not is_exist(map_path):
        print(f"Error: Map file not found at '{map_path}'. Exiting.", file=sys.stderr)
        sys.exit(1)

    if not is_exist(paircat_path):
        print(
            f"Error: Pair catalog not found at '{paircat_path}'. Exiting.",
            file=sys.stderr,
        )
        sys.exit(1)

    # --- Validate Output Path (Check if already exists, prevent overwrite) ---
    if is_exist(output_path):
        if SKIP_EXIST:
            print(
                "┌" + "-" * 62 + "┐",
                "│" + " " * 27 + "WARNING" + " " * 28 + "│",
                "├" + "-" * 62 + "┤",
                "| This job is currently trying to add new stacking results to  |",
                "| the output file that already exists.                         |",
                "└" + "-" * 62 + "┘",
                sep="\n",
            )
        else:
            print(
                "┌" + "-" * 62 + "┐",
                "│" + " " * 28 + "ERROR" + " " * 29 + "│",
                "├" + "-" * 62 + "┤",
                "| Output file already exists. To prevent accidental overwrite, |",
                "| this operation has been stopped.                             |",
                "|" + " " * 62 + "|",
                "| If you want to add new results to the existing file, instead |",
                "| of creating a new one, please set SKIP_EXIST = True.         |",
                "└" + "-" * 62 + "┘",
                sep="\n",
                file=sys.stderr,
            )
            sys.exit(1)

    # --- Print loaded configuration ---
    print("Processing with the following configuration:")
    print("--------------")
    print(f" Map: {map_path}")
    print(f" Map Masked: {INPUT_MAP_MASKED}")
    print(f" Pair Catalog: {paircat_path}")
    print(f" Output Stack: {output_path}")
    print(f" Output Stack Data Keys: {OUTPUT_STACK_DATA_KEYS}")
    print(f" NFS (frequency slices to stack): {NFS}")
    print(f" NWORKER (number of workers): {nworker}")
    print(f" SSIZE (catalog split size for processing): {SSIZE}")
    print(f" RANDOM_FLIP (randomly flip individual pair map): {RANDOM_FLIP}")
    print(f" COMPRESSION (compression method): {COMPRESSION}")
    print(f" SKIP_EXIST (skip existing stacks): {SKIP_EXIST}")
    print(f" HALFWIDTH (stack result map half-width): {HALFWIDTH}")
    print(f" NPIX_X (stack result map X pixels): {NPIX_X}")
    print(f" NPIX_Y (stack result map Y pixels): {NPIX_Y}")
    print(
        f" Stacked result map size: {2 * HALFWIDTH}x{2 * HALFWIDTH} with {NPIX_X}x{NPIX_Y} pixels"
    )
    print("--------------")

    # --- Loading Data ---
    print("\n--- Loading Data ---")

    ## paircat prepare
    # The pair catalog from the previous step saves 'is_ra' and 'pos'.
    # Ensure read_h5 correctly returns these.
    try:
        is_pra, pos = read_h5(paircat_path, INPUT_PAIECAT_KEYS)
    except Exception as e:
        print(f"Failed to load pair catalog from {paircat_path}: {e}", file=sys.stderr)
        sys.exit(1)

    # Shuffle the data to get more randomization
    # Concatenate is_pra (boolean) and pos_data (float) for shuffling.
    # We need to ensure 'is_pra' is numeric for column_stack, then cast back.
    pro_pos = np.column_stack(
        [is_pra.astype(np.int8), pos]
    )  # Convert bool to int8 for shuffling
    np.random.shuffle(pro_pos)

    is_pra = pro_pos[:, 0].astype(bool)  # Convert back to bool
    paircat_raw = pro_pos[:, 1:]  # Remaining columns are the position data

    # Split paircat to match the function input format: [[gal1_info], [gal2_info]]
    # Each inner list has [ra, dec, frequency]
    # 'pos' from previous script is [ra1, dec1, freq1, ra2, dec2, freq2]
    # So we need to reshape it into (N, 2, 3) where N is number of pairs.
    paircat = paircat_raw.reshape(
        -1, 2, 3
    )  # Reshape to N rows, 2 elements (galaxy), 3 values (ra, dec, freq)

    print(f"Loaded {len(is_pra)} galaxy pairs.")
    print(f"Pair catalog reshaped to: {paircat.shape}")

    ## mapfile prepare
    try:
        key_map, key_mask, key_fbin, key_xbin, key_ybin = INPUT_MAP_KEYS
        # load mapbins
        mapbin_keys_choice = [key_xbin, key_ybin, key_fbin]
        mapbins = read_h5(map_path, mapbin_keys_choice)
        # load map
        global mask_map # Make it global for multiprocessing access in get_ipm
        if INPUT_MAP_MASKED:
            map_array, mask_array = read_h5(map_path, [key_map, key_mask])
            mask_map = np.ma.masked_array(map_array, mask=mask_array, dtype=np.float32)
        else:
            map_array = read_h5(map_path, key_map)
            mask_map = np.ma.masked_array(
                map_array, mask=map_array == 0, dtype=np.float32
            )
            print("Masked map with zeros.")
        print(f"Loaded map from {map_path} with shape {map_array.shape}")  # pyright: ignore[reportAttributeAccessIssue]
    except Exception as e:
        print(f"Failed to load map data from {map_path}: {e}", file=sys.stderr)
        sys.exit(1)

    # Calculate resbins based on default args or env vars
    resbins = set_resbins(HALFWIDTH, NPIX_X, NPIX_Y, 2)

    # --- Stack ---
    print("\n---- Processing Stacking ----")
    stack_run(
        output=output_path,
        pair_catalog=paircat,
        is_pro_ra=is_pra,
        map_bins=mapbins,  # pyright: ignore[reportArgumentType]
        hist_bins=resbins,
        nfreqslice=NFS,
        split_size=SSIZE,
        nworker=nworker,
        random_flip=RANDOM_FLIP,
        savekeys=OUTPUT_STACK_DATA_KEYS,
        compression=COMPRESSION,
        skip_exist=SKIP_EXIST,
    )

    gc.collect()  # Trigger garbage collection
    print("\n---- Done ----")


if __name__ == "__main__":
    main()
