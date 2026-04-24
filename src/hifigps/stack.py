"""
Functions used to stack galaxy-pairs signal.

"""

from typing import Sequence, Tuple

import numpy as np


def get_projection_lengths(
    p1: Sequence[int], p2: Sequence[int], shape: Tuple[int, int], scale: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    get horizontal and vertical projection lengths for each grid of data, using p2 -p1 as unit vector.

    Args:
        p1 (list[int, int]): the coodinates of point one.
        p2 (list[int, int]): the coodinates of point two.
        shape (tuple[int, int]): shape of the map.
        scale (bool): scale the result with the length of vector.
    Returns:
        len_h (np.ndarray):  horizontal lengths
        len_v (np.ndarray):  vertical lengths
    """
    # calculate unit vector
    p1_array = np.array(p1)
    p2_array = np.array(p2)
    m = np.linalg.norm(p2_array - p1_array)

    u = (p2_array - p1_array) / m  # the unit vector in the direction of p2 - p1
    v = np.array([-u[1], u[0]])  # the perpendicular vector to u

    # vetorize the grids
    s1, s2 = shape
    grids = np.zeros([s1, s2, 2])
    coor_1 = np.arange(0, s1)
    coor_2 = np.arange(0, s2)
    grids[:, :, 0] += coor_1[:, None]
    grids[:, :, 1] += coor_2[None, :]
    vec = grids - p1_array[None, :]

    # calculate the horizontal and vertical lengths of the projections
    len_h = np.dot(vec, u)
    len_v = np.dot(vec, v)

    # scale the lengths
    if scale:
        len_h = len_h / m * 2
        len_v = len_v / m * 2

    return len_h, len_v


def hist_data_2d(
    mask_data: np.ma.MaskedArray,
    h: np.ndarray,
    v: np.ndarray,
    histbins: Sequence[np.ndarray],
) -> np.ma.MaskedArray:
    """
    using histogram2d to rotate and rescale the 2d data.

    Args:
        mask_data (np.ma.masked_array): the 2d data with masked invalid data points.
        h (np.array): the horizontal projection lengths for each grid of data.
        v (np.array): the vertical projection lengths for each grid of data.
        histbins (sequence[np.ndarray]): the bins of histogram.

    Return:
        the averaged histogram result, the bins which no points fall in were masked.
    """

    mf = mask_data.mask.flatten()

    # genearate a array with 1 for valid data points and 0 for invalid data points to count the total valid points.
    valid = np.ones_like(mask_data.flatten())
    valid[mf] = 0

    # unmasked data values
    data = mask_data.compressed()

    # Compute the 2D histogram with the specified bins
    n = np.histogram2d(h, v, histbins, weights=valid)[
        0
    ]  # Count of the valid data points in each bin
    valid_sum = np.histogram2d(h[~mf], v[~mf], histbins, weights=data)[
        0
    ]  # Sum of valid data values in each bin

    # mask the zero counts
    masked_sum = np.ma.array(valid_sum, mask=n == 0.0)

    # Normalize the data values with respect to counts and rotate the result
    s = (masked_sum / n).T

    return s


def hist_data_3d(
    mask_data: np.ma.MaskedArray,
    p1: Sequence,  # pyright: ignore[reportMissingTypeArgument]
    p2: Sequence,  # pyright: ignore[reportMissingTypeArgument]
    hist_bins: Sequence[np.ndarray],
) -> np.ma.MaskedArray:
    """
    For each frequency slices, use histogram2d to rotate and rescale it.

    Args:
        mask_data (np.ma.array): the 3d data with masked invalid data points.
        p1 (list[int, int]): the indices of point one.
        p2 (list[int, int]): the indices of point two.
        hist_bins (sequence[np.ndarray]): the bins of histogram.
    Return:
        the averaged histogram results for each 3d data, the bins which no points fall in were masked.
    """
    if mask_data.ndim == 3:
        trip = True
        nf, nx, ny = mask_data.shape
    elif mask_data.ndim == 2:
        trip = False
        nx, ny = mask_data.shape
    else:
        raise ValueError("input mask_data must not be 1-dimensional")

    # get the horizontal and vertical projection lengths for each grid of data
    len_h, len_v = get_projection_lengths(p1, p2, (nx, ny))

    # Flatten the arrays
    h = len_h.flatten() - 1.0  # minus one to move two points to (-1,0) and (1,0)
    v = len_v.flatten()

    if trip:
        signals = np.ma.array(
            [hist_data_2d(mask_data[i, :, :], h, v, hist_bins) for i in range(nf)]  # pyright: ignore[reportPossiblyUnboundVariable]
        )
    else:
        signals = hist_data_2d(mask_data, h, v, hist_bins)

    return signals


def cut_freq(
    freqs: np.ndarray, flags: np.ndarray, is_pro_ra: bool, freq_min: int, freq_max: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    cut frequency slices outside the range, and flags

    Args:
        freqs (np.ndarray): the indices of frequency slices.
        flags (np.ndarray): the flags for each frequency slices.
        is_pro_ra (bool): if True, this pairs will be projected along R.A. axis.
        freq_min (int): the minimum frequency index.
        freq_max (int): the maximum frequency index.
    """
    valid = (freqs >= freq_min) & (freqs <= freq_max)

    if is_pro_ra:
        flags[~valid, :] = True
    else:
        flags[:, ~valid] = True

    _ = np.clip(freqs, freq_min, freq_max, out=freqs)

    return freqs, flags
