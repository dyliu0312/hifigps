"""
Functions for bins
"""

from typing import List, Sequence, Tuple, Union  # pyright: ignore[reportDeprecated]

import numpy as np


def edge2center(edgebins: np.ndarray) -> np.ndarray:
    """
    turn edge bins to center bins.
    """
    return edgebins[:-1] + 0.5 * (edgebins[1:] - edgebins[:-1])


def get_linbin(
    st: float, ed: float, num: int, offset: int = 1, center: bool = True
) -> np.ndarray:
    """
    Get linear center/edge bins.

    Args:
        st (float): The starting value of the bins.
        ed (float): The ending value of the bins.
        num (int): The number of bins.
        offset (int, optional): The offset value for the bins. Defaults to 1.
        center (bool, optional): Determines if the bins are centered. Defaults to True.

    Returns:
        numpy.ndarray: An array representing the linear bins.
    """
    lbin = np.linspace(st, ed, num + offset)

    if center:
        lbin = edge2center(lbin)

    return lbin


def set_resbins(
    width: float, nx: int, ny: int, offset: int = 1, center: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    set stack result map bins.

    Args:
        width (float): The width of the bins.
        nx (int): The number of bins in the horizontal direction.
        ny (int): The number of bins in the vertical direction.
        offset (int, optional): The offset value for the bins. Defaults to 1.
        center (bool, optional): Determines if the bins are centered. Defaults to True.

    Returns:
        tuple: A tuple containing the horizontal and vertical bins.
    """

    # set horizental and vertical bins
    hbin = get_linbin(-width, width, nx, offset=offset, center=center)
    vbin = get_linbin(-width, width, ny, offset=offset, center=center)

    return hbin, vbin


def get_id_edge(
    coordinate: Union[float, Sequence[float], np.ndarray],
    bins: Union[List[float], np.ndarray],
    right: bool = False,
) -> Union[int, np.ndarray]:
    """
    Get the bin indices for given coordinates based on specified bin edges.

    Uses numpy.digitize() to find the indices of the bins that the coordinates fall into.
    The returned indices are adjusted to be zero-based.

    Note:
        - Returned indices less than 0 or greater than the size of bins are invalid
        - The first valid index is 0 after adjustment
        - The `right` parameter controls whether intervals include the right or left bin edge

    Args:
        coordinate: Input coordinate(s). Can be a single float or sequence of floats.
        bins: Array-like of bin edges. Must be monotonically increasing.
        right: If True, the bin intervals are [a, b) (right edge excluded).
               If False, the bin intervals are (a, b] (left edge excluded).
               Defaults to False.

    Returns:
        Bin indices as integer or numpy array. Zero-based indices indicating which bin
        each coordinate falls into. Invalid coordinates return indices outside the
        valid range [0, len(bins)-1].

    Example:
        >>> get_id_edge(2.5, [1, 2, 3, 4])
        np.int64(1)
        >>> get_id_edge([1.5, 2.5, 3.5], [1, 2, 3, 4])
        array([0, 1, 2])
    """
    inds = np.digitize(coordinate, bins=bins, right=right) - 1
    return inds


def get_ids_edge(
    coordinates: Tuple,  # pyright: ignore[reportMissingTypeArgument]
    bins: Tuple,  # pyright: ignore[reportMissingTypeArgument]
    right: bool = False,
) -> Tuple[Union[int, np.ndarray], Union[int, np.ndarray], Union[int, np.ndarray]]:
    """
    Get bin indices for 3D coordinates in X, Y, and frequency dimensions.

    Computes the bin indices for each coordinate component (X, Y, frequency) using
    their respective bin definitions.

    Args:
        coordinates: Tuple of three sequences containing X, Y, and frequency coordinates.
        bins: Tuple of three lists defining bin edges for X, Y, and frequency dimensions.
              Each must be monotonically increasing.
        right: If True, the bin intervals are [a, b) (right edge excluded).
               If False, the bin intervals are (a, b] (left edge excluded).
               Defaults to False.

    Returns:
        Tuple of three arrays (x_id, y_id, f_id) containing the bin indices for
        X, Y, and frequency coordinates respectively. All indices are zero-based.

    Example:
        >>> coords = ([1.5, 2.5], [2.5, 3.5], [0.5, 1.5])
        >>> bins = ([1, 2, 3], [2, 3, 4], [0, 1, 2])
        >>> get_ids_edge(coords, bins)
        (array([0, 1]), array([0, 1]), array([0, 1]))
    """
    x, y, f = coordinates
    xbin, ybin, fbin = bins

    x_id = get_id_edge(x, xbin, right)
    y_id = get_id_edge(y, ybin, right)
    f_id = get_id_edge(f, fbin, right)

    return x_id, y_id, f_id
