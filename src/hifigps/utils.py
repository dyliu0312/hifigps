"""
A collection of utility functions for fitting.

Including:

1.  Coordinates and Masking Utilities:
    - Generate 2D coordinate grids (`get_coord`).
    - Convert Cartesian coordinates to polar coordinates (`get_r_theta`).
    - Create various shapes of masks (square, circular sectors) for region
      selection (`get_mask_square`, `get_mask_sector`).

2.  Data Handling:
    - A generator to efficiently apply masks to multiple data arrays and
      optionally flatten them (`yield_mask_data`).

3.  Goodness of Fit:
    - A function to compute and display common goodness-of-fit statistics
      like chi-squared, reduced chi-squared, and RMS residual (`info_fitness`).
"""

from typing import Generator, Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import NDArray

MaskType = Union[NDArray[np.bool_], np.bool_, None]


# --- 1. Coordinates and Masking Utilities ---
def get_coord(
    n: int = 120,
    xlim: Sequence[float] = (-3, 3),
    ylim: Optional[Sequence[float]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate coordinate grid for 2D modeling.

    Parameters
    ----------
    n : int, optional
        Number of grid points in each dimension, by default 120
    xlim : Sequence[float], optional
        x-coordinate limits, by default (-3, 3)
    ylim : Sequence[float], optional
        y-coordinate limits, by default None (same as xlim)

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        x and y coordinate arrays
    """
    if ylim is None:
        ylim = xlim
    x = np.linspace(xlim[0], xlim[1], n)
    y = np.linspace(ylim[0], ylim[1], n)
    X, Y = np.meshgrid(x, y)
    return X, Y


def get_r_theta(
    center: Sequence[float] = [0, 0],
    coord: Optional[Sequence[np.ndarray]] = None,
    get_theta: bool = False,
) -> Union[np.ndarray, Sequence[np.ndarray]]:
    """
    Compute radial (and angular) coordinates from a 2D meshgrid.

    This function calculates the radial distance `r` and, optionally, the
    angle `theta` for each point on a coordinate grid with respect to a
    specified center.

    Parameters
    ----------
    center : Sequence[float], optional
        The [x, y] coordinates of the center, by default [0, 0].
    coord : Optional[Sequence[np.ndarray]], optional
        A tuple of (X, Y) coordinate arrays. If None, a default grid is
        generated using `get_coord()`. By default None.
    get_theta : bool, optional
        If True, return both `r` and `theta`. If False, return only `r`.
        By default False.

    Returns
    -------
    Union[np.ndarray, Sequence[np.ndarray]]
        - If `get_theta` is False: The radial coordinate array `r`.
        - If `get_theta` is True: A tuple of (r, theta) coordinate arrays.
    """
    if coord is None:
        coord = get_coord()
    x, y = coord
    x = x - center[0]
    y = y - center[1]
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    if get_theta:
        return r, theta
    else:
        return r


def get_mask_square(
    coord: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    xlim: Tuple[float, float] = (-1, 1),
    ylim: Tuple[float, float] = (-1, 1),
) -> np.ndarray:
    """
    Create a square boolean mask for region selection.

    This function generates a 2D boolean array where `True` values correspond
    to the region inside the specified rectangular limits.

    Parameters
    ----------
    coord : Optional[Tuple[np.ndarray, np.ndarray]], optional
        A tuple of (X, Y) coordinate arrays. If None, a default grid is
        generated using `get_coord()`. By default None.
    xlim : Tuple[float, float], optional
        The (min, max) limits for the x-coordinate, by default (-1, 1).
    ylim : Tuple[float, float], optional
        The (min, max) limits for the y-coordinate, by default (-1, 1).

    Returns
    -------
    np.ndarray
        A 2D boolean mask array.
    """
    if coord is None:
        coord = get_coord()
    x, y = coord
    x1, x2 = xlim
    y1, y2 = ylim
    mask = (x1 < x) & (x < x2) & (y1 < y) & (y < y2)
    return mask


def _get_mask_sector1(
    r: np.ndarray,
    theta: np.ndarray,
    rlim: Sequence[float] = (0.5, 1),
    thetalim: Sequence[float] = (-np.pi / 4, np.pi / 4),
) -> np.ndarray:
    """
    Create a single circular sector mask.

    This is an internal utility function that generates a boolean mask for a
    region defined by radial and angular limits.

    Parameters
    ----------
    r : np.ndarray
        The radial coordinate array.
    theta : np.ndarray
        The angular coordinate array (in radians).
    rlim : Sequence[float], optional
        The (min, max) limits for the radius, by default (0.5, 1).
    thetalim : Sequence[float], optional
        The (min, max) limits for the angle, by default (-np.pi / 4, np.pi / 4).

    Returns
    -------
    np.ndarray
        A 2D boolean mask for the specified sector.
    """
    mask = (rlim[0] < r) & (r < rlim[1]) & (thetalim[0] < theta) & (theta < thetalim[1])
    return mask


def _get_mask_sector2(
    r: np.ndarray,
    theta: np.ndarray,
    rlim: Sequence[float] = (0.5, 1),
    thetalim: Sequence[float] = (-np.pi / 4, np.pi / 4),
) -> np.ndarray:
    """
    Create a symmetric two-sector mask rotated by ±90 degrees.

    This internal utility function generates a mask composed of two identical
    sectors oriented perpendicular to each other.

    Parameters
    ----------
    r : np.ndarray
        The radial coordinate array.
    theta : np.ndarray
        The angular coordinate array (in radians).
    rlim : Sequence[float], optional
        The (min, max) limits for the radius of each sector, by default (0.5, 1).
    thetalim : Sequence[float], optional
        The (min, max) angular limits defining the opening of each sector,
        by default (-np.pi / 4, np.pi / 4).

    Returns
    -------
    np.ndarray
        A 2D boolean mask representing the two combined sectors.
    """

    # Create two symmetric sectors rotated by ±90 degrees
    mask1 = _get_mask_sector1(
        r,
        theta,
        rlim,
        (thetalim[0] - np.pi / 2, thetalim[1] - np.pi / 2),
    )
    mask2 = _get_mask_sector1(
        r,
        theta,
        rlim,
        (thetalim[0] + np.pi / 2, thetalim[1] + np.pi / 2),
    )
    mask = mask1 | mask2
    return mask


def _get_mask_sector4(
    r1: np.ndarray,
    theta1: np.ndarray,
    r2: np.ndarray,
    theta2: np.ndarray,
    rlim: Sequence[float] = (0.5, 1),
    thetalim: Sequence[float] = (-np.pi / 4, np.pi / 4),
) -> np.ndarray:
    """
    Create a four-sector mask from two pairs of sectors with offset centers.

    This internal utility combines two two-sector masks (generated by
    `_get_mask_sector2`) centered at different locations.

    Parameters
    ----------
    r1 : np.ndarray
        The radial coordinate array relative to the first center.
    theta1 : np.ndarray
        The angular coordinate array relative to the first center.
    r2 : np.ndarray
        The radial coordinate array relative to the second center.
    theta2 : np.ndarray
        The angular coordinate array relative to the second center.
    rlim : Sequence[float], optional
        The (min, max) radial limits for all sectors, by default (0.5, 1).
    thetalim : Sequence[float], optional
        The (min, max) angular limits for all sectors, by default (-np.pi / 4, np.pi / 4).

    Returns
    -------
    np.ndarray
        A 2D boolean mask representing the four combined sectors.
    """
    mask1 = _get_mask_sector2(r1, theta1, rlim, thetalim)
    mask2 = _get_mask_sector2(r2, theta2, rlim, thetalim)
    return mask1 | mask2


def get_mask_sector(
    rlim: Sequence[float] = (0.2, 5),
    thetalim: Sequence[float] = (-np.pi / 4, np.pi / 4),
    rc: Sequence[float] = (1, 0),
    lc: Sequence[float] = (-1, 0),
    coord: Optional[Sequence[np.ndarray]] = None,
) -> np.ndarray:
    """
    Create a complex four-sector mask with two offset centers.

    This function generates a mask composed of four sectors. Two sectors are
    centered at `rc` and two at `lc`. This shape is often used to select
    outflow regions in astrophysical maps.

    Parameters
    ----------
    rlim : Sequence[float], optional
        The (min, max) radial limits for the sectors, by default (0.2, 5).
    thetalim : Sequence[float], optional
        The (min, max) angular limits for the sectors, by default (-np.pi / 4, np.pi / 4).
        To get a full circle, use (-np.pi/2, np.pi/2).
    rc : Sequence[float], optional
        The [x, y] center for the first pair of sectors, by default (1, 0).
    lc : Sequence[float], optional
        The [x, y] center for the second pair of sectors, by default (-1, 0).
    coord : Optional[Sequence[np.ndarray]], optional
        A tuple of (X, Y) coordinate arrays. If None, a default grid is
        generated using `get_coord()`. By default None.

    Returns
    -------
    np.ndarray
        A 2D boolean mask array representing the four combined sectors.
    """
    if coord is None:
        coord = get_coord()

    r1, theta1 = get_r_theta(rc, coord, get_theta=True)
    r2, theta2 = get_r_theta(lc, coord, get_theta=True)

    return _get_mask_sector4(r1, theta1, r2, theta2, rlim, thetalim)


# --- 2. Yield masked datasets ---
def yield_mask_data(
    *args: np.ndarray,
    mask: MaskType = None,
    flat: bool = False,
) -> Generator[np.ndarray]:
    """
    Apply a mask to one or more arrays and yield the results.

    This function provides a generator that efficiently applies a boolean mask
    to a sequence of arrays. It can optionally flatten the masked output.

    Parameters
    ----------
    *args : np.ndarray
        One or more NumPy arrays to be masked.
    mask : MaskType, optional
        A boolean array for masking. If a mask is provided, only elements
        corresponding to `True` values are yielded. If `None`, the entire
        array is yielded. Defaults to None.
    flat : bool, optional
        If True, the masked output arrays will be flattened to 1D.
        By default False.

    Yields
    ------
    Generator[np.ndarray]
        A generator that yields each masked (and optionally flattened) array.
    """
    for arg in args:
        if mask is not None:
            if flat:
                yield arg[mask].flatten()
            else:
                yield arg[mask]
        else:
            if flat:
                yield arg.flatten()
            else:
                yield arg


def info_fitness(
    data: np.ndarray,
    fit_data: np.ndarray,
    n_params: int,
    return_res: bool = False,
) -> Union[np.ndarray, None]:
    """
    Calculate and print goodness-of-fit statistics.

    This function computes and prints key metrics for evaluating the quality
    of a model fit, including chi-squared, reduced chi-squared, and the RMS
    of the residuals.

    Parameters
    ----------
    data : np.ndarray
        The original data array.
    fit_data : np.ndarray
        The model-fitted data array.
    n_params : int
        The number of free parameters in the model.
    return_res : bool, optional
        If True, the function returns the residuals array. If False, it
        returns None. By default False.

    Returns
    -------
    Union[np.ndarray, None]
        The array of residuals (data - fit_data) if `return_res` is True,
        otherwise None.
    """
    # Calculate goodness of fit
    residuals = data - fit_data
    chi2 = np.sum(residuals**2)
    dof = data.size - n_params
    reduced_chi2 = chi2 / dof if dof > 0 else np.inf

    print("\nGoodness of fit:")
    print(f"Chi-squared: {chi2:.3f}")
    print(f"Reduced chi-squared: {reduced_chi2:.3f}")
    print(f"RMS residual: {np.sqrt(np.mean(residuals**2)):.3f}\n")

    if return_res:
        return residuals
