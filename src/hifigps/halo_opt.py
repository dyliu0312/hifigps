"""
Optimized verison of numerical halo-contribution fitting code.
"""

from typing import Callable, List, Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import NDArray
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares

from hifigps.bins import get_id_edge
from hifigps.utils import (
    get_coord,
    get_r_theta,
    info_fitness,
    yield_mask_data,
)


def halo_fit(
    data: NDArray,
    rbin_e: NDArray = np.arange(0, 6, 0.04),
    lc: List[float] = [-1, 0],
    rc: List[float] = [1, 0],
    mask: Optional[NDArray] = None,
    weight: Optional[NDArray] = None,
    coord: Optional[Sequence[NDArray]] = None,
    info_fit: bool = True,
    smooth_sigma: Optional[float] = None,
) -> Tuple[NDArray, NDArray]:
    """
    Advanced halo fitting function with weighted fitting and profile smoothing.

    This function models the data as a superposition of two radially symmetric
    profiles (halos) centered at `lc` and `rc`. It numerically fits the
    1D profile that best describes the 2D data and supports weighted least squares.

    Args:
        data (np.ndarray): 2D data to be fitted.
        rbin_e (np.ndarray): The edges of radius bins for the 1D halo profile.
        lc (list): The center of the left halo.
        rc (list): The center of the right halo.
        mask (np.ndarray, optional): A boolean mask to exclude regions from the
                                     fitting process. Defaults to None.
        weight (np.ndarray, optional): Weights for each data point. If None,
                                        equal weights are used. Defaults to None.
        coord (Sequence[np.ndarray], optional): Coordinates of the data. If None,
                                                they are generated based on the
                                                data shape. Defaults to None.
        info_fit (bool, optional): If True, print the fitting goodness-of-fit
                                   statistics. Defaults to True.
        smooth_sigma (Optional[float], optional): The sigma of the Gaussian
            filter to smooth the final fitted 2D halo map. If None, no
            smoothing is applied. Defaults to None.

    Returns:
        Tuple[np.ndarray, np.ndarray]:
            - fit_map: The 2D map of the fitted halo model.
            - fit_paras: The fitted 1D halo profile.
    """
    # Generate coordinates and weights if not provided
    if coord is None:
        coord = get_coord(data.shape[0])
    if weight is None:
        weight = np.ones_like(data)

    # Calculate radial distance of each pixel from the two halo centers
    r1 = get_r_theta(lc, coord)
    r2 = get_r_theta(rc, coord)

    # Assign each pixel to a radial bin index based on its distance
    r1_ids = get_id_edge(r1, rbin_e)  # pyright: ignore[reportArgumentType]
    r2_ids = get_id_edge(r2, rbin_e)  # pyright: ignore[reportArgumentType]

    # Apply the mask to the data, weights, and radial bin IDs
    # This flattens the arrays to 1D for the fitting process
    m_data, m_weight, m_r1_ids, m_r2_ids = yield_mask_data(
        data,
        weight,
        r1_ids,  # pyright: ignore[reportArgumentType]
        r2_ids,  # pyright: ignore[reportArgumentType]
        mask=mask,
        flat=True,  # pyright: ignore[reportArgumentType]
    )

    # --- Define helper functions for the optimization ---
    def get_fitmap(profile: NDArray, r1_ids, r2_ids) -> NDArray:
        """Constructs a 2D halo map from a 1D radial profile."""
        # Pad the profile to prevent index out-of-bounds errors
        padded_profile = np.pad(
            profile,
            (0, max(0, max(r1_ids.max(), r2_ids.max()) - len(profile) + 1)),
            "edge",
        )
        # The full model is the sum of two halos
        return padded_profile[r1_ids] + padded_profile[r2_ids]

    def error(
        paras: NDArray, data: NDArray, r1_ids, r2_ids, weight: NDArray
    ) -> NDArray:
        """Calculates the weighted residual for the least-squares fit."""
        fitmap = get_fitmap(paras, r1_ids, r2_ids)
        residual = data - fitmap
        return weight * residual

    # --- Perform the fitting ---
    num_bins = len(rbin_e) - 1
    # Initialize the profile with random values scaled by the data maximum
    init_paras = np.random.random(size=num_bins) * np.nanmax(m_data)

    # Run the non-linear least-squares optimization
    result = least_squares(
        error, init_paras, args=(m_data, m_r1_ids, m_r2_ids, m_weight)
    )
    fit_paras: NDArray = result.x

    # Optionally print goodness-of-fit statistics
    if info_fit:
        fit_map_masked = get_fitmap(fit_paras, m_r1_ids, m_r2_ids)
        info_fitness(m_data, fit_map_masked, num_bins)

    # Reconstruct the full 2D fitted map using the original (unmasked) bin IDs
    fit_map = get_fitmap(fit_paras, r1_ids, r2_ids)

    # Smooth the final 2D map if requested
    if smooth_sigma is not None:
        fit_map = gaussian_filter(fit_map, sigma=smooth_sigma)

    if isinstance(data, np.ma.MaskedArray):
        fit_map: NDArray = np.ma.array(fit_map, mask=data.mask)

    return fit_map, fit_paras


def halo_subtract(
    data: Union[NDArray, Sequence[NDArray]],
    func_halo_fit: Callable[..., Tuple[NDArray, NDArray]] = halo_fit,
    rbin_e: NDArray = np.arange(0, 6, 0.04),
    lc: List[float] = [-1, 0],
    rc: List[float] = [1, 0],
    mask: Optional[NDArray] = None,
    weight: Optional[NDArray] = None,
    coord: Optional[Sequence[NDArray]] = None,
    info_fit: bool = True,
    smooth_sigma: Optional[float] = None,
) -> Union[
    Tuple[NDArray, NDArray, NDArray],
    Tuple[List[NDArray], List[NDArray], NDArray],
]:
    """
    For each dataset provided, fit a halo profile and subtract it.

    This function acts as a wrapper around a halo fitting function (like `halo_fit`).
    It can process a single 2D array or a list of 2D arrays, making it
    convenient for batch processing.

    Args:
        data (Union[np.ndarray, Sequence[np.ndarray]]): A single 2D data array
            or a list of arrays to be fitted and subtracted.
        func_halo_fit (Callable): The function to use for fitting the halo.
            Defaults to `halo_fit`.
        rbin_e (np.ndarray): The edges of the radius bins for the halo profile.
        lc (list): The center of the left halo.
        rc (list): The center of the right halo.
        mask (np.ndarray, optional): Mask to exclude regions from fitting.
            Defaults to None.
        weight (np.ndarray, optional): Weights for the fitting. Defaults to None.
        coord (Sequence[np.ndarray], optional): Coordinate grid. Defaults to None.
        info_fit (bool, optional): If True, print fitting results. Defaults to True.
        smooth_sigma (Optional[float], optional): The sigma for Gaussian smoothing,
            passed to `func_halo_fit`. If None, no smoothing is applied.
            Defaults to None.

    Returns:
        Tuple containing the fitted map(s), the residual(s), and the mask.
        - If a single array was input: (fitted_map, residual, mask).
        - If a list of arrays was input: (list_of_fitted_maps, list_of_residuals, mask).
    """
    # Check if the input is a single np.ndarray or a list of them
    is_single = isinstance(data, np.ndarray)
    if is_single:
        # Wrap the single array in a list for consistent processing
        data = [data]

    fit: List[NDArray] = []
    res: List[NDArray] = []
    # Iterate through each dataset and perform halo fitting and subtraction
    for d in data:
        (
            fitmap,
            _,
        ) = func_halo_fit(
            d, rbin_e, lc, rc, weight, mask, coord, info_fit, smooth_sigma
        )

        fit.append(fitmap)
        res.append(d - fitmap)

    # Return the results in the same format as the input (single or list)
    if is_single:
        return fit[0], res[0], mask  # pyright: ignore[reportReturnType]
    else:
        return fit, res, mask  # pyright: ignore[reportReturnType]
