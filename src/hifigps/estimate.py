"""
estimate the signal level
"""

from typing import Optional, Sequence

import numpy as np
from astropy.modeling import fitting, models
from astropy.stats import sigma_clip

from hifigps.plot import plot_line_diff
from hifigps.utils import info_fitness, yield_mask_data


def get_int(x):
    """
    Round a float or array of floats to the nearest integer.

    Parameters
    ----------
    x : Union[float, np.ndarray]
        The input number or array.

    Returns
    -------
    Union[int, np.ndarray]
        The rounded integer or array of integers.
    """
    x = np.array(x)
    return np.rint(x).astype(int)


def get_start_end_pixel_index(
    x_min: float = -0.5,
    x_max: float = 0.5,
    shape: int = 120,
    lim: Sequence[float] = (-3, 3),
    offset: int = 1,
):
    """
    Get the start and end pixel indices for a given coordinate range.

    Parameters
    ----------
    x_min : float, optional
        The minimum world coordinate. Defaults to -0.5.
    x_max : float, optional
        The maximum world coordinate. Defaults to 0.5.
    shape : int, optional
        The number of pixels along the axis. Defaults to 120.
    lim : Sequence[float], optional
        The world coordinate limits of the axis, as (min, max).
        Defaults to (-3, 3).
    offset : int, optional
        An offset to add to the calculated end index. Defaults to 1.

    Returns
    -------
    Tuple[int, int]
        The calculated start and end pixel indices.
    """
    xx = np.linspace(lim[0], lim[1], shape)

    x_start = np.argmin(np.abs(xx - x_min))
    x_end = np.argmin(np.abs(xx - x_max)) + offset

    return x_start, x_end


def get_gaussian_fit(
    x: np.ndarray,
    y: np.ndarray,
    weights: Optional[np.ndarray] = None,
    sigma_clip_sigma: Optional[float] = None,
    bounds: Optional[dict] = None,
    print_info: bool = True,
    **kwargs,
):
    """
    Fit a Gaussian + constant model to 1D data with optional sigma clipping.

    Parameters
    ----------
    x : np.ndarray
        The independent variable.
    y : np.ndarray
        The dependent variable.
    weights : Optional[np.ndarray], optional
        Weights for the fitting. Defaults to None.
    sigma_clip_sigma : Optional[float], optional
        The number of standard deviations for sigma clipping. If None, no
        clipping is performed. Defaults to None.
    bounds : Optional[dict], optional
        A dictionary of bounds for the model parameters. Defaults to None.
    print_info : bool, optional
        Whether to print the fit information. Defaults to True.
    **kwargs
        Additional keyword arguments passed to the fitter.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], Optional[models.CompoundModel]]
        - The fitted y values.
        - The parameters of the fitted model.
        - The covariance matrix of the parameters (if available).
        - The fitted astropy model instance.
    """
    # Sigma clipping on the data
    clipped_points_count = 0
    if sigma_clip_sigma is not None:
        y_masked: np.ma.MaskedArray = sigma_clip(y, sigma=sigma_clip_sigma, masked=True)  # pyright: ignore[reportAssignmentType]
        mask_array = ~y_masked.mask
        clipped_points_count = np.sum(y_masked.mask)
    else:
        mask_array = np.ones_like(y, dtype=bool)

    # Apply the mask to the data arrays
    data_gen = yield_mask_data(x, y, mask=mask_array)
    x_clipped, y_clipped = next(data_gen), next(data_gen)

    weights_clipped = None
    if weights is not None:
        weights_clipped = next(yield_mask_data(weights, mask=mask_array))

    if len(x_clipped) < 3:  # Not enough points to fit for 3 parameters
        print("Warning: Not enough data points to fit after sigma clipping.")
        return np.zeros_like(x), np.array([0, 0, 0]), None, None

    # Initial guess for the parameters
    amp_init = np.max(y_clipped) - np.ma.median(y_clipped)
    c_init = np.ma.median(y_clipped)

    # Define the model: Gaussian + constant, with fixed mean=0
    model_init = models.Gaussian1D(  # pyright: ignore[reportOperatorIssue]
        amplitude=amp_init, mean=0, stddev=0.3
    ) + models.Const1D(amplitude=c_init)
    model_init.mean_0.fixed = True

    # Apply bounds if provided
    if bounds:
        for param_name, bound_tuple in bounds.items():
            param = getattr(model_init, param_name)
            param.bounds = bound_tuple

    # Fitter
    fitter = fitting.TRFLSQFitter()

    # Fit the model to the clipped data
    fitted_model = fitter(
        model_init,
        x_clipped,
        y_clipped,
        weights=weights_clipped,
        **kwargs,
    )

    # Get results
    fit_y = fitted_model(x)
    para = fitted_model.parameters
    cov = None
    if fitter.fit_info.get("param_cov") is not None:
        cov = fitter.fit_info["param_cov"]

    if print_info:
        print("Parameter names:")
        print(fitted_model.param_names)
        print("Gauss fit parameters: %s" % list(para))
        if cov is not None:
            print("Gauss fit covariance (without fix mean): \n%s" % cov)
        print(f"{clipped_points_count} points clipped during sigma clipping.")
        info_fitness(y, fit_y, 3)

    return fit_y, para, cov, fitted_model


def fit_yprofile(
    data: np.ndarray,
    x_range: Sequence[float] = [-0.5, 0.5],
    xlim: Sequence[float] = [-3, 3],
    ylim: Sequence[float] = [-3, 3],
    show_fit: bool = True,
    weights: Optional[np.ndarray] = None,
    sigma_clip_sigma: Optional[float] = None,
    **kwargs,
):
    """
    Fit the y-profile of the data, averaged over a given x-range.

    Parameters
    ----------
    data : np.ndarray
        The 2D data array to fit.
    x_range : Sequence[float], optional
        The x-range over which to average the data to create the y-profile.
        Defaults to [-0.5, 0.5].
    xlim : Sequence[float], optional
        The x-limits of the data array. Defaults to [-3, 3].
    ylim : Sequence[float], optional
        The y-limits of the data array. Defaults to [-3, 3].
    show_fit : bool, optional
        If True, plot the y-profile and the fit. Defaults to True.
    weights : Optional[np.ndarray], optional
        Weights for the data, used in averaging and fitting. Defaults to None.
    sigma_clip_sigma : Optional[float], optional
        The sigma for sigma clipping during the Gaussian fit. Defaults to None.
    **kwargs
        Additional keyword arguments passed to `get_gaussian_fit`.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]
        - y-coordinates array.
        - The averaged y-profile.
        - The fitted y-profile.
        - The parameters of the fitted model.
        - The covariance matrix of the parameters.
    """
    sy, sx = data.shape

    # cut the data to the given x range to estimate the width
    x_st, x_ed = get_start_end_pixel_index(*x_range, shape=sx, lim=xlim, offset=1)
    data_cut = data[:, x_st:x_ed]

    # fit the y_profile to estimate the width
    y_profile = data_cut.mean(axis=1)
    yy = np.linspace(*ylim, num=sy)

    # get weights for y_profile fitting
    fit_weights = None
    if weights is not None:
        if weights.shape != data.shape:
            raise ValueError(
                f"Shape of weights {weights.shape} must match shape of data {data.shape}."
            )
        weights_cut = weights[:, x_st:x_ed]
        fit_weights = weights_cut.mean(axis=1)

    # Default bounds for astropy fitting
    bounds = {
        "amplitude_0": (0, np.inf),  # amp > 0
        "stddev_0": (0, 0.6),  # sigma
    }

    # Update with any user-provided bounds
    if "bounds" in kwargs:
        bounds.update(kwargs["bounds"])

    kwargs["bounds"] = bounds

    fit_y, para, cov, _ = get_gaussian_fit(
        yy, y_profile, weights=fit_weights, sigma_clip_sigma=sigma_clip_sigma, **kwargs
    )

    if show_fit:
        _ = plot_line_diff(
            yy,
            y_profile,
            fit_y,
            labels=["Y-profile", "Gaussian fit"],
            xlabel="Y",
            ylabel=r"T [$\mu$K]",
        )

    return yy, y_profile, fit_y, para, cov


def get_center_outer_background(
    data, width, x_range=[-0.5, 0.5], xlim=[-3, 3], ylim=[-3, 3], shift_unit=2
):
    """
    Extract center, outer, and background regions from the data.

    The 'center' is defined by the `x_range` and a given `width` in y.
    The 'outer' regions are shifted versions of the center region along the x-axis (left and right).
    The 'background' is taken from regions in y far from the center.

    Parameters
    ----------
    data : np.ndarray
        The 2D data array.
    width : float
        The width of the central region (i.e., 1-sigma value) .
    x_range : list, optional
        The x-range of the central region. Defaults to [-0.5, 0.5].
    xlim : list, optional
        The x-limits of the data array. Defaults to [-3, 3].
    ylim : list, optional
        The y-limits of the data array. Defaults to [-3, 3].
    shift_unit : int, optional
        The shift distance for the outer regions. Defaults to 2.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        A tuple containing the extracted regions:
        (center, outer_left, outer_right, background).
    """
    sy, sx = data.shape
    npix_unit_x = sx / (xlim[1] - xlim[0])
    npix_unit_y = sy / (ylim[1] - ylim[0])
    width_pix = get_int(width * npix_unit_y)
    print("width: %s" % width)
    print("width_pix: %s" % width_pix)

    # cut the data to the given x range to estimate the width
    x_st, x_ed = get_start_end_pixel_index(*x_range, shape=sx, lim=xlim, offset=1)
    data_cut = data[:, x_st:x_ed]

    # get the center, outer and background region
    # The center area uses the Gaussian 1-sigma width
    y_st = sy // 2 - width_pix
    y_ed = sy // 2 + width_pix
    center = data_cut[y_st:y_ed, :]

    # The background area uses the 3-4 sigma range
    bk_top_st = y_ed + width_pix * 3
    bk_top_ed = y_ed + width_pix * 4
    bk_bottom_st = y_st - width_pix * 4
    bk_bottom_ed = y_st - width_pix * 3
    bk_list = list(range(bk_top_st, bk_top_ed)) + list(
        range(bk_bottom_st, bk_bottom_ed)
    )
    background = data[bk_list, x_st:x_ed]

    # The outer areas are the center region shifted along the x-axis
    left_ids = get_int(
        [x_st - shift_unit * npix_unit_x, x_ed - shift_unit * npix_unit_x]
    )
    right_ids = get_int(
        [x_st + shift_unit * npix_unit_x, x_ed + shift_unit * npix_unit_x]
    )
    left = data[y_st:y_ed, left_ids[0] : left_ids[1]]
    right = data[y_st:y_ed, right_ids[0] : right_ids[1]]

    return center, left, right, background


def get_signal_level(
    data,
    weights: Optional[np.ndarray] = None,
    x_range=[-0.5, 0.5],
    x_num=20,
    xlim=[-3, 3],
    ylim=[-3, 3],
    shift_unit=2,
    show_fit=True,
    sigma_clip_sigma: Optional[float] = None,
    **kwargs,
):
    """
    Get the signal level by fitting a y-profile and analyzing regions.

    This function orchestrates the fitting of a y-profile to get a width,
    then uses that width to define and extract center, outer, and background
    regions to determine signal levels.

    Parameters
    ----------
    data : np.ndarray
        The 2D data array.
    weights : Optional[np.ndarray], optional
        Weights for the fitting. Defaults to None.
    x_range : list, optional
        The x-range for profile averaging and analysis. Defaults to [-0.5, 0.5].
    x_num : int, optional
        Number of points for the output x-axis. Defaults to 20.
    xlim : list, optional
        The x-limits of the data array. Defaults to [-3, 3].
    ylim : list, optional
        The y-limits of the data array. Defaults to [-3, 3].
    shift_unit : int, optional
        The shift units to get the outer regions. Defaults to 2.
    show_fit : bool, optional
        If True, plot the y-profile fit. Defaults to True.
    sigma_clip_sigma : Optional[float], optional
        Sigma for sigma clipping in the fit. Defaults to None.
    **kwargs
        Additional keyword arguments for `get_gaussian_fit`.

    Returns
    -------
    Tuple[List, List, List, np.ndarray, Optional[np.ndarray]]
        - A list of x-arrays for plotting profiles.
        - A list of y-arrays for plotting profiles.
        - A list of flattened data from the center, left, right, and background regions.
        - The parameters of the fitted y-profile model.
        - The covariance matrix of the parameters.
    """

    # fit the y profile
    yy, yprofile, yfit, paras, cov = fit_yprofile(
        data,
        x_range,
        xlim,
        ylim,
        show_fit,
        weights=weights,
        sigma_clip_sigma=sigma_clip_sigma,
        **kwargs,
    )

    # For Gaussian1D + Const1D model, parameters are [amplitude, mean, stddev, constant].
    # The width is the standard deviation, which is paras[2].
    width = paras[2] if paras is not None else 0.1
    # get the center, left, right and background
    clrb = get_center_outer_background(data, width, x_range, xlim, ylim, shift_unit)

    # get the filament signal level and background level
    tf = np.mean(clrb[0], axis=0)
    tbg = np.mean(clrb[-1], axis=0)

    # get the x axis
    xx = np.linspace(*x_range, num=x_num)

    # create the line data
    line_x = [yy, yy, xx, xx]
    line_y = [yprofile, yfit, tf, tbg]

    # get flattened data
    clrb_flat = [x.flatten() for x in clrb]

    return line_x, line_y, clrb_flat, paras, cov
