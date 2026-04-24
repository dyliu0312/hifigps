"""
estimate the signal level baed on the given width.
"""

from typing import List, Optional, Sequence, Tuple

import numpy as np


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


def get_center_outer_background(
    data, data_err, width, x_range=[-0.5, 0.5], xlim=[-3, 3], ylim=[-3, 3], shift_unit=2
):
    """
    Extract center, outer, and background regions from the data.

    The 'center' is defined by the `x_range` and a given `width` in y.
    The 'outer' regions are shifted versions of the center region along the x-axis.
    The 'background' is taken from regions in y far from the center.

    Parameters
    ----------
    data : np.ndarray
        The 2D data array.
    data_err : np.ndarray
        The 1-sigma error of the 2D data array.
    width : float
        The half-width (e.g., 1-sigma) of the central region in world coordinates.
    x_range : list, optional
        The x-range of the central region. Defaults to [-0.5, 0.5].
    xlim : list, optional
        The x-limits of the data array. Defaults to [-3, 3].
    ylim : list, optional
        The y-limits of the data array. Defaults to [-3, 3].
    shift_unit : int, optional
        The shift distance for the outer regions, in world coordinates.
        Defaults to 2.

    Returns
    -------
    Tuple[Tuple, Tuple]
        A tuple containing the extracted regions and their errors:
        ((center, left, right, background), (center_err, left_err, right_err, background_err)).
    """
    sy, sx = data.shape
    npix_unit_x = sx / (xlim[1] - xlim[0])
    npix_unit_y = sy / (ylim[1] - ylim[0])
    width_pix = get_int(width * npix_unit_y)
    # print("width: %s" % width)
    # print("width_pix: %s" % width_pix)

    # cut the data to the given x range to estimate the width
    x_st, x_ed = get_start_end_pixel_index(*x_range, shape=sx, lim=xlim, offset=1)
    data_cut = data[:, x_st:x_ed]
    data_err_cut = data_err[:, x_st:x_ed]

    # get the center, outer and background region
    # The center area uses the Gaussian 1-sigma width
    y_st = sy // 2 - width_pix
    y_ed = sy // 2 + width_pix
    center = data_cut[y_st:y_ed, :]
    center_err = data_err_cut[y_st:y_ed, :]

    # The background area uses the 3-4 sigma range
    bk_top_st = y_ed + width_pix * 3
    bk_top_ed = y_ed + width_pix * 4
    bk_bottom_st = y_st - width_pix * 4
    bk_bottom_ed = y_st - width_pix * 3
    bk_list = list(range(bk_top_st, bk_top_ed)) + list(
        range(bk_bottom_st, bk_bottom_ed)
    )
    background = data[bk_list, x_st:x_ed]
    background_err = data_err[bk_list, x_st:x_ed]

    # The outer areas are the center region shifted along the x-axis
    left_ids = get_int(
        [x_st - shift_unit * npix_unit_x, x_ed - shift_unit * npix_unit_x]
    )
    right_ids = get_int(
        [x_st + shift_unit * npix_unit_x, x_ed + shift_unit * npix_unit_x]
    )
    left = data[y_st:y_ed, left_ids[0] : left_ids[1]]
    left_err = data_err[y_st:y_ed, left_ids[0] : left_ids[1]]
    right = data[y_st:y_ed, right_ids[0] : right_ids[1]]
    right_err = data_err[y_st:y_ed, right_ids[0] : right_ids[1]]

    clrb = (center, left, right, background)
    clrb_err = (center_err, left_err, right_err, background_err)

    return clrb, clrb_err


def get_signal_level_fixwidth(
    data: np.ndarray,
    width: float,
    data_err: Optional[np.ndarray] = None,
    x_range: Sequence[float] = [-0.5, 0.5],
    x_num: int = 20,
    xlim: Sequence[float] = [-3, 3],
    ylim: Sequence[float] = [-3, 3],
    shift_unit: int = 2,
):
    """
    Get the signal level with a given width.

    This function uses a given width to define and extract center, outer,
    and background regions to determine signal levels. It does not perform
    any fitting.

    Parameters
    ----------
    data : np.ndarray
        The 2D data array.
    width : float
        The half-width (e.g., 1-sigma) of the central region in world coordinates.
    data_err : Optional[np.ndarray], optional
        1-sigma error of the data. Defaults to None.
    x_range : Sequence[float], optional
        The x-range for profile averaging and analysis. Defaults to [-0.5, 0.5].
    x_num : int, optional
        Number of points for the output x-axis. Defaults to 20.
    xlim : Sequence[float], optional
        The x-limits of the data array. Defaults to [-3, 3].
    ylim : Sequence[float], optional
        The y-limits of the data array. Defaults to [-3, 3].
    shift_unit : int, optional
        The shift distance for the outer regions. Defaults to 2.
    show_profile : bool, optional
        If True, plot the y-profile. Defaults to True.

    Returns
    -------
    Tuple
        - A list of x-arrays for plotting profiles.
        - A list of y-arrays for plotting profiles.
        - A list of y-error arrays for plotting profiles.
        - A list of flattened data from the center, left, right, and background regions.
        - A list of flattened error data from the CLRB regions.
    """
    has_err = data_err is not None
    if not has_err:
        data_err = np.zeros_like(data)

    sy, sx = data.shape

    # Calculate y-profile for inspection
    x_st, x_ed = get_start_end_pixel_index(*x_range, shape=sx, lim=xlim, offset=1)
    data_cut = data[:, x_st:x_ed]
    y_profile = data_cut.mean(axis=1)
    yy = np.linspace(*ylim, num=sy)

    y_profile_err = None
    if has_err:
        if data_err.shape != data.shape:
            raise ValueError(
                f"Shape of data_err {data_err.shape} must match shape of data {data.shape}."
            )
        data_err_cut = data_err[:, x_st:x_ed]
        y_profile_err = np.sqrt(np.sum(data_err_cut**2, axis=1)) / data_err_cut.shape[1]

    # get the center, left, right and background
    clrb, clrb_err = get_center_outer_background(
        data, data_err, width, x_range, xlim, ylim, shift_unit
    )

    # get the filament signal level and background level
    tf = np.mean(clrb[0], axis=0)
    tbg = np.mean(clrb[-1], axis=0)

    # propagate errors for tf and tbg
    tf_err = np.sqrt(np.sum(clrb_err[0] ** 2, axis=0)) / clrb_err[0].shape[0]
    tbg_err = np.sqrt(np.sum(clrb_err[-1] ** 2, axis=0)) / clrb_err[-1].shape[0]

    # get the x axis
    xx = np.linspace(*x_range, num=x_num)

    # create the line data
    line_x = [yy, xx, xx]
    line_y = [y_profile, tf, tbg]
    line_y_err = [y_profile_err, tf_err, tbg_err]

    # get flattened data
    clrb_flat = [x.flatten() for x in clrb]
    clrb_err_flat = [x.flatten() for x in clrb_err]

    return line_x, line_y, line_y_err, clrb_flat, clrb_err_flat


def get_signal_level_fixwidth_multi(
    datas: Sequence[np.ndarray],
    width: float,
    data_errs: Optional[Sequence[np.ndarray]] = None,
    x_range: Sequence[float] = [-0.5, 0.5],
    x_num: int = 20,
    xlim: Sequence[float] = [-3, 3],
    ylim: Sequence[float] = [-3, 3],
    shift_unit: int = 2,
):
    """
    Get the signal level with a given width for multiple datasets.

    This function processes each 2D data array in `datas` independently
    using `get_signal_level_fixwidth` and returns a list of results
    for each dataset.

    Parameters
    ----------
    datas : Sequence[np.ndarray]
        A list of 2D data arrays.
    width : float
        The half-width (e.g., 1-sigma) of the central region in world coordinates.
    data_errs : Optional[Sequence[np.ndarray]], optional
        A list of 1-sigma error arrays for the data. Must have the same
        length as `datas`. Defaults to None.
    x_range : Sequence[float], optional
        The x-range for profile averaging and analysis. Defaults to [-0.5, 0.5].
    x_num : int, optional
        Number of points for the output x-axis. Defaults to 20.
    xlim : Sequence[float], optional
        The x-limits of the data array. Defaults to [-3, 3].
    ylim : Sequence[float], optional
        The y-limits of the data array. Defaults to [-3, 3].
    shift_unit : int, optional
        The shift distance for the outer regions. Defaults to 2.

    Returns
    -------
    Tuple[list, list, list, list, list]
        A tuple of lists, where each inner list contains the corresponding
        return value from `get_signal_level_fixwidth` for each input dataset.
        - List of line_x arrays.
        - List of line_y arrays.
        - List of line_y_err arrays.
        - List of clrb_flat arrays.
        - List of clrb_err_flat arrays.
    """
    results = []

    if data_errs is None:
        processed_data_errs = [None] * len(datas)
    else:
        processed_data_errs = data_errs

    if len(datas) != len(processed_data_errs):
        raise ValueError("datas and data_errs must have the same length.")

    if not datas:
        return [], [], [], [], []

    for i, data in enumerate(datas):
        err_data = processed_data_errs[i]

        line_x, line_y, line_y_err, clrb_flat, clrb_err_flat = (
            get_signal_level_fixwidth(
                data=data,
                width=width,
                data_err=err_data,
                x_range=x_range,
                x_num=x_num,
                xlim=xlim,
                ylim=ylim,
                shift_unit=shift_unit,
            )
        )
        results.append((line_x, line_y, line_y_err, clrb_flat, clrb_err_flat))

    transposed_results = list(zip(*results))
    return (
        transposed_results[0],
        transposed_results[1],
        transposed_results[2],
        transposed_results[3],
        transposed_results[4],
    )


def get_averaged_signal_levels(
    line_y_list: Sequence[Sequence[np.ndarray]],
    line_y_err_list: Sequence[Sequence[Optional[np.ndarray]]],
) -> Tuple[List[float], List[float], List[float], List[float]]:
    """
    Extracts and averages the signal level (tf) and background level (tbg)
    for each dataset from the outputs of `get_signal_level_fixwidth_multi`.

    Parameters
    ----------
    line_y_list : Sequence[Sequence[np.ndarray]]
        A list, where each element is a sequence containing [y_profile, tf, tbg]
        for a single dataset, and tf/tbg are 1D arrays.
    line_y_err_list : Sequence[Sequence[Optional[np.ndarray]]]
        A list, where each element is a sequence containing [y_profile_err, tf_err, tbg_err]
        for a single dataset, and tf_err/tbg_err are 1D arrays representing errors.
        y_profile_err can be None.

    Returns
    -------
    Tuple[List[float], List[float], List[float], List[float]]
        A tuple containing four lists:
        - all_tfs: A list of averaged signal levels, one for each dataset.
        - all_tf_errs: A list of errors for the averaged signal levels, one for each dataset.
        - all_tbg: A list of averaged background levels, one for each dataset.
        - all_tbg_errs: A list of errors for the averaged background levels, one for each dataset.
    """

    all_tfs = []
    all_tbg = []
    all_tf_errs = []
    all_tbg_errs = []

    for line_y, line_y_err in zip(line_y_list, line_y_err_list):
        tf = line_y[1]
        tbg = line_y[2]
        tf_err = line_y_err[1]
        tbg_err = line_y_err[2]

        all_tfs.append(np.mean(tf))
        all_tbg.append(np.mean(tbg))

        if tf_err is not None:
            if np.all(tf_err == 0):
                all_tf_errs.append(0.0)
            else:
                all_tf_errs.append(np.sqrt(np.sum(tf_err**2)) / len(tf_err))
        else:
            all_tf_errs.append(np.nan)

        if tbg_err is not None:
            if np.all(tbg_err == 0):
                all_tbg_errs.append(0.0)
            else:
                all_tbg_errs.append(np.sqrt(np.sum(tbg_err**2)) / len(tbg_err))
        else:
            all_tbg_errs.append(np.nan)

    return all_tfs, all_tf_errs, all_tbg, all_tbg_errs


def get_null_test_results(
    datas: Sequence[np.ndarray],
    width: float,
    x_range: Sequence[float] = [-0.5, 0.5],
    xlim: Sequence[float] = [-3, 3],
    ylim: Sequence[float] = [-3, 3],
    shift_unit: int = 2,
) -> Tuple[List[float], List[float]]:
    """
    Get the overall values of null-test results corresponding to the
    signal (tf) and background (tbg) regions across multiple realizations.

    This function processes each dataset to extract the center (signal) and
    background regions, combines all pixel values from these regions across
    all datasets.

    Parameters
    ----------
    datas : Sequence[np.ndarray]
        A list of 2D data arrays of the null tests.
    width : float
        The width of the central region (i.e., 1-sigma value).
    x_range : Sequence[float], optional
        The x-range of the central region. Defaults to [-0.5, 0.5].
    xlim : Sequence[float], optional
        The x-limits of the data array. Defaults to [-3, 3].
    ylim : Sequence[float], optional
        The y-limits of the data array. Defaults to [-3, 3].
    shift_unit : int, optional
        The shift distance for the outer regions. Defaults to 2.

    Returns
    -------
    Tuple[List[float], List[float]]
        A tuple containing:
        - tf_mean: The mean of all pixel values in the signal region across all datasets.
        - tf_std: The standard deviation of all pixel values in the signal region.
        - tbg_mean: The mean of all pixel values in the background region.
        - tbg_std: The standard deviation of all pixel values in the background region.
    """
    all_tf_pixels = []
    all_tbg_pixels = []

    for data in datas:
        clrb, _ = get_center_outer_background(
            data,
            data,  # placeholder of data_error
            width,
            x_range=x_range,
            xlim=xlim,
            ylim=ylim,
            shift_unit=shift_unit,
        )
        # clrb[0] is the center (tf) region, clrb[3] is the background (tbg) region
        all_tf_pixels.extend(clrb[0].flatten())
        all_tbg_pixels.extend(clrb[3].flatten())

    return all_tf_pixels, all_tbg_pixels
