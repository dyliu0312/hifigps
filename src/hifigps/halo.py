"""
Functions to fit and subtract halo contributions.

"""

from typing import Callable, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.optimize import least_squares

from hifigps.bins import edge2center, get_id_edge, get_linbin


def cal_r(
    point: Sequence[float],
    grids: np.ndarray,
    runit: np.ndarray,
    costheta: Optional[float] = None,
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """
    Given a point, calculate the radius of each grid,
    and mask outside the given costheta range

    Args:
        point (list[int, int]): the indices of the origin point.
        grids (np.ndarray): 2D grids that need to calculate it's distance to the origin point.
        runit (np.ndarray): the unit vector used to calculate the angles of grids.
        costheta (float, Optional): the valid area are selected within -costheta to costheta.
    Returns:
        r (np.ndarray): the radius of each grid.
        v (np.ndarray): the boolen of the valid area.
    """
    point_array = np.array(point)
    vector = grids - point_array[None, :]

    r = np.sqrt(np.sum(vector**2, axis=-1))  # radius

    if costheta is not None:
        dot = vector * runit[None, :]  # dot product
        ct = dot[:, :, 0] / r  # cosine theta
        ct.ravel()  # compress into 1D array-

        # select inside the area
        valid = (ct < costheta) & (ct > -costheta)

        return r, valid

    return r


def halo_fit(
    data: np.ndarray,
    rbin_e: np.ndarray = np.arange(0, 6, 0.04),
    xlim: list = [-3, 3],  # pyright: ignore[reportMissingTypeArgument]
    ylim: Optional[list] = None,  # pyright: ignore[reportMissingTypeArgument]
    lc: list = [-1, 0],  # pyright: ignore[reportMissingTypeArgument]
    rc: list = [1, 0],  # pyright: ignore[reportMissingTypeArgument]
    costheta: float = 1 / 2,
    print_fit_res: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    use a set of random numbers to be initial profile,
    use two triangle sector only to fit,
    get best fit profile by least_squares.

    **Assuming two halo contributed equally**.

    Args:
        data (np.ndarray): 2D data to be fitted.
        rbin_e (np.ndarray): the edge of radius bins.
        xlim (list): the x range of data.
        ylim (list, Optional): the y range of data.
        lc (list, Optional): the left corner of the triangle sector.
        rc (list, Optional): the right corner of the triangle sector.
        costheta (float, Optional): the valid area are selected within -costheta to costheta.
        print_fit_res (bool, Optional): whether to print the fit result.
    Returns:
        fit_map, fit_paras, mask_bool_array

    """

    s1, s2 = data.shape
    if ylim is None:
        ylim = xlim

    # use center bin
    x_bin_c = get_linbin(xlim[0], xlim[1], s1)
    y_bin_c = get_linbin(ylim[0], ylim[1], s2)

    grids = np.zeros([s1, s2, 2])
    grids[:, :, 0] += x_bin_c[:, None]
    grids[:, :, 1] += y_bin_c[None, :]

    runit = np.array([1, 0])  # horizental right unit vector (y,x)

    # calculate radius of grids
    r1, v1 = cal_r(lc, grids, runit, costheta)
    r2, v2 = cal_r(rc, grids, runit, costheta)

    r1_indices = get_id_edge(r1.T, rbin_e)
    r2_indices = get_id_edge(r2.T, rbin_e)

    rbin_c = edge2center(rbin_e)
    r = rbin_c[: r1_indices.max() + 1]  # pyright: ignore[reportAttributeAccessIssue]
    init_paras = np.random.random(size=r.shape).astype("float32") * data.max()

    def get_fitmap(profile, r1_indices, r2_indices):
        fitmap = profile[r1_indices] + profile[r2_indices]

        return fitmap

    def error(paras, data, r1_ids, r2_ids, v1, v2):
        fitmap = get_fitmap(paras, r1_ids, r2_ids)[v1 + v2]

        return data[v1 + v2] - fitmap

    res = least_squares(
        error, init_paras, args=(data, r1_indices, r2_indices, v1.T, v2.T)
    )

    if print_fit_res:
        print(res)

    fit_paras = res.x
    fit_map = get_fitmap(fit_paras, r1_indices, r2_indices)

    return fit_map, fit_paras, ~(v1 + v2).T


def halo_fit_seprate(
    data: np.ndarray,
    rbin_e: np.ndarray = np.arange(0, 6, 0.04),
    xlim: list = [-3, 3],  # pyright: ignore[reportMissingTypeArgument]
    ylim: Optional[list] = None,  # pyright: ignore[reportMissingTypeArgument]
    lc: list = [-1, 0],  # pyright: ignore[reportMissingTypeArgument]
    rc: list = [1, 0],  # pyright: ignore[reportMissingTypeArgument]
    costheta: float = 1 / 2,
    print_fit_res: bool = True,
):
    """
    use a set of random numbers to be initial profile, use two triangle sector only to fit, get best fit profile by least_squares.
    Assuming two halo are NOT equally contributed .

    Args:
        data (np.ndarray): the data to be fitted.
        rbin_e (np.ndarray): the edge of the radius bins.
        xlim (list): the x range of the data.
        ylim (list, Optional): the y range of the data.
        lc (list): the left corner of the left triangle sector.
        rc (list): the right corner of the right triangle sector.
        costheta (float, Optional): the valid area are selected within -costheta to costheta.
        print_fit_res (bool): whether to print the result of the fitting.
    Return:
        fit_map, fit_paras, mask_bool_array
    """

    s1, s2 = data.shape
    if ylim is None:
        ylim = xlim

    # use center bin
    x_bin_c = get_linbin(xlim[0], xlim[1], s1)
    y_bin_c = get_linbin(ylim[0], ylim[1], s2)

    grids = np.zeros([s1, s2, 2])
    grids[:, :, 0] += x_bin_c[:, None]
    grids[:, :, 1] += y_bin_c[None, :]

    runit = np.array([1, 0])  # horizental right unit vector (y,x)

    # calculate radius of grids
    r1, v1 = cal_r(lc, grids, runit, costheta)
    r2, v2 = cal_r(rc, grids, runit, costheta)

    r1_indices = get_id_edge(r1.T, rbin_e)
    r2_indices = get_id_edge(r2.T, rbin_e)

    rbin_c = edge2center(rbin_e)
    # left
    r1 = rbin_c[: r1_indices.max() + 1]  # pyright: ignore[reportAttributeAccessIssue]
    init_paras1 = np.random.random(size=r1.shape).astype("float32") * data[v1].max()
    # right
    r2 = rbin_c[: r2_indices.max() + 1]  # pyright: ignore[reportAttributeAccessIssue]
    init_paras2 = np.random.random(size=r2.shape).astype("float32") * data[v2].max()

    init_paras = [*init_paras1, *init_paras2]

    def get_fitmap(profiles, r1_indices, r2_indices):
        sep = len(profiles) // 2

        profile1 = profiles[:sep]
        profile2 = profiles[sep:]

        fitmap = profile1[r1_indices] + profile2[r2_indices]

        return fitmap

    def error(paras, data, r1_ids, r2_ids, v1, v2):
        fitmap = get_fitmap(paras, r1_ids, r2_ids)

        return data[v1 + v2] - fitmap[v1 + v2]

    res = least_squares(
        error, init_paras, args=(data, r1_indices, r2_indices, v1.T, v2.T)
    )

    if print_fit_res:
        print(res)

    fit_paras = res.x
    fit_map = get_fitmap(fit_paras, r1_indices, r2_indices)

    return fit_map, fit_paras, ~(v1 + v2).T


def halo_subtract(
    data: Union[np.ndarray, Sequence[np.ndarray]],
    func_halo_fit: Callable[..., Tuple[np.ndarray, np.ndarray, np.ndarray]] = halo_fit,
    rbin_e: np.ndarray = np.arange(0, 6, 0.04),
    lim: Sequence[float] = [-3, 3],
    costheta: float = np.cos(np.pi / 4),
    print_fit_res: bool = False,
    **kwargs,
):
    """
    for each data, fit a halo profile, and subtract it from the data.
    Args:
        data (list or np.ndarray): the data to be fitted and subtracted.
        func_halo_fit (Callable): the function to fit the halo profile.
        rbin_e (np.ndarray): the edge of the radius bins.
        lim (list[float,float]): the x range of the data.
        costheta (float): the valid area are selected within -costheta to costheta.
        print_fit_res (bool): whether to print the result of the fitting.
    Return:
        fitted_map_list, residuals_list, mask_bool_array
    """
    is_single = isinstance(data, np.ndarray)
    if is_single:
        data = [data]

    fit = []
    res = []
    fit_mask = None

    for d in data:
        fitmap, _, curr_fit_mask = func_halo_fit(
            d, rbin_e, lim, costheta=costheta, print_fit_res=print_fit_res, **kwargs
        )
        fit.append(fitmap)
        res.append(d - fitmap)
        if fit_mask is None:
            fit_mask = curr_fit_mask

    if is_single:
        return fit[0], res[0], fit_mask
    else:
        return fit, res, fit_mask
