"""
Custom functions to plot the stack results.
"""

from typing import (
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from numpy.typing import NDArray

from hifigps.plot import (
    arg_list,
    make_figure,
    plot_heatmaps,
    plot_hist,
    plot_line,
    plot_line_diff,
)


def plot_stack_fit_res(
    data: List[NDArray],
    axes: Optional[List[Axes]] = None,
    cmap: Union[str, List[str]] = ["viridis", "viridis", "RdBu_r"],
    norm: Union[str, List[str]] = "linear",
    vmin: Union[float, List[float], None] = None,
    vmax: Union[float, List[float], None] = None,
    q: float = 5,
    show_cbar: bool = True,
    cbar_label: Union[str, List[str]] = r"$T\ [\mu$K]",
    kw_makefigure: Optional[dict] = dict(sharey=True),  # pyright: ignore[reportMissingTypeArgument]
    title: Union[str, List[str]] = [
        "Pairwise-stacked map",
        "Fitted halo contribution",
        "Filament signal",
    ],
    xlabel: Union[str, List[str]] = "X",
    ylabel: Union[str, List[str]] = "Y",
) -> List[Axes]:
    """
    plot stack, fit and residual map for comparison.
    """

    if axes is None:
        _, axes = make_figure(1, 3, figsize=(12, 3), **kw_makefigure)  # pyright: ignore[reportCallIssue]

    axes = plot_heatmaps(
        data,
        axes=axes,
        cmap=cmap,
        norm=norm,
        vmin=vmin,
        vmax=vmax,
        q=q,
        show_cbar=show_cbar,
        cbar_label=cbar_label,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
    )

    return axes


def plot_res(
    data: List[NDArray],
    axes: Optional[List[Axes]] = None,
    norm: Union[str, List[str]] = "linear",
    vmin: Union[float, List[float]] = -10,
    vmax: Union[float, List[float]] = 40,
    q: float = 5,
    tick_in: bool = True,
    show_cbar: Union[bool, List[bool]] = [False, True],
    cmap: Union[str, List[str]] = "viridis",
    cbar_loc: str = "right",
    cbar_label: Union[str, List[str]] = r"$T\ [\mu$K]",
    title: Union[str, List[str]] = ["HI only", "HI + noise"],
    xlabel: Union[str, List[str]] = "X",
    ylabel: Union[str, List[str]] = "Y",
    kw_makefigure: Optional[dict] = dict(  # pyright: ignore[reportMissingTypeArgument]
        figsize=(6, 3),
        wspace=0.0,
        hspace=0.0,
    ),
    kw_pcolormesh: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
    kw_cbar: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
) -> List[Axes]:
    """
    Plot residual maps (HI only and HI + noise) for comparison.
    """

    if axes is None:
        _, axes = make_figure(1, 2, **kw_makefigure)  # pyright: ignore[reportCallIssue]
    if axes is None:
        raise ValueError("axes must be provided")

    plot_heatmaps(
        data,
        axes,
        cmap,
        norm,
        vmin,
        vmax,
        q,
        tick_in,
        show_cbar,
        cbar_loc,  # pyright: ignore[reportArgumentType]
        cbar_label,  # type: ignore
        title,
        xlabel,
        ylabel,
        kw_pcolormesh=kw_pcolormesh,
        kw_cbar=kw_cbar,
    )

    return axes


def plot_profile_2c(
    x: List[Union[NDArray, Sequence[float]]],
    y: List[Union[NDArray, Sequence[float]]],
    cut: int = 2,
    axes: Optional[List[Axes]] = None,
    y_err: Optional[List[Optional[Union[NDArray, Sequence[float]]]]] = None,
    x_err: Optional[List[Optional[Union[NDArray, Sequence[float]]]]] = None,
    text_pos: Optional[List[List[float]]] = None,
    width: float = 0.32,
    fontsize: Union[float, List[float]] = 12,
    xlabel: Union[str, List[Union[str, None]], None] = ["Y", "X"],
    ylabel: Union[str, List[Union[str, None]], None] = [r"$T\ [\mu$K]", None],
    title: Union[str, List[Union[str, None]], None] = [
        "Transverse-section profile",
        "Lengthwise-section profile",
    ],
    color: Union[str, List[str], None] = ["b", "r", "r", "k"],
    marker: Union[str, List[str], None] = [".", "", "o", "o"],
    linestyle: Union[str, List[str], None] = ["", "--", "--", "--"],
    label: Union[str, List[Union[str, None]], None] = [
        r"$|{\rm X}|< 0.5$",
        "Gaussian fit",
        r"$T_{\rm f}$",
        r"$T_{\rm bg}$",
    ],
    tick_in: bool = True,
) -> List[Axes]:
    """
    plot profiles in two column subplots.

    Parameters:
        x (list): x data points.
        y (list): y data points.
        cut (int): separate the data points into two parts, \
            e.g. cut=2 means the data points in the first two rows are plotted in the first subplot, \
                and rest are plotted in the second subplot.
        axes (Axes, optional): axes object. Defaults to None.
        y_err (list, optional): y error of the data. Defaults to None.
        x_err (list, optional): x error of the data. Defaults to None.
        text_pos (List[list], optional): text position. Defaults to None, and the text position will be auto.
        width (float): width of the filament to show on the left panel. Defaults to 0.32.
        fontsize (str): fontsize of the text. Defaults to 12.
        xlabel (str, list): x-axis labels. Defaults to ["Y", "X"].
        ylabel (str, list): y-axis labels. Defaults to [r"$T$ [$\\mu$K]", None].
        title (str, list): subplot titles. Defaults to ["Transverse-section profile", \
            "Lengthwise-section profile"].
        color (str, list): line colors. Defaults to ["b", "r", "r", "k"].
        marker (str, list): line markers. Defaults to [".", "", "o", "o"].
        linestyle (str, list): line styles. Defaults to ["", "--", "--", "--"].
        label (str, list): line labels. Defaults to [r"$|{\\rm X}|< 0.5$", "Gaussian fit", \
            r"$T_{\rm f}$", r"$T_{\rm bg$"].
        tick_in (bool): whether to show ticks inside the axes. Defaults to True.
    Returns:
        axes (list): axes object.

    """
    if axes is None:
        _, axes = make_figure(
            1,
            2,
            figsize=(8, 3),
            sharey=True,
            sharex=False,
            wspace=0,
            aspect=None,
        )
    if axes is None:
        raise ValueError("axes is None.")

    try:
        x_left, x_right = x[:cut], x[cut:]
        y_left, y_right = y[:cut], y[cut:]
    except IndexError as ie:
        raise ValueError("cut index is out of the range.") from ie

    y_err_left, y_err_right = None, None
    if y_err:
        y_err_left, y_err_right = y_err[:cut], y_err[cut:]

    x_err_left, x_err_right = None, None
    if x_err:
        x_err_left, x_err_right = x_err[:cut], x_err[cut:]

    nline = len(x)
    xlabel = arg_list(xlabel, 2)
    ylabel = arg_list(ylabel, 2)
    title = arg_list(title, 2)

    color = arg_list(color, nline)
    marker = arg_list(marker, nline)
    linestyle = arg_list(linestyle, nline)
    label = arg_list(label, nline)
    fontsize = arg_list(fontsize, 3)

    for i, ax in enumerate(list(axes)):  # type: ignore
        c = color[:cut] if i == 0 else color[cut:]
        m = marker[:cut] if i == 0 else marker[cut:]
        ls = linestyle[:cut] if i == 0 else linestyle[cut:]
        lb = label[:cut] if i == 0 else label[cut:]
        plot_line(
            x_left if i == 0 else x_right,
            y_left if i == 0 else y_right,
            ax=ax,
            y_err=y_err_left if i == 0 else y_err_right,
            x_err=x_err_left if i == 0 else x_err_right,
            xlabel=xlabel[i],
            ylabel=ylabel[i],
            title=title[i],
            color=c,
            marker=m,
            linestyle=ls,
            label=lb,
            tick_in=tick_in,
        )
        # add horizontal lines that marker T = 0
        ax.axhline(0, linestyle="--", c="gray")
        # add text that marker the mean values of filament and background.
        if i == 1:
            # filament
            mf = np.mean(y_right[0])
            if y_err_right and y_err_right[0] is not None:
                mf_err = np.sqrt(np.sum(y_err_right[0] ** 2)) / len(y_right[0])  # pyright: ignore[reportOperatorIssue]
            else:
                mf_err = np.std(y_right[0]) / np.sqrt(len(y_right[0]))

            strf = r"$<T_{\rm f}>$ = %.2f $\pm$ %.2f $\mu$K" % (mf, mf_err)
            if text_pos is None:
                ax.text(
                    0.05,
                    0.9,
                    strf,
                    color="r",
                    fontsize=fontsize[0],
                    transform=ax.transAxes,
                )
            else:
                ax.text(*text_pos[0], s=strf, color="r", fontsize=fontsize[0])

            # background
            mb = np.mean(y_right[1])
            if y_err_right and y_err_right[1] is not None:
                mb_err = np.sqrt(np.sum(y_err_right[1] ** 2)) / len(y_right[1])  # pyright: ignore[reportOperatorIssue]
            else:
                mb_err = np.std(y_right[1]) / np.sqrt(len(y_right[1]))

            strb = r"$<T_{\rm bg}>$ = %.2f $\pm$ %.2f $\mu$K" % (mb, mb_err)
            if text_pos is None:
                ax.text(
                    0.05,
                    0.8,
                    strb,
                    color="k",
                    fontsize=fontsize[1],
                    transform=ax.transAxes,
                )
            else:
                ax.text(*text_pos[1], s=strb, color="k", fontsize=fontsize[1])

        # add text that marker the width of filament
        else:
            strw = r"$w_{\rm t} = %.2f$" % width
            if text_pos is None:
                ax.text(
                    0.05,
                    0.9,
                    strw,
                    color="r",
                    fontsize=fontsize[2],
                    transform=ax.transAxes,
                )
            else:
                ax.text(
                    *text_pos[2],
                    s=strw,
                    color="r",
                    fontsize=fontsize[2],
                )

    return axes


def plot_profile_2r(
    x: List[List[Union[NDArray, Sequence[float]]]],
    y: List[List[Union[NDArray, Sequence[float]]]],
    cut: int = 2,
    axes: Union[Tuple[Axes], List[Axes], None] = None,
    y_err: Optional[List[List[Optional[Union[NDArray, Sequence[float]]]]]] = None,
    x_err: Optional[List[List[Optional[Union[NDArray, Sequence[float]]]]]] = None,
    text_pos: Optional[List[List[float]]] = None,
    width: Union[float, List[float]] = 0.32,
    fontsize: Union[float, List[float]] = 12,
    title: Union[str, List[str], None] = None,
    xlabel: Union[str, List[Union[str, None]], None] = ["Y", "X"],
    ylabel: Union[str, List[Union[str, None]], None] = r"$T\ [\mu$K]",
    color: Union[str, List[str], None] = None,
    marker: Union[str, List[str], None] = None,
    linestyle: Union[str, List[str]] = "--",
    label: Union[str, List[str], None] = None,
    tick_in: bool = True,
) -> Union[List[Axes], Tuple[Axes]]:
    """
    plot lines in two row subplots.

    Args:
        x (List[list]): x data points.
        y (List[list]): y data points.
        cut (int): cut index to separate the data into upper and bottom two parts. Defaults to 2. \
            i.e. The upper and bottom panels of the first subplot is showing x[0, :cut] and x[0, cut:] respectively.
        axes (List[Axes], Tuple[Axes], None): axes object. Defaults to None.
        y_err (list, optional): y error of the data. Defaults to None.
        x_err (list, optional): x error of the data. Defaults to None.
        text_pos (list, optional): text positions. Defaults to None, and the text position will be auto.
        width (float, list): width of the filament to show on the upper pannel. Defaults to 0.32.
        fontsize (float): fontsize of the text. Defaults to 12.
        title (str, list): titles of the subplots. Defaults to None.
        xlabel_up (str, list): x labels of the upper pannels. Defaults to "Y".
        xlabel_down (str, list): x labels of the lower pannels. Defaults to "X". Set to None to same as xlabel_up.
        ylabel_up (str, list): y labels of the upper pannels. Defaults to r"T [$\\mu$K]".
        ylabel_down (str, list): y labels of the lower pannels. Defaults to None to same as ylabel_up.
        color (str, list): colors of the lines. Defaults to None.
        marker (str, list): markers of the lines. Defaults to None.
        linestyle (str, list): linestyles of the lines. Defaults to '--'.
        label (str, list): labels of the lines. Defaults to None.
        tick_in (bool): whether to show the ticks inside the subplots. Defaults to True.
    """

    ncol = len(x)
    nline = len(x[0])
    nrow = 2
    ntext = 3

    if axes is None:
        _, axes = make_figure(
            2,
            ncol,
            figsize=(4 * ncol + 2, 8),
            sharey="row",
            sharex="row",
            wspace=0,
            hspace=0.4,
            aspect=None,
        )
    if axes is None:
        raise ValueError("axes is None.")

    ## chenck line parameters
    color = arg_list(color, nline, ncol)
    marker = arg_list(marker, nline, ncol)
    linestyle = arg_list(linestyle, nline, ncol)
    label = arg_list(label, nline, ncol)

    # check text parameters
    if text_pos is not None:
        text_pos = arg_list(text_pos, ntext, ncol)
    fontsize = arg_list(fontsize, ntext, ncol)
    width = arg_list(width, ncol)

    # check xlabels
    xlabel = arg_list(xlabel, nrow * ncol)
    ylabel = arg_list(ylabel, nrow * ncol)
    title = arg_list(title, nrow * ncol)

    for i in range(ncol):
        st, ed = i * nline, (i + 1) * nline
        try:
            x_up, x_down = x[i][:cut], x[i][cut:]
            y_up, y_down = y[i][:cut], y[i][cut:]
            c_up, c_down = color[st:ed][:cut], color[st:ed][cut:]
            m_up, m_down = marker[st:ed][:cut], marker[st:ed][cut:]
            ls_up, ls_down = linestyle[st:ed][:cut], linestyle[st:ed][cut:]
            la_up, la_down = label[st:ed][:cut], label[st:ed][cut:]

            y_err_up, y_err_down = None, None
            if y_err and y_err[i]:
                y_err_up, y_err_down = y_err[i][:cut], y_err[i][cut:]

            x_err_up, x_err_down = None, None
            if x_err and x_err[i]:
                x_err_up, x_err_down = x_err[i][:cut], x_err[i][cut:]

            _text_pos = (
                text_pos[i * ntext : (i + 1) * ntext] if text_pos is not None else None
            )
            _xlabel = xlabel[i * nrow : (i + 1) * nrow]
            _ylabel = ylabel[i * nrow : (i + 1) * nrow]
            _title = title[i * nrow : (i + 1) * nrow]

        except IndexError as ie:
            raise ValueError("cut index is out of the range.") from ie

        # plot upper pannel
        ax_up = axes[0, i] if ncol > 1 else axes[0]  # pyright: ignore[reportArgumentType, reportCallIssue]
        plot_line(
            x_up,
            y_up,
            ax_up,
            y_err=y_err_up,
            x_err=x_err_up,
            xlabel=_xlabel[0],
            ylabel=_ylabel[0],
            title=_title[0],
            color=c_up,
            marker=m_up,
            linestyle=ls_up,
            label=la_up,
            tick_in=tick_in,
        )
        ax_up.axhline(0, linestyle="--", c="gray")
        strw = r"$w_{\rm t} = %.2f$" % width[i]
        if _text_pos is None:
            ax_up.text(
                0.05,
                0.9,
                strw,
                color="r",
                fontsize=fontsize[0],
                transform=ax_up.transAxes,
            )
        else:
            ax_up.text(*_text_pos[2], s=strw, color="r", fontsize=fontsize[0])

        # plot lower pannel
        ax_down = axes[1, i] if ncol > 1 else axes[1]  # pyright: ignore[reportGeneralTypeIssues, reportArgumentType, reportCallIssue]
        plot_line(
            x_down,
            y_down,
            ax_down,
            y_err=y_err_down,
            x_err=x_err_down,
            xlabel=_xlabel[1],
            ylabel=_ylabel[1],
            title=_title[1],
            color=c_down,
            marker=m_down,
            linestyle=ls_down,
            label=la_down,
            tick_in=tick_in,
        )
        ax_down.axhline(0, linestyle="--", c="gray")

        # filament
        mf = np.mean(y_down[0])
        if y_err_down and y_err_down[0] is not None:
            mf_err = np.sqrt(np.sum(y_err_down[0] ** 2)) / len(y_down[0])  # pyright: ignore[reportOperatorIssue]
        else:
            mf_err = np.std(y_down[0]) / np.sqrt(len(y_down[0]))

        strf = r"$<T_{\rm f}>$ = %.2f $\pm$ %.2f $\mu$K" % (mf, mf_err)
        if _text_pos is None:
            ax_down.text(
                0.05,
                0.9,
                strf,
                color="r",
                fontsize=fontsize[1],
                transform=ax_down.transAxes,
            )
        else:
            ax_down.text(*_text_pos[0], s=strf, color="r", fontsize=fontsize[1])

        # background
        mb = np.mean(y_down[1])
        if y_err_down and y_err_down[1] is not None:
            mb_err = np.sqrt(np.sum(y_err_down[1] ** 2)) / len(y_down[1])  # pyright: ignore[reportOperatorIssue]
        else:
            mb_err = np.std(y_down[1]) / np.sqrt(len(y_down[1]))

        strb = r"$<T_{\rm bg}>$ = %.2f $\pm$ %.2f $\mu$K" % (mb, mb_err)
        if _text_pos is None:
            ax_down.text(
                0.05,
                0.8,
                strb,
                color="k",
                fontsize=fontsize[2],
                transform=ax_down.transAxes,
            )
        else:
            ax_down.text(*_text_pos[1], s=strb, color="k", fontsize=fontsize[2])

    return axes


def plot_profile_2c2r(
    x: List[List[Union[NDArray, Sequence[float]]]],
    y: List[List[Union[NDArray, Sequence[float]]]],
    cut: int = 2,
    axes: Union[Tuple[Axes], None] = None,
    y_err: Optional[List[List[Optional[Union[NDArray, Sequence[float]]]]]] = None,
    x_err: Optional[List[List[Optional[Union[NDArray, Sequence[float]]]]]] = None,
    text_pos: Optional[List[List[float]]] = None,
    width: Union[float, List[float]] = [0.32] * 2,
    fontsize: Union[float, List[float]] = [12.0] * 3,
    xlabel: Union[str, List[Union[str, None]], None] = [None, None, "Y", "X"],
    ylabel: Union[str, List[Union[str, None]], None] = [
        r"$T\ [\mu$K] (HI only)",
        None,
        r"$T\ [\mu$K] (HI + noise)",
        None,
    ],
    title: List[Union[str, None]] = [
        "Transverse-section profile",
        "Lengthwise-section profile",
        None,
    ],
    color: Union[str, List[str]] = ["b", "r", "r", "k"],
    marker: Union[str, List[str]] = [".", "", "o", "o"],
    linestyle: Union[str, List[str]] = ["", "--", "--", "--"],
    label: Union[str, List[Union[str, None]], None] = [
        r"$|{\rm X}|< 0.5$",
        "Gaussian fit",
        r"$T_{\rm f}$",
        r"$T_{\rm bg}$",
        None,
    ],
    tick_in: bool = True,
) -> Tuple[Axes]:
    """
    plot lines in two column two row subplots. The top and bottom rows are shared the same args.

    Args:
        x (List[list]): x data points.
        y (List[list]): y data points.
        cut (int): cut index. Defaults to 2.
        axes (Axes, optional): axes object. Defaults to None to create a new figure.
        y_err (list, optional): y error of the data. Defaults to None.
        x_err (list, optional): x error of the data. Defaults to None.
        text_pos (List[list], optional): text position. Defaults to None, and the text position will be auto.
        width (float): width of the filament.
        fontsize (str): fontsize of the text.
        xlabel (List[list]): x-axis labels.
        ylabel (List[list]): y-axis labels.
        title (List[list]): subplot titles.
        color (list): line colors.
        marker (list): line markers.
        linestyle (list): line styles.
        label (list): line labels.
        tick_in (bool): whether to show the ticks inside the subplots. Defaults to True.
    Returns:
        axes (Tuple[Axes]): axes object.
    """

    if axes is None:
        _, axes = make_figure(
            2,
            2,
            figsize=(8, 6),
            sharex="col",
            sharey=True,
            wspace=0,
            hspace=0,
            aspect=None,
        )
    if axes is None:
        raise ValueError("axes must be specified.")

    ncol = nrow = len(x)
    nline = len(x[0])
    ntext = 3

    if nrow != 2:
        raise ValueError("This function only supports two column two row subplots.")

    ## chenck line parameters
    color = arg_list(color, nline, nrow)
    marker = arg_list(marker, nline, nrow)
    linestyle = arg_list(linestyle, nline, nrow)

    label = arg_list(label, nline * nrow)

    # check text parameters
    if text_pos is not None:
        text_pos = arg_list(text_pos, ntext, nrow)
    width = arg_list(width, nrow)
    fontsize = arg_list(fontsize, nrow)

    # check xlabels
    xlabel = arg_list(xlabel, ncol, nrow)
    ylabel = arg_list(ylabel, ncol, nrow)
    title = arg_list(title, ncol * nrow)

    for i in range(nrow):
        y_err_row = y_err[i] if y_err else None
        x_err_row = x_err[i] if x_err else None
        plot_profile_2c(
            x[i],
            y[i],
            cut=cut,
            axes=list(axes)[i],  # pyright: ignore[reportArgumentType]
            y_err=y_err_row,
            x_err=x_err_row,
            text_pos=(text_pos[:ntext] if i == 0 else text_pos[ntext:])
            if text_pos is not None
            else None,
            width=width[i],
            fontsize=fontsize[i],
            xlabel=xlabel[:ncol] if i == 0 else xlabel[ncol:],
            ylabel=ylabel[:ncol] if i == 0 else ylabel[ncol:],
            title=title[:ncol] if i == 0 else title[ncol:],
            color=color[:nline] if i == 0 else color[nline:],
            marker=marker[:nline] if i == 0 else marker[nline:],
            linestyle=linestyle[:nline] if i == 0 else linestyle[nline:],
            label=label[:nline] if i == 0 else label[nline:],
            tick_in=tick_in,
        )

    return axes


def plot_hist_2c(
    data: List[NDArray],
    bins: Union[int, Sequence[float], str, None] = None,
    cut: int = 1,
    axes: Optional[List[Axes]] = None,
    sharey: bool = True,
    label: Union[str, List[Union[str, None]], None] = None,
    color: Union[str, List[str], None] = None,
    density: Union[bool, List[bool]] = True,
    histtype: Union[str, List[str]] = "step",
    title: Union[str, List[Union[str, None]], None] = None,
    xlabel: Union[str, List[Union[str, None]], None] = None,
    ylabel: Union[str, List[Union[str, None]], None] = ["PDF", None],
    **kwargs,
) -> List[Axes]:
    """
    plot the histogram in two column subplots, defualt is step density histogram
    """

    if axes is None:
        fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=sharey)
    if axes is None:
        raise ValueError("ax must be specified.")
    try:
        data_left, data_right = data[:cut], data[cut:]
    except IndexError as ie:
        raise ValueError("cut index is out of the range.") from ie

    n = len(data)
    xlabel = arg_list(xlabel, 2)
    ylabel = arg_list(ylabel, 2)
    title = arg_list(title, 2)
    bins = arg_list(bins, n)
    label = arg_list(label, n)
    color = arg_list(color, n)
    density = arg_list(density, n)
    histtype = arg_list(histtype, n)

    for i, ax in enumerate(list(axes)):  # type: ignore
        b = bins[:cut] if i == 0 else bins[cut:]
        c = color[:cut] if i == 0 else color[cut:]
        d = density[:cut] if i == 0 else density[cut:]
        ht = histtype[:cut] if i == 0 else histtype[cut:]
        lb = label[:cut] if i == 0 else label[cut:]

        plot_hist(
            data_left if i == 0 else data_right,
            bins=b,
            ax=ax,
            label=lb,
            color=c,
            density=d,
            histtype=ht,
            title=title[i],
            xlabel=xlabel[i],
            ylabel=ylabel[i],
            **kwargs,
        )

    if sharey:
        fig = axes[0].get_figure()
        if fig is not None:
            fig.subplots_adjust(wspace=0)

    return axes


def plot_hist_result(
    data: List[NDArray],
    bins: Union[int, Sequence[float], str, None] = None,
    cut: int = 4,
    axes: Optional[List[Axes]] = None,
    sharey: bool = True,
    label: Union[str, List[Union[str, None]], None] = [None] * 4 + ["C", "L", "R", "B"],
    color: Union[str, List[str], None] = ["r", "b", "g", "k"] * 2,
    density: Union[bool, List[bool]] = True,
    histtype: Union[str, List[str]] = "step",
    title: Union[str, List[Union[str, None]], None] = ["HI only", "HI + noise"],
    xlabel: Union[str, List[Union[str, None]], None] = r"$T\ [\mu$K]",
    ylabel: Union[str, List[Union[str, None]], None] = ["PDF", None],
    **kwargs,
) -> List[Axes]:
    """
    plot the histogram results in two column subplots, defualt is step density histogram.
    """

    axes = plot_hist_2c(
        data,
        bins=bins,
        cut=cut,
        sharey=sharey,
        axes=axes,
        label=label,
        color=color,
        density=density,
        histtype=histtype,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        **kwargs,
    )

    for ax in axes:  # type: ignore
        ax.axvline(0, linestyle="--", c="gray", zorder=-1000)

    return axes


def compare_profiles(
    data_list: List[NDArray],
    err_list: Optional[List[NDArray]] = None,
    hinds: List[int] = [60],
    vinds: List[int] = [40, 80],
    titles: Optional[List[str]] = ["Y=0", "X=-1", "X=1"],
    axes: Optional[NDArray] = None,
    figsize: Tuple[float, float] = (12, 4),
    label_outer: bool = True,
    labels: List[str] = ["Stacked", "Fitted"],
    **kwargs,
) -> NDArray:
    """
    Compare datasets by plotting their horizontal and vertical profiles.

    This function uses `plot_line_diff` to create plots comparing horizontal
    and vertical slices of the datasets provided in `data_list`.

    Args:
        data_list (List[NDArray]): A list of numpy arrays to compare.
            Typically, this would be a list of two arrays, e.g., [data1, data2].
        err_list (Optional[List[NDArray]]): A list of numpy arrays for the errors.
        hinds (List[int]): A list of indices for horizontal profiles (Y profiles).
        vinds (List[int]): A list of indices for vertical profiles (X profiles).
        titles (Optional[List[str]]): Titles for the profile plots. The number
            of titles should match the total number of profiles (len(hinds) + len(vinds)).
        axes (Optional[NDArray]): Existing axes to plot on. Should be a
            NumPy array of shape (2, n_plots), where n_plots is the total
            number of profiles. If None, new axes are created. Defaults to None.
        figsize (Tuple[float, float]): Figure size for the entire figure if
            axes are not provided. Defaults to (12, 4).
        label_outer (bool): If True, only show outer labels and tick labels.
            Defaults to True.
        labels (List[str]): Labels for the lines being compared in each plot.
            Passed to `plot_line_diff`. Defaults to ["Stacked", "Fitted"].
        **kwargs: Additional keyword arguments to be passed to `plot_line_diff`.

    Returns:
        NDArray: The array of axes used for plotting.
    """

    def get_yprofile(dlist, ind):
        return [x[ind, :] for x in dlist]

    def get_xprofile(dlist, ind):
        return [x[:, ind] for x in dlist]

    lines_list = [get_yprofile(data_list, ind) for ind in hinds] + [
        get_xprofile(data_list, ind) for ind in vinds
    ]

    errs_list = []
    if err_list is not None:
        errs_list = [get_yprofile(err_list, ind) for ind in hinds] + [
            get_xprofile(err_list, ind) for ind in vinds
        ]

    n_plots = len(lines_list)

    if axes is None:
        _, axes = make_figure(
            nrows=2,
            ncols=n_plots,
            figsize=figsize,
            sharex=True,
            sharey="row",
            wspace=0,
            hspace=0,
            aspect=None,
            gridspec_kw={"height_ratios": [2, 1]},
        )
    if axes is None:
        raise ValueError("Failed to create axes.")

    if n_plots == 1:
        axes = axes[:, np.newaxis]

    for i, lines in enumerate(lines_list):
        title = titles[i] if titles is not None else None
        _axes = axes[:, i]

        plot_kwargs = {
            "axes": _axes,
            "title": title,
            "labels": labels,
            **kwargs,
        }

        if err_list is not None and errs_list:
            errors = errs_list[i]
            plot_kwargs["y1_err"] = errors[0]
            plot_kwargs["y2_err"] = errors[1]

        plot_line_diff(
            *lines,
            **plot_kwargs,
        )

        if label_outer:
            _axes[0].label_outer()
            _axes[1].label_outer()

    return axes


def shade_width_regions(
    ax: Axes,
    width: float = 0.5,
    color: Sequence[str] = ["pink", "lightblue", "lightblue"],
    alpha: float = 0.5,
    **kwargs,
):
    """
    Add shaded vertical regions to a matplotlib axes.

    This function creates three adjacent shaded vertical bands on the plot:
    1. A central region symmetric around zero.
    2. Two outer regions positioned further from the center.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object to which the shaded regions will be added.
    width : float, default=0.5
        Base width parameter that defines the size of the regions.
        The central region spans from `-width` to `width`.
        Each outer region has the same width as the central region
        and is positioned three widths away from the center.
    color : Sequence[str], default=["pink", "lightblue", "lightblue"]
        A sequence of three colors for the shaded regions.
        The first color is used for the central region.
        The second and third colors are used for the left and right outer regions respectively.
    alpha : float, default=0.5
        Transparency level of the shaded regions (0=transparent, 1=opaque).
    **kwargs : dict
        Additional keyword arguments passed to `ax.axvspan`.

    Notes
    -----
    The regions are defined as:
        - Central: from `-width` to `width`
        - Left outer: from `-width*4` to `-width*3`
        - Right outer: from `width*3` to `width*4`

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> shade_width_regions(ax, width=0.3, color=["yellow", "gray", "gray"])
    >>> ax.set_xlim(-2, 2)
    >>> plt.show()
    """
    ax.axvspan(-width, width, color=color[0], alpha=alpha, **kwargs)
    ax.axvspan(-width * 4, -width * 3, color=color[1], alpha=alpha, **kwargs)
    ax.axvspan(width * 3, width * 4, color=color[2], alpha=alpha, **kwargs)


def plot_profile_fixwidth_2r(
    x: List[List[Union[NDArray, Sequence[float]]]],
    y: List[List[Union[NDArray, Sequence[float]]]],
    cut: int = 1,
    axes: Union[Tuple[Axes], List[Axes], None] = None,
    y_err: Optional[List[List[Optional[Union[NDArray, Sequence[float]]]]]] = None,
    x_err: Optional[List[List[Optional[Union[NDArray, Sequence[float]]]]]] = None,
    text_pos: Optional[List[List[float]]] = None,
    width: Union[float, List[float]] = 0.5,
    fontsize: Union[float, List[float]] = 12,
    title: Union[str, List[str], None] = [
        "Transverse-section profile",
        "Lengthwise-section profile",
    ],
    xlabel: Union[str, List[Union[str, None]], None] = ["Y", "X"],
    ylabel: Union[str, List[Union[str, None]], None] = r"$T\ [\mu$K]",
    color: Union[str, List[str], None] = ["b", "r", "k"],
    marker: Union[str, List[str], None] = [".", "o", "o"],
    linestyle: Union[str, List[str]] = "--",
    label: Union[str, List[str], None] = [
        r"$|{\rm X}|< 0.5$",
        r"$T_{\rm f}$",
        r"$T_{\rm bg}$",
    ],
    tick_in: bool = True,
    add_shade: bool = True,
) -> Union[List[Axes], Tuple[Axes]]:
    """
    Plot profiles in a (2-row, n-col) subplots.

    Args:
        x,y (List[List[NDArray]]): x, y data points.
        cut (int): cut index to separate the data into upper and bottom two parts. Defaults to 1. \
            i.e. The upper and bottom panels of the first subplot is showing x[0, :cut] and x[0, cut:] respectively.
        axes (List[Axes], Tuple[Axes], None): axes object. Defaults to None.
        y_err, x_err (list, optional): errors of the data. Defaults to None.
        text_pos (list, optional): text positions. Defaults to None, to use default relative positions.
        width (float, list): width of the filament to show on the upper pannel. Defaults to 0.5.
        fontsize (float): fontsize of the text. Defaults to 12.
        title (str, list): titles of the subplots. Defaults to None.
        xlabel_up (str, list): x labels of the upper pannels. Defaults to "Y".
        xlabel_down (str, list): x labels of the lower pannels. Defaults to "X". Set to None to same as xlabel_up.
        ylabel_up (str, list): y labels of the upper pannels. Defaults to r"T[$\\mu$K]".
        ylabel_down (str, list): y labels of the lower pannels. Defaults to None to same as ylabel_up.
        color (str, list): colors of the lines. Defaults to None.
        marker (str, list): markers of the lines. Defaults to None.
        linestyle (str, list): linestyles of the lines. Defaults to '--'.
        label (str, list): labels of the lines. Defaults to None.
        tick_in (bool): whether to show the ticks inside the subplots. Defaults to True.
        add_shade (bool):  whether to add shaded width regions to the upper subplots. Defaults to True.
    """

    ncol = len(x)
    nline = len(x[0])
    nrow = 2
    ntext = 3

    if axes is None:
        _, axes = make_figure(
            2,
            ncol,
            figsize=(4 * ncol, 6),
            sharey="col",
            sharex="row",
            wspace=0,
            hspace=0.4,
            aspect=None,
        )
    if axes is None:
        raise ValueError("axes is None.")

    ## chenck line parameters
    color = arg_list(color, nline, ncol)
    marker = arg_list(marker, nline, ncol)
    linestyle = arg_list(linestyle, nline, ncol)
    label = arg_list(label, nline, ncol)

    # check text parameters
    if text_pos is not None:
        text_pos = arg_list(text_pos, ntext, ncol)
    width = arg_list(width, ncol)

    # check xlabels
    xlabel = arg_list(xlabel, nrow * ncol)
    ylabel = arg_list(ylabel, nrow * ncol)
    title = arg_list(title, nrow * ncol)

    for i in range(ncol):
        st, ed = i * nline, (i + 1) * nline
        try:
            x_up, x_down = x[i][:cut], x[i][cut:]
            y_up, y_down = y[i][:cut], y[i][cut:]
            c_up, c_down = color[st:ed][:cut], color[st:ed][cut:]
            m_up, m_down = marker[st:ed][:cut], marker[st:ed][cut:]
            ls_up, ls_down = linestyle[st:ed][:cut], linestyle[st:ed][cut:]
            la_up, la_down = label[st:ed][:cut], label[st:ed][cut:]

            y_err_up, y_err_down = None, None
            if y_err and y_err[i]:
                y_err_up, y_err_down = y_err[i][:cut], y_err[i][cut:]

            x_err_up, x_err_down = None, None
            if x_err and x_err[i]:
                x_err_up, x_err_down = x_err[i][:cut], x_err[i][cut:]

            _text_pos = (
                text_pos[i * ntext : (i + 1) * ntext] if text_pos is not None else None
            )
            _xlabel = xlabel[i * nrow : (i + 1) * nrow]
            _ylabel = ylabel[i * nrow : (i + 1) * nrow]
            _title = title[i * nrow : (i + 1) * nrow]

        except IndexError as ie:
            raise ValueError("cut index is out of the range.") from ie

        # plot upper pannel
        ax_up = axes[0, i] if ncol > 1 else axes[0]  # pyright: ignore[reportArgumentType, reportCallIssue]
        plot_line(
            x_up,
            y_up,
            ax_up,
            y_err=y_err_up,
            x_err=x_err_up,
            xlabel=_xlabel[0],
            ylabel=_ylabel[0],
            title=_title[0],
            color=c_up,
            marker=m_up,
            linestyle=ls_up,
            label=la_up,
            tick_in=tick_in,
        )
        ax_up.axhline(0, linestyle="--", c="gray")
        # add text that marker the width of filament
        strw = r"$w_{\rm t} = %.2f$" % width[i]
        if _text_pos is None:
            ax_up.text(
                0.05,
                0.9,
                strw,
                color="r",
                fontsize=fontsize,
                transform=ax_up.transAxes,
            )
        else:
            ax_up.text(*_text_pos[2], s=strw, color="r", fontsize=fontsize)
        if add_shade:
            shade_width_regions(ax_up, width=width[i])

        ax_down = axes[1, i] if ncol > 1 else axes[1]  # pyright: ignore[reportGeneralTypeIssues, reportArgumentType, reportCallIssue]
        plot_line(
            x_down,
            y_down,
            ax_down,
            y_err=y_err_down,
            x_err=x_err_down,
            xlabel=_xlabel[1],
            ylabel=_ylabel[1],
            title=_title[1],
            color=c_down,
            marker=m_down,
            linestyle=ls_down,
            label=la_down,
            tick_in=tick_in,
        )
        ax_down.axhline(0, linestyle="--", c="gray")

        # filament
        mf = np.mean(y_down[0])
        if y_err_down and y_err_down[0] is not None:
            mf_err = np.sqrt(np.sum(y_err_down[0] ** 2)) / len(y_down[0])  # pyright: ignore[reportOperatorIssue]
        else:
            mf_err = np.std(y_down[0]) / np.sqrt(len(y_down[0]))

        strf = r"$<T_{\rm f}>$ = %.2f $\pm$ %.2f $\mu$K" % (mf, mf_err)
        if _text_pos is None:
            ax_down.text(
                0.05,
                0.9,
                strf,
                color="r",
                fontsize=fontsize,
                transform=ax_down.transAxes,
            )
        else:
            ax_down.text(*_text_pos[0], s=strf, color="r", fontsize=fontsize)

        # background
        mb = np.mean(y_down[1])
        if y_err_down and y_err_down[1] is not None:
            mb_err = np.sqrt(np.sum(y_err_down[1] ** 2)) / len(y_down[1])  # pyright: ignore[reportOperatorIssue]
        else:
            mb_err = np.std(y_down[1]) / np.sqrt(len(y_down[1]))

        strb = r"$<T_{\rm bg}>$ = %.2f $\pm$ %.2f $\mu$K" % (mb, mb_err)
        if _text_pos is None:
            ax_down.text(
                0.05,
                0.8,
                strb,
                color="k",
                fontsize=fontsize,
                transform=ax_down.transAxes,
            )
        else:
            ax_down.text(*_text_pos[1], s=strb, color="k", fontsize=fontsize)

    return axes


def plot_profile_fixwidth_2c(
    x: List[Union[NDArray, Sequence[float]]],
    y: List[Union[NDArray, Sequence[float]]],
    cut: int = 1,
    axes: Optional[List[Axes]] = None,
    y_err: Optional[List[Optional[Union[NDArray, Sequence[float]]]]] = None,
    x_err: Optional[List[Optional[Union[NDArray, Sequence[float]]]]] = None,
    text_pos: Optional[List[List[float]]] = None,
    width: float = 0.5,
    fontsize: Union[float, List[float]] = 12,
    xlabel: Union[str, List[Union[str, None]], None] = ["Y", "X"],
    ylabel: Union[str, List[Union[str, None]], None] = [r"$T\ [\mu$K]", None],
    title: Union[str, List[Union[str, None]], None] = [
        "Transverse-section profile",
        "Lengthwise-section profile",
    ],
    color: Union[str, List[str], None] = ["b", "r", "k"],
    marker: Union[str, List[str], None] = [".", "o", "o"],
    linestyle: Union[str, List[str], None] = "--",
    label: Union[str, List[Union[str, None]], None] = [
        r"$|{\rm X}|< 0.5$",
        r"$T_{\rm f}$",
        r"$T_{\rm bg}$",
    ],
    tick_in: bool = True,
    add_shade: bool = True,
) -> List[Axes]:
    """
    plot profiles in two column subplots.

    Parameters:
        x (list): x data points.
        y (list): y data points.
        cut (int): separate the data points into two parts, \
            e.g. cut=2 means the data points in the first two rows are plotted in the first subplot, \
                and rest are plotted in the second subplot.
        axes (Axes, optional): axes object. Defaults to None.
        y_err (list, optional): y error of the data. Defaults to None.
        x_err (list, optional): x error of the data. Defaults to None.
        text_pos (List[list], optional): text position. Defaults to None, and the text position will be auto.
        width (float): width of the filament to show on the left panel. Defaults to 0.5.
        fontsize (str): fontsize of the text. Defaults to 12.
        xlabel (str, list): x-axis labels. Defaults to ["Y", "X"].
        ylabel (str, list): y-axis labels. Defaults to [r"$T$ [$\\mu$K]", None].
        title (str, list): subplot titles. Defaults to ["Transverse-section profile", \
            "Lengthwise-section profile"].
        color (str, list): line colors. Defaults to ["b", "r", "k"].
        marker (str, list): line markers. Defaults to [".", "o", "o"].
        linestyle (str, list): line styles. Defaults to "--".
        label (str, list): line labels. Defaults to [r"$T_{\rm f}$", r"$T_{\rm bg$"].
        tick_in (bool): whether to show ticks inside the axes. Defaults to True.
        add_shade (bool):  whether to add shaded width regions to the upper subplots. Defaults to True.
    Returns:
        axes (list): axes object.

    """
    if axes is None:
        _, axes = make_figure(
            1,
            2,
            figsize=(8, 3),
            sharey=True,
            sharex=False,
            wspace=0,
            aspect=None,
        )
    if axes is None:
        raise ValueError("axes is None.")

    try:
        x_left, x_right = x[:cut], x[cut:]
        y_left, y_right = y[:cut], y[cut:]
    except IndexError as ie:
        raise ValueError("cut index is out of the range.") from ie

    y_err_left, y_err_right = None, None
    if y_err:
        y_err_left, y_err_right = y_err[:cut], y_err[cut:]

    x_err_left, x_err_right = None, None
    if x_err:
        x_err_left, x_err_right = x_err[:cut], x_err[cut:]

    nline = len(x)
    xlabel = arg_list(xlabel, 2)
    ylabel = arg_list(ylabel, 2)
    title = arg_list(title, 2)

    color = arg_list(color, nline)
    marker = arg_list(marker, nline)
    linestyle = arg_list(linestyle, nline)
    label = arg_list(label, nline)

    for i, ax in enumerate(list(axes)):  # type: ignore
        c = color[:cut] if i == 0 else color[cut:]
        m = marker[:cut] if i == 0 else marker[cut:]
        ls = linestyle[:cut] if i == 0 else linestyle[cut:]
        lb = label[:cut] if i == 0 else label[cut:]
        plot_line(
            x_left if i == 0 else x_right,
            y_left if i == 0 else y_right,
            ax=ax,
            y_err=y_err_left if i == 0 else y_err_right,
            x_err=x_err_left if i == 0 else x_err_right,
            xlabel=xlabel[i],
            ylabel=ylabel[i],
            title=title[i],
            color=c,
            marker=m,
            linestyle=ls,
            label=lb,
            tick_in=tick_in,
        )
        # add horizontal lines that marker T = 0
        ax.axhline(0, linestyle="--", c="gray")
        # add text that marker the mean values of filament and background.
        if i == 1:
            # filament
            mf = np.mean(y_right[0])
            if y_err_right and y_err_right[0] is not None:
                mf_err = np.sqrt(np.sum(y_err_right[0] ** 2)) / len(y_right[0])  # pyright: ignore[reportOperatorIssue]
            else:
                mf_err = np.std(y_right[0]) / np.sqrt(len(y_right[0]))

            strf = r"$<T_{\rm f}>$ = %.2f $\pm$ %.2f $\mu$K" % (mf, mf_err)
            if text_pos is None:
                ax.text(
                    0.05,
                    0.9,
                    strf,
                    color="r",
                    fontsize=fontsize,
                    transform=ax.transAxes,
                )
            else:
                ax.text(*text_pos[0], s=strf, color="r", fontsize=fontsize)

            # background
            mb = np.mean(y_right[1])
            if y_err_right and y_err_right[1] is not None:
                mb_err = np.sqrt(np.sum(y_err_right[1] ** 2)) / len(y_right[1])  # pyright: ignore[reportOperatorIssue]
            else:
                mb_err = np.std(y_right[1]) / np.sqrt(len(y_right[1]))

            strb = r"$<T_{\rm bg}>$ = %.2f $\pm$ %.2f $\mu$K" % (mb, mb_err)
            if text_pos is None:
                ax.text(
                    0.05,
                    0.8,
                    strb,
                    color="k",
                    fontsize=fontsize,
                    transform=ax.transAxes,
                )
            else:
                ax.text(*text_pos[1], s=strb, color="k", fontsize=fontsize)

        # add text that marker the width of filament
        else:
            strw = r"$w_{\rm t} = %.2f$" % width
            if text_pos is None:
                ax.text(
                    0.05,
                    0.9,
                    strw,
                    color="r",
                    fontsize=fontsize,
                    transform=ax.transAxes,
                )
            else:
                ax.text(
                    *text_pos[2],
                    s=strw,
                    color="r",
                    fontsize=fontsize,
                )
            # add shades vertical span to show width regions
            if add_shade:
                shade_width_regions(ax, width=width)
    return axes


def plot_pair_sep(
    data,
    bins=20,
    ax=None,
    label=["Transversal", "Radial", "Spatial"],
    color=["r", "b", "k"],
    linewidth=1.4,
    xlabel=r"Separation [$h^{-1}$ Mpc]",
    show_mean: bool = True,
    **kwargs,
):
    """
    Plot the separation distribution of selected galaxy pairs.

    Parameters
    ----------
    data : Sequence[np.ndarray]
        A sequence of data arrays representing pair separations. For example,
        this could be a list of [transversal, radial, spatial] separations.
    bins : int or sequence, default=20
        The number of bins for the histogram or a sequence of bin edges.
    ax : matplotlib.axes.Axes, optional
        The axes object to plot on. If None, a new one is created by `plot_hist`.
    label : list of str, default=["Transversal", "Radial", "Spatial"]
        Labels for each dataset in `data`.
    color : list of str, default=["r", "b", "k"]
        Colors for each histogram.
    linewidth : float, default=1.4
        The width of the histogram lines.
    xlabel : str, default=r"Separation [$h^{-1}$ Mpc]"
        The label for the x-axis.
    show_mean : bool, default=True
        If True, a vertical dashed line is drawn at the mean of each dataset.
    **kwargs
        Additional keyword arguments passed to `plot_hist`.

    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the plot.
    """

    ax = plot_hist(
        data,
        bins=bins,
        ax=ax,
        label=label,
        color=color,
        linewidth=linewidth,
        xlabel=xlabel,
        **kwargs,
    )

    if show_mean:
        for d, c, la in zip(data, color, label):
            m = d.mean()
            ax.axvline(m, linestyle="--", c=c, zorder=-100)
            print("%s mean separation: %.2f Mpc/h" % (la, m))

    return ax
