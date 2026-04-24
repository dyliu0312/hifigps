"""
General and useful functions for plotting.
"""

from typing import (
    Any,
    Iterable,
    List,
    Literal,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from matplotlib.figure import Figure
from matplotlib.patches import Arc, Ellipse
from numpy.typing import NDArray


def arg_list(arg: Any, n: int = 1, m: Optional[int] = None) -> List[Any]:
    """
    check the arg list.
    If arg is not a list, return a list with arg repeated n times.
    If arg is a list, return a list with length n, if the length of arg is less than n, append the last element of arg to the list.

    Args:
        arg (Any): the arg to check.
        n (int, optional): the length of the list. Defaults to 1.
        m (Optional[int], optional): the number of times to repeat the list. Defaults to None.

    Returns:
        List[Any]: the list of args.
    """
    if not isinstance(arg, list):
        arg = [arg] * n
    if isinstance(arg, list) and len(arg) < n:
        arg = arg + [arg[-1]] * (n - len(arg))

    if m is not None and m > 1:
        arg = arg * m

    return arg


def get_minmax(data: NDArray, q: float = 5) -> NDArray:
    """
    get the min and max values of data with percentile.
    """
    if isinstance(data, np.ma.MaskedArray):
        data = data.compressed()

    return np.percentile(data, [q, 100 - q])


def get_ticks_labels(
    ticks: List[int] = [20, 60, 100],
    shape: int = 121,
    lim: Union[Tuple[float, float], List[float]] = (-3, 3),
    digits: int = 0,
) -> Tuple[List[int], List[str]]:
    """
    get the ticks and tick labels for a range.

    Args:
        ticks (List[int], optional): the number of ticks. Defaults to [20, 60, 100].
        shape (int, optional): the shape of the ticks. Defaults to 120.
        lim (Tuple[float, float]): the range of ticks. Defaults to (-3, 3).
        digits (int, optional): the number of digits of tick labels. Defaults to 0.
    Returns:
        Tuple(List[int], List[str]): the ticks and tick labels.
    """
    tick_labels = np.linspace(lim[0], lim[1], shape)
    tick_labels = [f"%.{digits}f" % tick_labels[t] for t in ticks]
    return ticks, tick_labels


def set_ticks(
    ax: Axes,
    xticks: Optional[List[int]] = None,
    xticklabels: Optional[List[str]] = None,
    yticks: Optional[List[int]] = None,
    yticklabels: Optional[List[str]] = None,
    share: bool = True,
    kw_get_ticks: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
    **kargs,
) -> None:
    """
    set the ticks and ticklabels for x and y axis.

    Args:
        ax (Axes): the axes to set ticks.
        xticks (List[int], optional): the ticks for x axis. Defaults to None.
        xticklabels (List[str], optional): the tick labels for x axis. Defaults to None.
        yticks (List[int], optional): the ticks for y axis. Defaults to None.
        yticklabels (List[str], optional): the tick labels for y axis. Defaults to None.
        share (bool, optional): whether to share the ticks with x axis. Defaults to True.
        kw_get_ticks (dict, optional): the kwargs for `get_ticks_labels`. Defaults to None. The default value is `dict(ticks=[20, 60, 100], shape=121, lim=(-3, 3), digits=0)`.
        **kargs (dict, optional): the kwargs for `Axes.set_xticks` and `Axes.set_yticks`. Defaults to None.
    """
    if xticks is None and yticks is None:
        default_kwargs = dict(ticks=[20, 60, 100], shape=121, lim=(-3, 3), digits=0)
        if kw_get_ticks is not None:
            default_kwargs.update(kw_get_ticks)
        xticks, xticklabels = get_ticks_labels(**default_kwargs)  # pyright: ignore[reportArgumentType]

    if xticks is not None:
        ax.set_xticks(xticks, xticklabels, **kargs)
        if share and yticks is None:
            ax.set_yticks(xticks, xticklabels, **kargs)
    if yticks is not None:
        ax.set_yticks(yticks, yticklabels, **kargs)
        if share and xticks is None:
            ax.set_xticks(yticks, yticklabels, **kargs)


def get_colorbar_cax(
    ax: Axes,
    loc: Literal["right", "top", "bottom"] = "right",
    w_pad: float = 0.01,
    w_factor: float = 0.05,
    h_pad: float = 0.01,
    h_factor: float = 0.05,
    **kwargs,
) -> Axes:
    """
    get colorbar axes.
    """
    # Get the bounding box of the axes
    bbox = ax.get_position()

    # Calculate the position of the colorbar
    if loc == "right":
        cax_width = bbox.width * w_factor
        cax_height = bbox.height
        cax_x0 = bbox.x1 + w_pad
        cax_y0 = bbox.y0
    elif loc == "top":
        cax_width = bbox.width
        cax_height = bbox.height * h_factor
        cax_x0 = bbox.x0
        cax_y0 = bbox.y1 + h_pad
    elif loc == "bottom":
        cax_width = bbox.width
        cax_height = bbox.height * h_factor
        cax_x0 = bbox.x0
        cax_y0 = bbox.y0 - h_pad - cax_height

    else:
        raise ValueError(f"colorbar location = {loc} is not supported")

    # Create the colorbar axes
    fig = ax.get_figure()

    if fig is None:
        raise ValueError("ax must be a matplotlib.axes.Axes instance")

    cax = fig.add_axes((cax_x0, cax_y0, cax_width, cax_height), **kwargs)
    return cax


def make_figure(
    nrows: int = 1,
    ncols: int = 1,
    figsize: Tuple[float, float] = (4, 3),
    sharex: Union[bool, Literal["none", "all", "row", "col"]] = True,
    sharey: Union[bool, Literal["none", "all", "row", "col"]] = True,
    wspace: float = 0.0,
    hspace: float = 0.0,
    aspect: Optional[str] = "equal",
    gridspec_kw: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
    subplot_kw: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
    **kwargs,
) -> Tuple[Figure, Any]:
    """
    make a figure and axes with subplots.

    Args:
        nrows (int, optional): number of rows. Defaults to 1.
        ncols (int, optional): number of columns. Defaults to 1.
        figsize (Tuple[float, float], optional): figure size. Defaults to (6, 4).
        sharex (bool, optional): share x axis or not. Defaults to True.
        sharey (bool, optional): share y axis or not. Defaults to True.
        wspace (float, optional): width space between subplots. Defaults to 0.
        hspace (float, optional): height space between subplots. Defaults to 0.
        aspect (str, optional): aspect of subplots. Defaults to "equal". Set to None for auto aspect.
        gridspec_kw (dict, optional): gridspec kwargs. Defaults to None.
        subplot_kw (dict, optional): subplot kwargs. Defaults to None.
    Returns:
        figure and axes.
    """

    g_kw = {"wspace": wspace, "hspace": hspace}
    if gridspec_kw is not None:
        g_kw.update(gridspec_kw)

    if aspect is not None:
        s_kw = {"aspect": aspect}
    else:
        s_kw = {}
    if subplot_kw is not None:
        s_kw.update(subplot_kw)

    return plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=figsize,
        sharex=sharex,
        sharey=sharey,
        subplot_kw=s_kw,
        gridspec_kw=g_kw,
        **kwargs,
    )


def get_norm(
    norm: Literal["linear", "log", "symlog"] = "linear",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
):
    """
    get norm for colorbar.

    Args:
        norm (str, optional): norm type. Defaults to "linear".
        vmin (float, optional): vmin of colorbar. Defaults to None.
        vmax (float, optional): vmax of colorbar. Defaults to None.
    Returns:
        norm.
    """
    if norm == "linear":
        return Normalize(vmin=vmin, vmax=vmax)
    elif norm == "log":
        return LogNorm(vmin=vmin, vmax=vmax)
    elif norm == "symlog":
        return SymLogNorm(linthresh=1e-3, vmin=vmin, vmax=vmax)


def plot_heatmap(
    *args,
    ax: Optional[Axes] = None,
    cmap: str = "viridis",
    norm: Union[Literal["linear", "log", "symlog"], Normalize] = "linear",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    q: Optional[float] = 5,
    tick_in: bool = True,
    show_cbar: bool = True,
    cbar_ax: Optional[Axes] = None,
    cbar_loc: Literal["right", "top", "bottom"] = "right",
    cbar_label: Optional[str] = r"$T\ [\mu$K]",
    title: Optional[str] = "Heatmap",
    xlabel: Optional[str] = "X",
    ylabel: Optional[str] = "Y",
    change_ticks: bool = True,
    kw_makefigure: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
    kw_pcolormesh: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
    kw_cbar: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
    kw_ticks: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
    sym_cbar: bool = False,
) -> Axes:
    """
    plot heatmap with colorbar.

    Args:
        x (array-like, optional): x axis. Defaults is None.
        y (array-like, optional): y axis. Defaults is None.
        data (array-like): data to plot.
        ax (matplotlib.axes.Axes, optional): axis to plot. Defaults is None, create a new axis with `make_figure`.
        cmap (str, optional): color map. Defaults is 'viridis'.
        norm (str, optional): colorbar norm. Defaults is None, use 'linear'.
        vmin (float, optional): vmin of colorbar. Defaults is None, to be calculated with `get_minmax`.
        vmax (float, optional): vmax of colorbar. Defaults is None, to be calculated with `get_minmax`.
        q (float, optional): percentage of data to be used to calculate vmin and vmax. Defaults is 5.
        tick_in (bool, optional): whether to turn the direction of ticks to in. Defaults is True.
        show_cbar (bool, optional): whether to show colorbar. Defaults is True.
        cbar_ax (matplotlib.axes.Axes, optional): colorbar axis. Defaults is None, create a new axis with `get_colorbar_cax`.
        cbar_loc (str, optional): colorbar location. Defaults is 'right'. Options are 'right', 'top', 'bottom'.
        cbar_label (str, optional): colorbar label. Defaults is 'T[$\\mu$K]'.
        title (str, optional): plot title. Defaults is None.
        xlabel (str, optional): x axis label. Defaults is 'X'.
        ylabel (str, optional): y axis label. Defaults is 'Y'.
        change_ticks (bool, optional): whether to change ticks. Defaults is True.
        kw_makefigure (dict, optional): kwargs for make_figure. Defaults is None.
        kw_pcolormesh (dict, optional): kwargs for pcolormesh. Defaults is None.
        kw_cbar (dict, optional): kwargs for colorbar. Defaults is None.
        kw_ticks (dict, optional): kwargs for `set_ticks`. Defaults is None to use default values `dict(share=True, kw_get_ticks = dict(ticks=[20, 60, 100], shape=121, lim=(-3, 3), digits=0))`.
        sym_cbar: bool, optional: whether to use symmetric colorbar limits based on the maximum value between vmax and vmin. Defaults is False.
    Returns:
        matplotlib.axes.Axes: axis with heatmap and colorbar.
    """
    if ax is None:
        kw_fig = {"figsize": (4, 3)}
        if kw_makefigure is not None:
            kw_fig.update(kw_makefigure)
        fig, ax = make_figure(**kw_fig)  # pyright: ignore[reportArgumentType]
    else:
        fig = ax.get_figure()

    if ax is None:
        raise ValueError("ax is None")
    if fig is None:
        raise ValueError("fig is None")

    n = len(args)

    if n == 1:
        data = args[0]
    elif n == 2:
        x, data = args
        y = x
    elif n == 3:
        x, y, data = args
    else:
        raise ValueError("Expected 1-3 arguments: [data], [x, data], or [x, y, data].")

    # if data is boolean array
    if data.dtype.kind == "b":
        q = None  # disable percentile calculation

    # --- Compute vmin/vmax ---
    if q is not None:
        _vmin, _vmax = get_minmax(data, q=q)
        vmin = _vmin if vmin is None else vmin
        vmax = _vmax if vmax is None else vmax
        if sym_cbar:
            max_abs = max(abs(vmin), abs(vmax))  # pyright: ignore[reportArgumentType]
            vmin = -max_abs
            vmax = max_abs

    # --- Set norm ---
    if isinstance(norm, str):
        norm = get_norm(norm, vmin=vmin, vmax=vmax)
        p_kw = {"cmap": cmap, "norm": norm}
    else:
        p_kw = {"cmap": cmap, "norm": norm, "vmin": vmin, "vmax": vmax}

    if kw_pcolormesh is not None:
        p_kw.update(kw_pcolormesh)

    if n == 1:
        pcm = ax.pcolormesh(data, **p_kw)  # pyright: ignore[reportArgumentType]
    else:
        pcm = ax.pcolormesh(x, y, data, **p_kw)  # pyright: ignore[reportPossiblyUnboundVariable, reportArgumentType]

    if tick_in:
        ax.tick_params(axis="both", direction="in")

    if title is not None:
        ax.set_title(title)  # type: ignore
    if xlabel is not None:
        ax.set_xlabel(xlabel)  # type: ignore
    if ylabel is not None:
        ax.set_ylabel(ylabel)  # type: ignore

    if show_cbar:
        if cbar_ax is None:
            cbar_ax = get_colorbar_cax(ax, loc=cbar_loc)  # type: ignore
        if cbar_ax is None:
            raise ValueError("cbar_ax is None")

        orientation = "horizontal" if cbar_loc in ["bottom", "top"] else "vertical"
        if cbar_label is None:
            cbar_label = ""
        cb_kw = {"label": cbar_label, "orientation": orientation}
        if kw_cbar is not None:
            cb_kw.update(kw_cbar)
        cb = fig.colorbar(pcm, cax=cbar_ax, **cb_kw)  # pyright: ignore[reportArgumentType]

        if cbar_loc == "top":
            cb.ax.xaxis.set_label_position("top")
            cb.ax.xaxis.set_ticks_position("top")
        if cbar_loc == "bottom":
            cb.ax.xaxis.set_label_position("bottom")
            cb.ax.xaxis.set_ticks_position("bottom")
        if tick_in:
            cb.ax.tick_params(axis="both", direction="in")

    if change_ticks:
        kw_set_ticks = dict(
            share=True,
            kw_get_ticks=dict(ticks=[20, 60, 100], shape=121, lim=(-3, 3), digits=0),
        )
        if kw_ticks is not None:
            kw_set_ticks.update(kw_ticks)
        set_ticks(ax, **kw_set_ticks)  # pyright: ignore[reportArgumentType]
    return ax


def plot_heatmaps(
    data: List[NDArray],
    axes: Optional[List[Axes]] = None,
    cmap: Union[str, List[str]] = "viridis",
    norm: Union[str, List[str]] = "linear",
    vmin: Union[float, List[float], None] = None,
    vmax: Union[float, List[float], None] = None,
    q: float = 5,
    tick_in: bool = True,
    show_cbar: Union[bool, List[bool]] = True,
    cbar_loc: Literal["right", "top", "bottom"] = "right",
    cbar_label: Union[str, List[str], None] = r"$T\ [\mu$K]",
    title: Union[str, List[str], None] = None,
    xlabel: Union[str, List[str]] = "X",
    ylabel: Union[str, List[str]] = "Y",
    change_ticks: bool = True,
    label_outer: bool = True,
    kw_makefigure: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
    kw_pcolormesh: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
    kw_cbar: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
    kw_ticks: Optional[dict] = None,  # pyright: ignore[reportMissingTypeArgument]
    sym_cbar: Union[bool, List[bool]] = False,
) -> List[Axes]:
    """
    plot multiple heatmaps in a row.

    Notes:
        All args (except axes) are list, and the length of the list is the number of heatmaps to plot.
        If an arg is not a list, it will be converted to a list with arg be repeated.

    Args:
        data: list of NDArray, the data to plot.
        axes: list of matplotlib.axes.Axes.
        cmap: str or a list of str, the colormap.
        norm: matplotlib.color.norm or str or a list of norm.
        vmin: float or a list of float, the minimum value of the colormap.
        vmax: float or a list of float, the maximum value of the colormap.
        q: float, the percentage used to calculate the minimum and maximum value of the colormap.
        tick_in: bool, whether to show ticks inside the plot.
        show_cbar: bool or a list of bool, whether to show colorbar.
        cbar_loc: str, the location of colorbar. Defaults is "right".
        cbar_label: str or a list of str, the label of colorbar. Defaults is "T [$\\mu$K]".
        title: str or a list of str, the title of the plot.
        xlabel: str or a list of str, the label of x axis.
        ylabel: str or a list of str, the label of y axis.
        change_ticks: bool, whether to change the ticks of the plot.
        label_outer: bool, whether to label the outer axes.
        kw_makefigure: dict, kwargs for make_figure.
        kw_pcolormesh: dict, kwargs for pcolormesh.
        kw_cbar: dict, kwargs for colorbar.
        kw_ticks: dict, kwargs for `set_ticks`.
        sym_cbar: bool or a list of bool, whether to use symmetric colorbar limits. Defaults is False.
    Returns:
        axes
    """
    n = len(data)
    cmap = arg_list(cmap, n)
    norm = arg_list(norm, n)
    vmin = arg_list(vmin, n)  # type: ignore
    vmax = arg_list(vmax, n)  # type: ignore

    show_cbar = arg_list(show_cbar, n)  # type: ignore
    cbar_label = arg_list(cbar_label, n)
    sym_cbar = arg_list(sym_cbar, n)
    title = arg_list(title, n)
    xlabel = arg_list(xlabel, n)
    ylabel = arg_list(ylabel, n)

    if axes is None:
        kw_fig = {"nrows": 1, "ncols": n, "figsize": (4 * n, 3)}
        if kw_makefigure is not None:
            kw_fig.update(kw_makefigure)
        _, axes = make_figure(**kw_fig)  # pyright: ignore[reportArgumentType]

    if axes is None:
        raise ValueError("axes is None")

    for i, ax in enumerate(list(axes)):  # type: ignore
        plot_heatmap(
            data[i],
            ax=ax,
            cmap=cmap[i],
            norm=norm[i],  # pyright: ignore[reportArgumentType]
            vmin=vmin[i],  # type: ignore
            vmax=vmax[i],  # type: ignore
            q=q,
            tick_in=tick_in,
            show_cbar=show_cbar[i],  # type: ignore
            cbar_loc=cbar_loc,
            cbar_label=cbar_label[i],
            title=title[i],
            xlabel=xlabel[i],
            ylabel=ylabel[i],
            change_ticks=change_ticks,
            kw_pcolormesh=kw_pcolormesh,
            kw_cbar=kw_cbar,
            kw_ticks=kw_ticks,
            sym_cbar=sym_cbar[i],
        )

        if label_outer:
            ax.label_outer()

    return axes


def plot_line(
    x: List[Union[NDArray, Sequence[float]]],
    y: List[Union[NDArray, Sequence[float]]],
    ax: Optional[Axes] = None,
    y_err: Optional[List[Optional[Union[NDArray, Sequence[float]]]]] = None,
    x_err: Optional[List[Optional[Union[NDArray, Sequence[float]]]]] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    color: Union[str, List[str], List[Union[str, None]], List[None], None] = None,
    marker: Union[str, List[str], List[Union[str, None]], List[None], None] = None,
    linestyle: Union[str, List[str], List[Union[str, None]], List[None], None] = None,
    label: Union[str, List[str], List[Union[str, None]], List[None], None] = None,
    tick_in: bool = True,
    **kwargs,
) -> Axes:
    """
    plot a line or lines.

    Args:
        x (list): x data points.
        y (list): y data points.
        ax (matplotlib.axes.Axes, optional): axes object. Defaults is None, use the axes of new created figure.
        y_err (list, optional): y error of the data. Defaults to None.
        x_err (list, optional): x error of the data. Defaults to None.
        xlabel (str, optional): x label. Defaults is None.
        ylabel (str, optional): y label. Defaults is None.
        title (str, optional): title. Defaults is None.
        color (str, list, optional): color of the line. Defaults is None.
        marker (str, list, optional): marker of the line. Defaults is None.
        linestyle (str, list, optional): linestyle of the line. Defaults is None.
        label (str, list, optional): label of the line. Defaults is None.
        tick_in (bool, optional): whether to show ticks inside the axes. Defaults is True.
        **kwargs: other arguments passed to ax.plot().
    Returns:
        ax (matplotlib.axes.Axes): axes object.
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(5, 4))

    nline = len(x)

    colors = arg_list(color, nline)
    markers = arg_list(marker, nline)
    linestyles = arg_list(linestyle, nline)
    labels = arg_list(label, nline)
    y_errs = arg_list(y_err, nline)
    x_errs = arg_list(x_err, nline)

    # Plot lines with labels first
    for i in range(nline):
        la = labels[i]
        if la is not None:
            if (y_errs[i] is not None) or (x_errs[i] is not None):
                ax.errorbar(
                    x[i],
                    y[i],
                    yerr=y_errs[i],
                    xerr=x_errs[i],
                    color=colors[i],
                    linestyle=linestyles[i],
                    marker=markers[i],
                    label=la,
                    **kwargs,
                )
            else:
                ax.plot(
                    x[i],
                    y[i],
                    color=colors[i],
                    linestyle=linestyles[i],
                    marker=markers[i],
                    label=la,
                    **kwargs,
                )

    if any(la is not None for la in labels):
        ax.legend()

    # Plot lines without labels
    for i in range(nline):
        la = labels[i]
        if la is None:
            if (y_errs[i] is not None) or (x_errs[i] is not None):
                ax.errorbar(
                    x[i],
                    y[i],
                    yerr=y_errs[i],
                    xerr=x_errs[i],
                    color=colors[i],
                    linestyle=linestyles[i],
                    marker=markers[i],
                    **kwargs,
                )
            else:
                ax.plot(
                    x[i],
                    y[i],
                    color=colors[i],
                    linestyle=linestyles[i],
                    marker=markers[i],
                    **kwargs,
                )

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)
    if tick_in:
        ax.tick_params(axis="both", which="both", direction="in")

    return ax


def plot_hist(
    data: List[NDArray],
    bins: Union[int, Sequence[float], str, None] = None,
    ax: Optional[Axes] = None,
    label: Union[str, List[Union[str, None]], None] = None,
    color: Union[str, List[str], None] = None,
    density: Union[bool, List[bool]] = True,
    histtype: Union[str, List[str]] = "step",
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = "PDF",
    **kwargs,
) -> Axes:
    """
    plot the histogram, defualt is step density histogram
    """

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    if ax is None:
        raise ValueError("ax must be specified.")
    data_array = np.array(data)
    if len(data_array.shape) == 1:
        data_array = data_array.reshape(-1, 1)
    n, _ = data_array.shape

    bins = arg_list(bins, n)
    label = arg_list(label, n)
    color = arg_list(color, n)
    density = arg_list(density, n)
    histtype = arg_list(histtype, n)

    for i, la in enumerate(label):
        if la is not None:
            ax.hist(
                data_array[i],
                bins=bins[i],
                label=la,
                color=color[i],
                density=density[i],
                histtype=histtype[i],  # pyright: ignore[reportArgumentType]
                **kwargs,
            )
    if any([la is not None for la in label]):
        ax.legend()

    for i, la in enumerate(label):
        if la is None:
            ax.hist(
                data_array[i],
                bins=bins[i],
                color=color[i],
                density=density[i],
                histtype=histtype[i],  # pyright: ignore[reportArgumentType]
                **kwargs,
            )

    if title is not None:
        ax.set_title(title)  # type: ignore
    if xlabel is not None:
        ax.set_xlabel(xlabel)  # type: ignore
    if ylabel is not None:
        ax.set_ylabel(ylabel)  # type: ignore

    return ax


def plot_ellipses(
    ax: Axes,
    xy_pos: List[Sequence[float]] = [(0, 0)],
    width: Union[float, List[float]] = 1,
    height: Union[float, List[float], None] = None,
    angle: Union[float, List[float]] = 0,
    **kwargs,
) -> None:
    """
    Plot ellipses (or circles) on the given axis.

    Args:
        ax (matplotlib.axes.Axes): The axis on which to plot the ellipses.
        xy_pos (list[tuple(float, float)]): xy positions.
        width  (float, list): widths of the ellipses.
        height (float, list): heights of the ellipses, default is None, use the width.
        angle  (float, list): rotation of the ellipse in degrees (**counterclockwise**)
        theta1,theta2 (float, list): starting and ending angles in degrees.
        **kwargs: additional keyword arguments.
    """
    if height is None:
        height = width

    # find the number of ellipse need to plot
    n = len(xy_pos)
    if n == 1:
        e = Ellipse(xy_pos[0], width, height, angle=angle, **kwargs)  # pyright: ignore[reportArgumentType]
        ax.add_patch(e)
    else:
        width_list = arg_list(width, n)
        height_list = arg_list(height, n)
        angle_list = arg_list(angle, n)
        for key, value in kwargs.items():
            kwargs[key] = arg_list(
                value, n
            )  # Ensure any list in kwargs is the correct length

        for i, (xy, w, h, a) in enumerate(
            zip(xy_pos, width_list, height_list, angle_list)
        ):
            # Extract other properties from kwargs
            ellipse_kwargs = {key: value[i] for key, value in kwargs.items()}
            e = Ellipse(xy, w, h, angle=a, **ellipse_kwargs)  # pyright: ignore[reportArgumentType]
            ax.add_patch(e)


def plot_arcs(
    ax: Axes,
    xy_pos: List[Sequence[float]] = [(0, 0)],
    width: Union[float, List[float]] = 1,
    height: Union[float, List[float], None] = None,
    angle: Union[float, List[float]] = 0,
    theta1: Union[float, List[float]] = 0,
    theta2: Union[float, List[float]] = 360,
    **kwargs,
) -> None:
    """
    Plot arcs on the given axis.

    Args:
        ax (matplotlib.axes.Axes): The axis on which to plot the ellipses.
        xy_pos (list[tuple(float, float)]): xy positions.
        width  (float, list): widths of the ellipses.
        height (float, list): heights of the ellipses, default is None, use the width.
        angle  (float, list): rotation of the ellipse in degrees (**counterclockwise**)
        theta1,theta2 (float, list): starting and ending angles in degrees.
        **kwargs: additional keyword arguments.
    """
    if height is None:
        height = width

    # find the number of ellipse need to plot
    n = len(xy_pos)
    if n == 1:
        a = Arc(
            xy_pos[0],  # pyright: ignore[reportArgumentType]
            width,  # pyright: ignore[reportArgumentType]
            height,  # pyright: ignore[reportArgumentType]
            angle=angle,  # pyright: ignore[reportArgumentType]
            theta1=theta1,  # pyright: ignore[reportArgumentType]
            theta2=theta2,  # pyright: ignore[reportArgumentType]
            **kwargs,
        )  # type: ignore
        ax.add_patch(a)
    else:
        width_list = arg_list(width, n)
        height_list = arg_list(height, n)
        angle_list = arg_list(angle, n)
        theta1_list = arg_list(theta1, n)
        theta2_list = arg_list(theta2, n)
        for key, value in kwargs.items():
            kwargs[key] = arg_list(
                value, n
            )  # Ensure any list in kwargs is the correct length

        for i, (xy, w, h, a, t1, t2) in enumerate(
            zip(xy_pos, width_list, height_list, angle_list, theta1_list, theta2_list)
        ):
            arg_kwargs = {key: value[i] for key, value in kwargs.items()}
            a = Arc(xy, w, h, angle=a, theta1=t1, theta2=t2, **arg_kwargs)
            ax.add_patch(a)


def plot_sector(
    ax: Axes,
    xy_pos: List[Tuple[float, float]] = [(40, 60), (80, 60)],
    width: List[float] = [10, 30, 50, 70, 90, 110],
    height: Optional[List[float]] = None,
    angle: Union[float, List[float]] = 0,
    theta1: Union[float, List[float]] = [45, 225],
    theta2: Union[float, List[float]] = [135, 315],
    linestyle="--",
    ec="k",
    linewidth=1.2,
    **kwargs,
) -> None:
    n_arcs = len(width)
    n_sectors = len(xy_pos)

    for i in range(n_sectors):
        for j in range(2):  # plot two arcs area for each sector
            plot_arcs(
                ax,
                xy_pos=[xy_pos[i]] * n_arcs,
                width=width,
                height=height,
                angle=angle,
                theta1=theta1[j],  # pyright: ignore[reportIndexIssue]
                theta2=theta2[j],  # pyright: ignore[reportIndexIssue]
                linestyle=linestyle,
                ec=ec,
                linewidth=linewidth,
                **kwargs,
            )


def plot_axlines(
    ax: Union[Axes, Sequence[Axes], NDArray],
    vl: Union[float, Iterable[float], None] = None,
    hl: Union[float, Iterable[float], None] = None,
    color: str = "r",
    linestyle: str = "--",
    linewidth: float = 1,
    alpha: float = 1,
    zorder: float = 1,
    **kwargs,
) -> None:
    """
    Plot vertical and horizontal lines on the given axis.

    Args:
        ax (matplotlib.axes.Axes, list, numpy.ndarray): The axis on which to plot the lines.
        vl (float, list): x positions of vertical lines.
        hl (float, list): y positions of horizontal lines.
        color (str): color of the lines.
        linestyle (str): linestyle of the lines.
        linewidth (float): linewidth of the lines.
        alpha (float): alpha of the lines.
        zorder (float): zorder of the lines.
        **kwargs: additional keyword arguments.
    """

    line_kwargs = dict(
        color=color,
        linestyle=linestyle,
        linewidth=linewidth,
        alpha=alpha,
        zorder=zorder,
    )
    if kwargs is not None:
        line_kwargs.update(kwargs)

    if isinstance(ax, np.ndarray):
        axes = ax.flatten()
    else:
        axes = ax if isinstance(ax, Sequence) else [ax]

    for current_ax in axes:
        # Vertical lines
        if vl is not None:
            if isinstance(vl, (int, float)):
                current_ax.axvline(x=vl, **line_kwargs)  # pyright: ignore[reportArgumentType]
            elif isinstance(vl, Iterable):
                for x in vl:
                    current_ax.axvline(x=x, **line_kwargs)  # pyright: ignore[reportArgumentType]

        # Horizontal lines
        if hl is not None:
            if isinstance(hl, (int, float)):
                current_ax.axhline(y=hl, **line_kwargs)  # pyright: ignore[reportArgumentType]
            elif isinstance(hl, Iterable):
                for y in hl:
                    current_ax.axhline(y=y, **line_kwargs)  # pyright: ignore[reportArgumentType]
    return None


def save_plot(
    plot_obj: Union[Figure, Axes, Sequence[Axes], NDArray],
    filename: str,
    dpi: int = 300,
    bbox_inches: str = "tight",
    pad_inches: float = 0.1,
    **kwargs,
) -> None:
    if isinstance(plot_obj, Figure):
        fig = plot_obj
    elif isinstance(plot_obj, Axes):
        fig = plot_obj.figure
    elif isinstance(plot_obj, Sequence) and len(plot_obj) > 0:
        first_ax = plot_obj[0]
        if isinstance(first_ax, Axes):
            fig = first_ax.figure
        else:
            raise ValueError("The first element of the plot_obj is not an Axes object")
    elif isinstance(plot_obj, np.ndarray):  # Assuming plot_obj is a numpy array
        if plot_obj.size == 0:
            raise ValueError("Axes array is empty")
        flat_array = plot_obj.flat
        first_ax = next((ax for ax in flat_array if isinstance(ax, Axes)), None)
        if first_ax is None:
            raise ValueError("No Axes object found in the array")
        fig = first_ax.figure
    else:
        raise TypeError(f"Unsupported type for plot_obj: {type(plot_obj)}")

    fig.savefig(  # pyright: ignore[reportAttributeAccessIssue]
        filename, dpi=dpi, bbox_inches=bbox_inches, pad_inches=pad_inches, **kwargs
    )

    print(f"Plot saved to {filename}")

    return None


def plot_line_diff(
    *args: Union[List[float], NDArray],
    y1_err: Optional[NDArray] = None,
    y2_err: Optional[NDArray] = None,
    axes: Optional[List[Axes]] = None,
    kw_line1: Optional[dict] = None,
    kw_line2: Optional[dict] = None,
    kw_diff: Optional[dict] = None,
    labels: Optional[Sequence[str]] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = r"$T\ [\mu$K]",
    ylabel2: Optional[str] = None,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (4, 4),
    grid: bool = True,
) -> List[Axes]:
    """
    Compare two lines and plot the difference.

    Parameters:
    -----------
    *args : x, y1, y2 or y1, y2
        Either three arrays (x, y1, y2) or two arrays (y1, y2)
    y1_err, y2_err : np.ndarray, optional
        Error on y1 and y2. Defaults to None.
    axes : list of matplotlib Axes, optional
        Existing axes to plot on. Defaults to None to create new axes.
    kw_line1, kw_line2, kw_diff : dict, optional
        Keyword arguments for matplotlib plot function for each line. Defaults are:
        - line1: marker='o', linestyle='', color='blue'
        - line2: marker='', linestyle='--', color='red'
        - diff: marker='s', linestyle='-', color='black'
    labels : tuple of str, optional
        Labels for the three lines. Defaults to None.
    xlabel, ylabel : str, optional
        Axis labels
    ylable2 : str, optional
        Label for the second y-axis. Defaults to None to use the same as ylabel.
    title : str, optional
        Title of the plot
    figsize : tuple of float, optional
        Figure size. Defaults to (4, 4).
    grid : bool
        Whether to show grid

    Returns:
    --------
    list of matplotlib Axes
    """
    # Input validation
    if (length := len(args)) == 1:
        raise ValueError("Two lines are required at least.")
    elif length == 2:
        y1, y2 = args
        # check length consistency
        if len(y1) != len(y2):
            raise ValueError("y1 and y2 must have the same length")
        x = np.arange(len(y1))
    elif length == 3:
        x, y1, y2 = args
        if not (len(x) == len(y1) == len(y2)):
            raise ValueError("x, y1, and y2 must have the same length")
    else:
        raise ValueError("Too many arguments. Expected 2 or 3 arrays.")

    # Create axes if not provided
    if axes is None:
        _, axes = make_figure(
            2,
            1,
            figsize=figsize,
            hspace=0,
            wspace=0,
            sharex=True,
            sharey=False,
            aspect=None,
            gridspec_kw={"height_ratios": [2, 1]},
        )
    if axes is None:
        raise ValueError("Failed to create axes.")

    # Set default plotting parameters
    dkw_line1 = {"marker": ".", "linestyle": "", "color": "royalblue"}
    dkw_line2 = {"marker": "", "linestyle": "--", "color": "tomato"}
    dkw_diff = {"marker": ".", "linestyle": "", "color": "black", "alpha": 0.6}

    # Update with user provided parameters
    if kw_line1 is not None:
        dkw_line1.update(kw_line1)
    if kw_line2 is not None:
        dkw_line2.update(kw_line2)
    if kw_diff is not None:
        dkw_diff.update(kw_diff)

    # Calculate difference
    if isinstance(y1, list):
        y1 = np.array(y1)
    if isinstance(y2, list):
        y2 = np.array(y2)
    diff = y1 - y2

    # Plot lines
    if y1_err is not None:
        line1 = axes[0].errorbar(x, y1, yerr=y1_err, **dkw_line1)  # pyright: ignore[reportArgumentType]
    else:
        line1 = axes[0].plot(x, y1, **dkw_line1)[0]  # pyright: ignore[reportArgumentType]

    if y2_err is not None:
        line2 = axes[0].errorbar(x, y2, yerr=y2_err, **dkw_line2)  # pyright: ignore[reportArgumentType]
    else:
        line2 = axes[0].plot(x, y2, **dkw_line2)[0]  # pyright: ignore[reportArgumentType]

    # Plot difference
    diff_err = None
    if y1_err is not None and y2_err is not None:
        diff_err = np.sqrt(y1_err**2 + y2_err**2)  # pyright: ignore[reportArgumentType]
    elif y1_err is not None:
        diff_err = y1_err
    elif y2_err is not None:
        diff_err = y2_err

    if diff_err is not None:
        linediff = axes[1].errorbar(x, diff, yerr=diff_err, **dkw_diff)  # pyright: ignore[reportArgumentType]
    else:
        linediff = axes[1].plot(x, diff, **dkw_diff)[0]  # pyright: ignore[reportArgumentType]

    # Add labels if provided
    if labels is not None:
        if len(labels) == 2:
            axes[0].legend([line1, line2], labels)
        elif len(labels) == 3:
            axes[0].legend([line1, line2], labels[:-1])
            axes[1].legend([linediff], [labels[-1]])
        else:
            raise ValueError("labels must have length 2 or 3")

    # Set axis labels
    if xlabel:
        axes[1].set_xlabel(xlabel)
    if ylabel:
        axes[0].set_ylabel(ylabel)
        ylabel2 = r"$\Delta$" + ylabel if ylabel2 is None else ylabel2
    if ylabel2:
        axes[1].set_ylabel(ylabel2)

    # Add grid
    if grid:
        axes[0].grid(True, alpha=0.3)
        axes[1].grid(True, alpha=0.3)

    # Add zero reference line to difference plot
    axes[1].axhline(y=0, color="gray", linestyle="--")

    # Set title
    if title:
        axes[0].set_title(title)

    return axes
