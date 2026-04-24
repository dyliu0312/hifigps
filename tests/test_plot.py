import numpy as np
import pytest
import matplotlib.pyplot as plt
from hifigps.plot import (
    plot_heatmap, plot_heatmaps, set_ticks, get_ticks_labels, make_figure, plot_axlines,
    plot_line, plot_hist, plot_sector, plot_ellipses, plot_arcs, save_plot
)
from hifigps.plot_custom import (
    plot_stack_fit_res, plot_res, plot_profile_2c, plot_profile_2c2r, plot_profile_2r,
    plot_hist_result
)

@pytest.fixture
def test_data_2d():
    shape = (120, 120)
    data = np.zeros(shape)
    x = np.arange(shape[0])
    y = np.arange(shape[1])
    X, Y = np.meshgrid(x, y)
    data += np.exp(-((X - 40)**2 + (Y - 60)**2) / 100) * 100
    data += np.exp(-((X - 80)**2 + (Y - 60)**2) / 100) * 100
    return data

def test_heatmap_plots(test_data_2d):
    # Test plot_heatmap
    ax = plot_heatmap(test_data_2d, title='Test Heatmap')
    assert ax is not None
    plt.close()

    # Test plot_heatmaps
    axes = plot_heatmaps([test_data_2d, test_data_2d], title=['H1', 'H2'])
    assert len(axes) == 2
    plt.close()

def test_line_plots():
    x = [np.linspace(0, 10, 100)]
    y = [np.sin(x[0])]
    ax = plot_line(x, y, title='Test Line')
    assert ax is not None
    plt.close()

def test_hist_plots():
    data = [np.random.normal(size=100)]
    ax = plot_hist(data, bins=10, title='Test Hist')
    assert ax is not None
    plt.close()

def test_custom_plots(test_data_2d):
    # plot_res
    axes = plot_res([test_data_2d, test_data_2d])
    assert len(axes) == 2
    plt.close()

    # plot_stack_fit_res
    axes = plot_stack_fit_res([test_data_2d, test_data_2d, test_data_2d])
    assert len(axes) == 3
    plt.close()

def test_profile_plots():
    x = np.linspace(-3, 3, 120)
    y = np.exp(-x**2)
    
    # plot_profile_2c
    fig, axes = make_figure(1, 2)
    # Provide 4 lines to satisfy default cut=2 (2 for left, 2 for right)
    plot_profile_2c(x=[x, x, x, x], y=[y, y, y, y], axes=axes)
    plt.close()


def test_geometric_plots():
    fig, ax = make_figure()
    plot_sector(ax)
    plot_ellipses(ax, [[0, 0]], [1])
    plot_arcs(ax, [[0, 0]], [1], [0], [90])
    plt.close()

def test_save_plot(test_data_2d, tmp_path):
    ax = plot_heatmap(test_data_2d)
    filename = str(tmp_path / "test_plot.png")
    save_plot(ax, filename)
    assert os.path.exists(filename)
    plt.close()

import os
