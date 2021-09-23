import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.interpolate import griddata
from tabulate import tabulate
import numpy as np


def show_table(Array, labels, name):
    print('=========================================================================')
    print(name)
    print('=========================================================================')
    A = []
    for i in range(len(labels)):
        line = [labels[i]]
        line.extend(Array[i])
        A.append(line)
    # A = np.vstack((names, tree))
    print(tabulate(A, labels))


def get_points_values(table):
    points = []
    values = []
    for i in range(0, len(table[0])):
        for j in range(0, len(table[0])):
            points.append((i, j))
            values.append(table[i][j])
    return points, values


# model type : 'nearest', 'linear', 'cubic'
def interpolation(table, model):
    grid_x, grid_y = np.mgrid[0:5:300j, 0:5:300j]
    points, values = get_points_values(table)
    grid_z = griddata(points, values, (grid_x, grid_y), method=model)
    return grid_x, grid_y, grid_z


def set_axis_labels(axis, x_label, y_label, z_label, ticker_positions, ticker_labels, figure_name):
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.set_zlabel(z_label)
    axis.set_title(figure_name)

    axis.xaxis.set_major_locator(ticker.FixedLocator(ticker_positions))
    axis.xaxis.set_major_formatter(ticker.FixedFormatter(ticker_labels))

    axis.yaxis.set_major_locator(ticker.FixedLocator(ticker_positions))
    axis.yaxis.set_major_formatter(ticker.FixedFormatter(ticker_labels))


def plot_points(t1, t2, labels):
    t1_xs = []
    t1_ys = []
    t1_zs = []
    t2_xs = []
    t2_ys = []
    t2_zs = []
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(projection='3d')
    for i in range(len(t1[0])):
        for j in range(len(t1[0])):
            t1_xs.append(i)
            t1_ys.append(j)
            t1_zs.append(t1[i][j])
            ax1.scatter(t1_xs, t1_ys, t1_zs, marker='o')

            t2_xs.append(i)
            t2_ys.append(j)
            t2_zs.append(t2[i][j])
            ax1.scatter(t2_xs, t2_ys, t2_zs, marker='^')

    t1_X, t1_Y, t1_Z = interpolation(t1, 'cubic')
    t2_X, t2_Y, t2_Z = interpolation(t2, 'cubic')

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(projection='3d')

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(projection='3d')

    # Plot a basic wireframe.
    ax2.plot_wireframe(t1_X, t1_Y, t1_Z, rstride=10, cstride=10, color='C0')
    ax3.plot_wireframe(t2_X, t2_Y, t2_Z, rstride=10, cstride=10, color='C1')

    positions = []
    for i in range(len(labels)):
        positions.append(i)

    set_axis_labels(ax1, 'seq names', 'seq names', 'distance', positions, labels, 'scattered points')
    set_axis_labels(ax2, 'seq names', 'seq names', 'distance', positions, labels, 'gtr distance')
    set_axis_labels(ax3, 'seq names', 'seq names', 'distance', positions, labels, 'tree distance')

    plt.show()
