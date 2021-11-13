import matplotlib.pyplot as plt
from Bio import Phylo
from matplotlib import ticker
from tabulate import tabulate
import numpy as np

from src.data import interpolation_3d, interpolation_2d
from src.path import PICTURE_PATH
from src.seq_operations import WINDOW_SIZE


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


def set_axis_labels(axis, x_label, y_label, z_label, ticker_labels, figure_name):
    ticker_positions = np.arange(len(ticker_labels))
    # for i in range(len(ticker_labels)):
    # ticker_positions.append(i)

    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.set_zlabel(z_label)
    axis.set_title(figure_name)

    axis.xaxis.set_major_locator(ticker.FixedLocator(ticker_positions))
    axis.xaxis.set_major_formatter(ticker.FixedFormatter(ticker_labels))

    axis.yaxis.set_major_locator(ticker.FixedLocator(ticker_positions))
    axis.yaxis.set_major_formatter(ticker.FixedFormatter(ticker_labels))


def plot_surface(data, labels, type: str):
    t1_X, t1_Y, t1_Z = interpolation_3d(data, type)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    # Plot a basic wireframe.
    ax.plot_wireframe(t1_X, t1_Y, t1_Z, rstride=10, cstride=10, color='C0')

    set_axis_labels(ax, 'seq names', 'seq names', 'distance', labels, 'gtr distance')
    plt.show()


def plot_points_scatter(t1, t2, labels):
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

    set_axis_labels(ax1, 'seq names', 'seq names', 'distance', labels, 'scattered points')

    plt.show()


def plot_histogram_3d(t1, labels):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    x_data, y_data = np.meshgrid(np.arange(t1.shape[1]) + 0.25,
                                 np.arange(t1.shape[0]) + 0.25, indexing='ij')
    # Flatten out the arrays so that they may be passed to "ax.bar3d".
    # Basically, ax.bar3d expects three one-dimensional arrays:
    # x_data, y_data, z_data. The following call boils down to picking
    # one entry from each array and plotting a bar to from
    # (x_data[i], y_data[i], 0) to (x_data[i], y_data[i], z_data[i]).
    #
    x_data = x_data.ravel()
    y_data = y_data.ravel()
    z_data = t1.flatten()

    dx = dy = 0.5 * np.ones(len(z_data))

    ax.bar3d(x_data,
             y_data,
             np.zeros(len(z_data)),
             dx, dy, z_data, color='#00ceaa')

    set_axis_labels(ax, 'seq names', 'seq names', 'distance', labels, '3d histogram')
    plt.show()


def build_pair_label(labels):
    # number of taxa entries
    result = []
    n = len(labels)
    # build new labels
    for i in range(0, n):
        for j in range(i + 1, n):
            tag = '(' + labels[i] + ', ' + labels[j] + ')'
            result.append(tag)
    return result


def build_dist_list(t):
    result = []
    # number of taxa entries
    n = len(t[0])
    # build new list
    for i in range(0, n):
        for j in range(i + 1, n):
            result.append(t[i][j])
    return result


def plot_histogram_2d(t1, labels):
    view_labels = build_pair_label(labels)
    pos = np.arange(len(view_labels))  # the label locations
    entry = build_dist_list(t1)
    # print(len(entry))
    # print(len(view_labels))
    plt.barh(pos, entry, align='center', alpha=0.5)
    plt.yticks(pos, view_labels)
    plt.xlabel('distance')
    plt.ylabel('paired seq names')
    plt.title('2d histogram')

    plt.show()


def plot_histogram_2d_onplanes(t1, labels):
    fig = plt.figure()
    entry = t1.reshape(-1)
    pos = np.arange(len(labels))
    ax = fig.add_subplot(projection='3d')
    for i in range(0, len(labels)):
        # Plot the bar graph given by xs and ys on the plane y=k with 80% opacity.
        ax.bar(pos, t1[i], zs=i, zdir='y', alpha=0.9)

    set_axis_labels(ax, 'seq names', 'seq names', 'distance', labels, '2d histogram on differenst planes')

    plt.show()


def show_tree(tree, out_path):
    Phylo.draw(tree, do_show=False)
    plt.savefig(out_path)


def plot_histogram_2d_group(t1, labels):
    n_groups = len(labels)
    # create plot
    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.10
    opacity = 0.8

    for i in range(0, n_groups):
        rect = plt.bar(index + i * bar_width, t1[i], bar_width,
                       alpha=opacity, label=labels[i])

    # rect1 = plt.bar(index+1*bar_width, t1[1], bar_width,
    #                  alpha=opacity, label = labels[1])

    # rect2 = plt.bar(index+2*bar_width, t1[2], bar_width,
    #                  alpha=opacity, label = labels[2])

    plt.xlabel('seq name')
    plt.ylabel('distance')
    plt.title('grouped 2d histogram')
    plt.xticks(index, labels)
    plt.legend()

    plt.tight_layout()
    plt.show()


def plot_histogram_2d_compare(t1, t2, labels, title: str):
    view_labels = build_pair_label(labels)
    pos = np.arange(len(view_labels))  # the label locations
    entry1 = build_dist_list(t1)
    entry2 = build_dist_list(t2)
    bar_width = 0.35

    plt.figure(figsize=(11, 6))
    plt.barh(pos, entry1, bar_width, alpha=0.5, label='model distance')
    plt.barh(pos + bar_width, entry2, bar_width, alpha=0.5, label='tree distance')
    plt.yticks(pos, view_labels)
    plt.xlabel('distance')
    plt.ylabel('paired seq names')
    plt.title(title)
    plt.legend()
    plt.savefig(PICTURE_PATH + title)
    plt.show()


def plot_sliding_window(ts, end, step, names, output):
    x_data = np.arange(0, end, step).tolist()
    n = len(ts[0])  # number of taxa
    titles = build_pair_label(names)
    count = 0

    # labels on the axis
    pos = np.arange(0, end, 50).tolist()  # change according to window size
    labels = []
    for p in pos:
        label = str(p) + '-' + str(p + WINDOW_SIZE)
        labels.append(label)

    # draw
    for i in range(0, n):
        for j in range(i + 1, n):
            fig, ax = plt.subplots(figsize=(8, 8))
            y_data = []
            for t in ts:
                y_data.append(t[j][i])
            f = interpolation_2d(x_data, y_data, 'cubic')
            ax.plot(x_data, f(x_data), '-')
            ax.set_title(titles[count])
            # axis setup
            plt.xticks(pos, labels)
            ax.set_xlabel('nucleotide position')
            ax.set_ylabel('model distance')
            for tick in ax.get_xticklabels():
                tick.set_rotation(55)
            # save fig
            plt.savefig(output + titles[count])  # save fig
            count = count + 1

    plt.show()
