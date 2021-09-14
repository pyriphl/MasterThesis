import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from tabulate import tabulate
import numpy as np


def show_table(gtr, tree, names):
    print(gtr)
    A = []
    for i in range(len(names)):
        line = [names[i]]
        line.extend(tree[i])
        A.append(line)
    # A = np.vstack((names, tree))
    print(tabulate(A, names))


def get_points_values(table):
    points = []
    values = []
    for i in range(0, len(table[0])):
        for j in range(0, len(table[0])):
            points.append((i, j))
            values.append(table[i][j])
    return points, values


# model type : 'nearest', 'linear', 'cubic'
def plot_graph(table, model):
    grid_x, grid_y = np.mgrid[0:5:100j, 0:5:100j]
    points, values = get_points_values(table)
    # print(points)
    # print((values))
    print(grid_x)
    grid_z = griddata(points, values, (grid_x, grid_y), method=model)
    plt.subplot(221)
    plt.imshow(grid_z.T, extent=(0, 1, 0, 1), origin='lower')
    plt.title(model)
    # plt.gcf().set_size_inches(6, 6)
    plt.show()


def plot_points(t1, t2):
    t1_xs = []
    t1_ys = []
    t1_zs = []
    t2_xs = []
    t2_ys = []
    t2_zs = []
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for i in range(len(t1[0])):
        for j in range(len(t1[0])):
            t1_xs.append(i)
            t1_ys.append(j)
            t1_zs.append(t1[i][j])
            ax.scatter(t1_xs, t1_ys, t1_zs, marker='o')

            t2_xs.append(i)
            t2_ys.append(j)
            t2_zs.append(t2[i][j])
            ax.scatter(t2_xs, t2_ys, t2_zs, marker='^')
    ax.set_xlabel('taxa')
    ax.set_ylabel('taxa')
    ax.set_zlabel('distance')

    plt.show()
