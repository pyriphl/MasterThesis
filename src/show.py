import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


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
            t2_zs.append(t1[i][j])
            ax.scatter(t1_xs, t1_ys, t1_zs, marker='^')
    ax.set_xlabel('taxa')
    ax.set_ylabel('taxa')
    ax.set_zlabel('distance')

    plt.show()
