import numpy as np
from scipy.interpolate import griddata, interp1d


# model type : 'nearest', 'linear', 'cubic'
def interpolation_3d(table, model):
    grid_x, grid_y = np.mgrid[0:5:300j, 0:5:300j]
    points, values = get_points_values(table)
    grid_z = griddata(points, values, (grid_x, grid_y), method=model)
    return grid_x, grid_y, grid_z


def interpolation_2d(xdata, ydata, model):
    x = xdata
    y = ydata
    f = interp1d(x, y, kind=model)
    return f


# x_scaled = (x - min(x))/(max(x) - min(x))
def normalize_minmax(min_x, max_x, xs):
    result = np.zeros(xs.shape)
    for i in range(0, len(xs[0])):
        for j in range(0, len(xs[0])):
            result[i][j] = (xs[i][j] - min_x) / (max_x - min_x)
    return result


def get_points_values(table):
    points = []
    values = []
    for i in range(0, len(table[0])):
        for j in range(0, len(table[0])):
            points.append((i, j))
            values.append(table[i][j])
    return points, values