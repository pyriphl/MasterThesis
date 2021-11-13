import numpy as np
from scipy.interpolate import griddata, interp1d

from src.show import get_points_values


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
def normalize_minmax(min, max, xs):
    # TODO
    return