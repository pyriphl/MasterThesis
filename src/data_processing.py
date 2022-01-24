import numpy
import numpy as np
from scipy.interpolate import griddata, interp1d
from sklearn.linear_model import LinearRegression, LogisticRegression

# model type : 'nearest', 'linear', 'cubic'
from sklearn.model_selection import train_test_split


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
def normalize_minmax_array(min_x, max_x, xs):
    result = np.zeros(xs.shape)
    for i in range(0, len(xs[0])):
        for j in range(0, len(xs[0])):
            result[i][j] = (xs[i][j] - min_x) / (max_x - min_x)
    return result


def normalize_minmax_list(min_x, max_x, xs):
    result = []
    for x in xs:
        result.append((x - min_x) / (max_x - min_x))
    return result


def get_points_values(table):
    points = []
    values = []
    for i in range(0, len(table[0])):
        for j in range(0, len(table[0])):
            points.append((i, j))
            values.append(table[i][j])
    return points, values


# E{[X-E(X)]*[Y-E(Y)]}/E{[X-E(X)]^2}*E{[Y-E(Y)]^2}
# here the expectation is the mean
def correlation(xs, ys):
    x_data = reduce_to_list(xs)
    y_data = reduce_to_list(ys)
    mean_x = numpy.mean(x_data)
    mean_y = numpy.mean(y_data)
    squared_delta_x = [(x - mean_x) ** 2 for x in x_data]
    squared_delta_y = [(y - mean_y) ** 2 for y in y_data]
    numerator = [numpy.mean((x - mean_x) * (y - mean_y)) for x, y in zip(x_data, y_data)]
    # print(numerator)
    denominator = [numpy.mean(x) * numpy.mean(y) for x, y in zip(squared_delta_x, squared_delta_y)]
    return [n / d for n, d in zip(numerator, denominator)]


def linear_regression(x_train, y_train, x_test):
    # Reshape your data either using array.reshape(-1, 1) if your data has a single feature or
    # array.reshape(1, -1) if it contains a single sample.
    x_data = np.array(x_train)
    y_data = np.array(y_train)
    regressor = LinearRegression()
    regressor.fit(x_data.reshape(-1, 1), y_data.reshape(-1, 1))
    y_pred = regressor.predict(x_test.reshape(-1, 1))
    return y_pred.flatten(), regressor.coef_[0][0]


def reduce_to_list(table):
    result = []
    for i in range(0, len(table[0])):
        for j in range(i, len(table[0])):
            result.append(table[i][j])
    return result


def split_model_distance_list(ds, end, step, names):
    n = len(ds[0])
    count = 0
    titles = build_pair_label(names)
    name_pairs = []
    result = []
    x_data = np.arange(0, end, step).tolist()
    for i in range(0, n):
        for j in range(i + 1, n):
            y_data = []
            for d in ds:
                y_data.append(d[j][i])
            # f = interpolation_2d(x_data, y_data, 'cubic')
            result.append(y_data)
            name_pairs.append(titles[count])
            # save fig
            count = count + 1
    return x_data, result, name_pairs


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


def classification(X, y, xtest, model_class, **model_params):
    # Instantiate the model object
    model = model_class(**model_params)
    # Fits the model with the reduced data
    model.fit(X, y)
    y_pred = model.predict(xtest)
    y_pred_proba_1 = model.predict_proba(xtest)[::, 1]
    y_pred_proba_0 = model.predict_proba(xtest)[::, 0]
    return y_pred, y_pred_proba_1, y_pred_proba_0

# def random_forest():
#

# split the training data uniformly
def training_data_split(Xs, ys, num_sample):
    counter = 0
    half = num_sample / 2
    ys_train = []
    Xs_train = []
    ys_test = []
    Xs_test = []
    # search for y = 1
    for i, y in enumerate(ys):
        if y == 1 and counter < half:
            ys_train.append(y)
            Xs_train.append(Xs.iloc[i])
            counter = counter + 1
        if y == 1 and counter >= half:
            ys_test.append(y)
            Xs_test.append(Xs.iloc[i])
    # search for y = 0
    for i, y in enumerate(ys):
        if y == 0 and len(ys_train) < num_sample:
            ys_train.append(y)
            Xs_train.append(Xs.iloc[i])
        if y == 0 and len(ys_train) >= num_sample:
            ys_test.append(y)
            Xs_test.append(Xs.iloc[i])
    return Xs_train, Xs_test, ys_train, ys_test
