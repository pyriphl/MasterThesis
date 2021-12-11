import numpy
from sklearn import metrics

from src.data_processing import split_model_distance_list, linear_regression, normalize_minmax_list

from src.show import show_compare_sw_results_linreg


def compare_sw_results_correlation(ds, end, step, names):
    x_data, data_list, name_pairs = split_model_distance_list(ds, end, step, names)
    n = len(data_list)
    result = []
    tags = []
    for i in range(0, n):
        for j in range(i + 1, n):
            d1 = data_list[i]
            d2 = data_list[j]
            corrcoef = numpy.corrcoef(d1, d2)
            # cov = numpy.cov(d1, d2)
            # print(name_pairs[i])
            # print(d1)
            # print(name_pairs[j])
            # print(d2)
            # print(corrcoef)
            # print(cov)
            tags.append((name_pairs[i], name_pairs[j]))
            result.append(corrcoef[0][1])
    return result, tags


def compare_sw_results_linreg(ds, end, step, names, output):
    x_data, data_list, name_pairs = split_model_distance_list(ds, end, step, names)
    n = len(data_list)
    errors = []
    slopes = []
    y_preds = []
    count = 0
    for y_test in data_list:
        norm_y_test = normalize_minmax_list(numpy.amin(y_test), numpy.amax(y_test), y_test)
        y_pred, slope = linear_regression(x_data, norm_y_test, numpy.array(x_data))
        error = metrics.mean_squared_error(norm_y_test, y_pred)
        # error = metrics.mean_absolute_error(y_test, y_pred)
        # error = np.sqrt(metrics.mean_squared_error(y_test, y_pred))
        y_preds.append(y_pred)
        errors.append(error)
        slopes.append(slope)
        show_compare_sw_results_linreg(x_data, norm_y_test, y_pred, (error, slope), name_pairs[count],
                                       output + name_pairs[count])
        count = count + 1
    return y_preds, errors, slopes, name_pairs
