import matplotlib.pyplot as plt
import numpy
import pandas
import sklearn
from Bio import Phylo
from matplotlib import ticker, pyplot as plt
from matplotlib.widgets import Slider
from sklearn import svm
from sklearn.linear_model import LogisticRegression
from tabulate import tabulate
import numpy as np
import seaborn as sns

from src.data_processing import interpolation_3d, interpolation_2d, reduce_to_list, linear_regression, build_pair_label, \
    build_dist_list
from src.path import PICTURE_PATH, WINDOW_SIZE


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


def plot_correlation_sw(xs, ys_list, min_x, max_x):
    fig, ax = plt.subplots()
    plt.title('correlation with sliding window')
    plt.subplots_adjust(left=0.25, bottom=0.25)
    x_data = reduce_to_list(xs)
    y_data_list = [reduce_to_list(ys) for ys in ys_list]
    x_test = np.arange(min_x, max_x, max_x / len(x_data))
    p = plt.scatter(x_data, y_data_list[0], color='gray')

    # x_data = np.arange(0,10,0.5)
    # y_data = 2 * x_data
    # x_test = x_data
    y_pred, = linear_regression(x_data, y_data_list[0], x_test)
    l, = plt.plot(x_test, y_pred, color='red', linewidth=2)
    # l, = plt.plot(x_test, linear_regression(x_data, y_data, x_test), color='red', linewidth=2)

    axslider = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    sposition = Slider(axslider, 'Sliding window pos', 0.0, 499.0, valinit=0.0, valstep=1.0)

    def update(val):
        position = sposition.val
        new_position = numpy.array((x_data, y_data_list[int(position)]))
        p.set_offsets(new_position.transpose())
        y_pred_new = linear_regression(x_data, y_data_list[int(position)])
        l.set_ydata(y_pred_new, x_test)
        fig.canvas.draw_idle()

    sposition.on_changed(update)
    ax.set_xlabel('species tree distance')
    ax.set_ylabel('model distance')
    plt.show()


def show_compare_sw_results_linreg(xs, ys, y_pred, text, title, output):
    fig, ax = plt.subplots()
    # x_data = reduce_to_list(xs)
    x_data = xs
    y_data = ys
    # y_data = reduce_to_list(ys)
    p = plt.scatter(x_data, y_data, color='gray')
    l, = plt.plot(x_data, y_pred, color='red', linewidth=2)
    # ax.set_xlabel('species tree distance')
    # ax.set_ylabel('model distance')
    text = "squared error = " + str(text[0]) + "\nslope = " + str(text[1])
    ax.text(0.5, 0.8, text, fontsize=8, transform=fig.transFigure)
    plt.title(title)
    plt.savefig(output)
    # plt.show()


def plot_boxplot(ts, names):
    # prepare the data
    n = len(names)  # number of taxa
    data = [[] for i in range(n)]
    for i in range(0, n):
        for t in ts:
            data[i] = numpy.concatenate((data[i], t[i]))

    # plot
    fig, ax = plt.subplots()
    ax.boxplot(data, labels=names)
    ax.set_title('Boxplot')
    plt.show()

def plot_confusion_matrix(y_test, y_pred, model):
    cnf_matrix = sklearn.metrics.confusion_matrix(y_test, y_pred)
    class_names = [0, 1]  # name  of classes
    fig1, ax = plt.subplots(figsize=(7, 5))
    tick_marks = np.arange(len(class_names))
    plt.xticks(tick_marks, class_names)
    plt.yticks(tick_marks, class_names)
    # create heatmap
    sns.heatmap(pandas.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu", fmt='g')
    # ax.xaxis.set_label_position("top")
    # plt.tight_layout()
    plt.title('Confusion matrix')
    plt.ylabel('Actual label')
    plt.xlabel('Predicted label')
    plt.savefig(PICTURE_PATH + model+"_heatmap.png")
    # plt.show()
def plot_ROC(y_test, y_pred_proba, model):
    fig2, ax = plt.subplots()
    fpr, tpr, _ = sklearn.metrics.roc_curve(y_test, y_pred_proba)
    auc = sklearn.metrics.roc_auc_score(y_test, y_pred_proba)
    plt.plot(fpr, tpr, label="data 1, auc=" + str(auc))
    plt.legend(loc=4)
    plt.title("Receiver Operating Characteristic(ROC) curve")
    plt.ylabel('True positive rate')
    plt.xlabel('False positive rate')
    plt.savefig(PICTURE_PATH + model +"roc.png")
    # plt.show()

def plot_result_distribution(pd_dataframe):
    fig1, ax = plt.subplots(figsize=(7, 5))
    colums = ['error', 'slope']
    df_y1 = pd_dataframe[pd_dataframe['y']==1]
    df_y0 = pd_dataframe[pd_dataframe['y']==0]
    os_xs = df_y0['error']
    os_ys = df_y0['slope']
    xs_xs = df_y1['error']
    xs_ys = df_y1['slope']
    ax.scatter(os_xs, os_ys, marker='o')
    ax.scatter(xs_xs, xs_ys, marker='x')
    plt.ylabel('slopes')
    plt.xlabel('errors')
    plt.savefig(PICTURE_PATH + "data.png")
def plot_decision_boundary_svm(X, y, kernel_type):
    try:
        X = np.array(X)
        y = np.array(y).flatten()
    except:
        print("Coercing input data to NumPy arrays failed")
    # Reduces to the first two columns of data
    reduced_data = X[:, :2]
    # Instantiate the model object
    model = svm.SVC(kernel=kernel_type)
    # Fits the model with the reduced data
    model.fit(reduced_data, y)


    # Step size of the mesh. Decrease to increase the quality of the VQ.
    h = .00002  # point in the mesh [x_min, m_max]x[y_min, y_max].
    p = 0.001

    # Plot the decision boundary. For that, we will assign a color to each
    x_min, x_max = reduced_data[:, 0].min() - p, reduced_data[:, 0].max() + p
    y_min, y_max = reduced_data[:, 1].min() - p, reduced_data[:, 1].max() + p
    # Meshgrid creation
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

    # Obtain labels for each point in mesh using the model.
    Z = model.predict(np.c_[xx.ravel(), yy.ravel()])

    x_min, x_max = X[:, 0].min() - p, X[:, 0].max() + p
    y_min, y_max = X[:, 1].min() - p, X[:, 1].max() + p
    xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.0001),
                         np.arange(y_min, y_max, 0.0001))

    # Predictions to obtain the classification results
    Z = model.predict(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)

    # Plotting
    plt.contourf(xx, yy, Z, alpha=0.4)
    plt.scatter(X[:, 0], X[:, 1], c=y, alpha=0.8)
    plt.xlabel("Feature-1", fontsize=15)
    plt.ylabel("Feature-2", fontsize=15)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.show()
    return plt
def plot_decision_boundary(X, y, name, model_class, **model_params):
    """
        Function to plot the decision boundaries of a classification model.
        This uses just the first two columns of the data for fitting
        the model as we need to find the predicted value for every point in
        scatter plot.
        Arguments:
                X: Feature data as a NumPy-type array.
                y: Label data as a NumPy-type array.
                model_class: A Scikit-learn ML estimator class
                e.g. GaussianNB (imported from sklearn.naive_bayes) or
                LogisticRegression (imported from sklearn.linear_model)
                **model_params: Model parameters to be passed on to the ML estimator

        Typical code example:
                plt.figure()
                plt.title("KNN decision boundary with neighbros: 5",fontsize=16)
                plot_decision_boundaries(X_train,y_train,KNeighborsClassifier,n_neighbors=5)
                plt.show()
        """
    try:
        X = np.array(X)
        y = np.array(y).flatten()
    except:
        print("Coercing input data to NumPy arrays failed")
    # Reduces to the first two columns of data
    reduced_data = X[:, :2]
    # Instantiate the model object
    model = model_class(**model_params)
    # Fits the model with the reduced data
    model.fit(reduced_data, y)


    # Step size of the mesh. Decrease to increase the quality of the VQ.
    h = .0001  # point in the mesh [x_min, m_max]x[y_min, y_max].
    p = 0.001

    # Plot the decision boundary. For that, we will assign a color to each
    x_min, x_max = reduced_data[:, 0].min() - p, reduced_data[:, 0].max() + p
    y_min, y_max = reduced_data[:, 1].min() - p, reduced_data[:, 1].max() + p
    # Meshgrid creation
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

    # Obtain labels for each point in mesh using the model.
    Z = model.predict(np.c_[xx.ravel(), yy.ravel()])

    x_min, x_max = X[:, 0].min() - p, X[:, 0].max() + p
    y_min, y_max = X[:, 1].min() - p, X[:, 1].max() + p
    xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.00001),
                         np.arange(y_min, y_max, 0.00001))

    # Predictions to obtain the classification results
    Z = model.predict(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)

    # Plotting
    fig1, ax = plt.subplots(figsize=(10, 5))
    plt.contourf(xx, yy, Z, alpha=0.4)
    plt.scatter(X[:, 0], X[:, 1], c=y, alpha=0.8)
    plt.xlabel("error", fontsize=15)
    plt.ylabel("slope", fontsize=15)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(name)
    plt.savefig(PICTURE_PATH + name+".png")
    # plt.show()
    return plt


# def plot_decision_boundary( X, y, X_train, y_train):
#     # X_array = X.values
#     X_array = np.array(X_train)
#     X = X.to_numpy()
#     # define bounds of the domain
#     min1, max1 = X_array[:, 0].min() - 0.001, X_array[:, 0].max() + 0.001
#     min2, max2 = X_array[:, 1].min() - 0.001, X_array[:, 1].max() + 0.001
#     # define the x and y scale
#     x1grid = np.arange(min1, max1, 0.0001)
#     x2grid = np.arange(min2, max2, 0.0001)
#     # create all of the lines and rows of the grid
#     xx, yy = np.meshgrid(x1grid, x2grid)
#     # flatten each grid to a vector
#     r1, r2 = xx.flatten(), yy.flatten()
#     r1, r2 = r1.reshape((len(r1), 1)), r2.reshape((len(r2), 1))
#     # horizontal stack vectors to create x1,x2 input for the model
#     grid = np.hstack((r1, r2))
#     # define the model
#     model = LogisticRegression()
#     # fit the model
#     model.fit(X_train, y_train)
#     # make predictions for the grid
#     yhat = model.predict(grid)
#     # yhat = model.predict_proba(grid)
#     # keep just the probabilities for class 0
#     # yhat = yhat[:, 1]
#     # reshape the predictions back into a grid
#     zz = yhat.reshape(xx.shape)
#     # plot the grid of x, y and z values as a surface
#     plt.contourf(xx, yy, zz, cmap='Paired')
#     # c = plt.contourf(xx, yy, zz, cmap='RdBu')
#     # add a legend, called a color bar
#     # plt.colorbar(c)
#     # create scatter plot for samples from each class
#     for class_value in range(2):
#         # get row indexes for samples with this class
#         row_ix = np.where(y == class_value)
#         # create scatter of these samples
#         plt.scatter(X[row_ix, 0], X[row_ix, 1], cmap='Paired')
#     plt.show()