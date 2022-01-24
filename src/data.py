import numpy
import pandas
from Bio import Phylo
from matplotlib import pyplot as plt

from src.distance_compare import compare_sw_results_linreg
from src.fileIO import create_dir, write_distances, write_distance, write_tags, load_tags, load_distance, load_distances
from src.path import WINDOW_SIZE, STEP_SIZE, PICTURE_PATH, FILE_PATH, CSV_PATH
from src.seq_operations import dist_window_average
from src.tree_operation import pairwise_terminal_dist


# a single dataset with multiple sequences from a single s_tree and one or more g_trees
class dataset:
    def __init__(self, id, sequences, g_trees, s_tree, model):
        self.id = id
        self.model = model
        self.sequences = sequences
        self.seq_length = len(list(sequences.values())[0])
        self.g_trees = g_trees
        self.s_tree = s_tree
        self.taxa_names = []
        self.taxa_numbers = -1
        self.model_distance = numpy.zeros(0)
        self.sliding_window_distance_list = []
        self.s_tree_names = []
        self.tree_distance = numpy.zeros(0)
        self.create_path()

    def __repr__(self):
        representation = 'id: ' + self.id + '\n' + 'model: ' + self.model + '\n' + 'seqs: ' + str(
            self.sequences) + '\n' + 'model dist: \n' + str(self.model_distance)
        return representation

    def calculate_distances_names(self):
        self.taxa_names, self.model_distance, self.sliding_window_distance_list = dist_window_average(self.sequences,
                                                                                                      self.model,
                                                                                                      WINDOW_SIZE)
        self.s_tree_names = [n[0] for n in self.taxa_names]
        self.taxa_numbers = len(self.taxa_names)
        self.tree_distance = pairwise_terminal_dist(self.s_tree_names, self.s_tree)

    def save_dataset(self, path):
        create_dir(path)
        write_distances(self.sliding_window_distance_list, path)
        write_distance(self.model_distance, 'model', path)
        write_distance(self.tree_distance, 'tree', path)
        write_tags(self.taxa_names, path)

    def load_dataset(self, path):
        self.taxa_names = load_tags(path + 'tags.txt')
        self.taxa_numbers = len(self.taxa_names)
        self.model_distance = load_distance(path + 'model_dist.txt')
        self.sliding_window_distance_list = load_distances(path)
        self.s_tree_names = [n[0] for n in self.taxa_names]
        self.tree_distance = load_distance(path + 'tree_dist.txt')

    def create_path(self):
        create_dir(FILE_PATH)
        create_dir(PICTURE_PATH)
        create_dir(CSV_PATH)
        create_dir(PICTURE_PATH+self.id+'/')

    def save_trees(self):
        Phylo.draw(self.s_tree, do_show=False)
        plt.savefig(PICTURE_PATH + self.id + '/s_tree.png')
        counter = 1
        for tree in self.g_trees:
            Phylo.draw(tree, do_show=False)
            plt.savefig(PICTURE_PATH + self.id + '/g_tree'+str(counter)+'.png')
            counter = counter+1

# a group of data that contains multiple datasets, may be used for ML methods
class dataframe:
    def __init__(self, datasets):
        self.datasets = datasets
        self.pd_dataframe = pandas.DataFrame()
        self.ys = []

    def calculate_pd_dataframe(self):
        Xs = numpy.array(())
        is_empty = True
        for data in self.datasets:
            create_dir(PICTURE_PATH + data.id + '/linear_regression/')
            y_preds, errors, slopes, names = compare_sw_results_linreg(data.sliding_window_distance_list,
                                                                       data.seq_length - WINDOW_SIZE, STEP_SIZE,
                                                                       data.taxa_names,
                                                                       PICTURE_PATH + data.id + '/linear_regression/')
            id_names = [data.id + name for name in names]
            slopes_abs = [abs(s) for s in slopes]
            xs = numpy.array((errors, slopes_abs, id_names))
            # print(Xs.shape)
            # print(xs.shape)
            if is_empty:
                Xs = xs
                is_empty = False
            else:
                Xs = numpy.concatenate((Xs, xs), axis=1)
        # print(Xs.shape)
        self.pd_dataframe = pandas.DataFrame(Xs.T, columns=['error', 'slope', 'name'])

    def save_trees(self):
        for data in self.datasets:
            data.save_trees()
    def write_pd_dataframe(self):
        create_dir(FILE_PATH + 'csv')
        self.pd_dataframe.to_csv(CSV_PATH + 'results.csv')

    def read_pd_dataframe(self):
        self.pd_dataframe = pandas.read_csv(CSV_PATH + 'results.csv')

    def read_y(self):
        dataframe = pandas.read_csv(CSV_PATH + 'results_y.csv')
        self.ys = dataframe.y
        self.pd_dataframe = dataframe