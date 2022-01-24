import itertools
from collections import defaultdict

import numpy


def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names


def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]


def generate_pairs(self):
    pairs = itertools.tee(self)
    pairs[1].next()
    return itertools.izip(pairs[0], pairs[1])


def terminal_neighbor_dists(self):
    """Return a list of distances between adjacent terminals."""
    return [self.distance(*i) for i in generate_pairs(self.find_clades(terminal=True))]


def pairwise_terminal_dist(names, tree):
    tree_distance = numpy.zeros((len(names), len(names)))
    name_tree = lookup_by_names(tree)
    for i in range(0, len(names)):
        for j in range(0, len(names)):
            tree_distance[i][j] = tree.distance(name_tree.get(names[i]), name_tree.get(names[j]))
    return tree_distance


def pairwise_node_dist(names, tree):
    dists = numpy.zeros((len(names), len(names)))
    # name_tree = lookup_by_names(tree)
    for i in range(0, len(names)):
        for j in range(0, len(names)):
            trace = tree.trace(names[i], names[j])
            # print(names[i] + ', ' + names[j])
            # print(trace)
            dist = len(trace)
            dists[i][j] = dist - 1 if dist > 1 else 0
    return dists

def get_partitions(tree, partition_all):
    terminals = tree.get_terminals()
    non_terminals = tree.get_nonterminals()
    partitions = {}
    # init partitions
    for pair in partition_all:
        partitions.setdefault(pair, 0)

    n = len(terminals) # number of taxa
    for i in range(0,n):
        for j in range(i+1,n):
            for clade in non_terminals:
                if not clade.is_preterminal():
                    continue
                if clade.is_parent_of(terminals[i]) and clade.is_parent_of(terminals[j]):
                    partitions[(terminals[i].name, terminals[j].name)] = clade.branch_length
                    partitions[(terminals[j].name, terminals[i].name)] = clade.branch_length
    for t in terminals:
        partitions[(t.name, 'none')] = t.branch_length
        partitions[('none', t.name)] = t.branch_length
    return partitions