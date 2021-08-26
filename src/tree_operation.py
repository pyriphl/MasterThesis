import itertools


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
    return [self.gtr_distance(*i) for i in generate_pairs(self.find_clades(terminal=True))]


def get_edit_distance():
    # TODO
    return
