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


def terminal_neighbor_dists(self):
    """Return a list of distances between adjacent terminals."""

    def generate_pairs(self):
        pairs = itertools.tee(self)
        pairs[1].next()
        return itertools.izip(pairs[0], pairs[1])

    return [self.distance(*i) for i in generate_pairs(self.find_clades(terminal=True))]

# l_tree = Phylo.read("./data/test/1/l_trees.trees", "newick")
# print(l_tree)
# Phylo.draw_ascii(l_tree)
# g_tree = Phylo.read("./data/test/1/g_trees1.trees", "newick")
# Phylo.draw_ascii(g_tree)
