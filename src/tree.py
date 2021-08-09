from ete3 import Tree
from Bio import Phylo
# t = Tree('((A,B),D);')
# print(t)
#A = t & "A"
# t.show()

l_tree = Phylo.read("./data/test/1/l_trees.trees", "newick")
print(l_tree)
Phylo.draw_ascii(l_tree)
g_tree = Phylo.read("./data/test/1/g_trees1.trees", "newick")
Phylo.draw_ascii(g_tree)