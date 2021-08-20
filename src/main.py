from fileIO import read_from, load_tree

if __name__ == '__main__':
    record = read_from('data/Simphy/HGT/1/data_1.fasta','fasta')
    tree = load_tree('data/Simphy/HGT/1/g_trees1.trees','newick')
    taxa = tree.get_terminals()
    print(taxa)


