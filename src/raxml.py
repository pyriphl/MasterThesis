from Bio.Phylo.Applications import RaxmlCommandline
import os
raxml_cline = RaxmlCommandline(sequences="101.phy",

                               model="GTRGAMMA", algorithm = "x", name="test")
# print(raxml_cline)
# raxml_cline()


os.system('raxmlHPC -f x -m GTRGAMMA -n test2 -p 10000 -s 101.phy')