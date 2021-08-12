import os


# -sb f:5e-09       // speciation rate (time units)
# -sd f:4.9e-09     // exctintion rate (time units)
# -gt f:4.9e-10     // global transfer rate
# -lt sl:0.0,1.0,gt // per-family transfer rate
# -lk 0             // transfer recipient does not depend on the evolutionary distance
# -SP f:1           // population size (1 implies not ILS at all)
# -rl f:100         // number of gene family
####input####
# -sl f: num_taxa
# -st ln:21.25,0.2 tree height ln(logNormal):location,scale or f:fixed_point
# -su ln:-21.9,0.1  substitution rate ln:location,scale or f:fixed_point
# output: relative path
def tree_generation_simphy(num_taxa: int, height: float, sub_rate: float, output: str):
    command: str = './SimPhy_1.0.2/bin/simphy_lnx64 -sb f:5e-09 -sd f:4.9e-09 -gt f:4.9e-10 -lt sl:0.0,1.0,gt -lk 0 ' \
                   '-SP f:1 -rl f:100 ' + '-o ../' + output + '-sl f:' + num_taxa + '-st f:' + height + '-su f:' + sub_rate
    # os.system('./SimPhy_1.0.2/bin/simphy_lnx64 -rs 1 -rg 1 -lt f:0.0000005 -o ../SimPhy_test')
    os.system(command)
