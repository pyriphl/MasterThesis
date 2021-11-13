import os
import subprocess

"""variables"""

# -sb f:5e-09       // speciation rate (time units)
# -sd f:4.9e-09     // exctintion rate (time units)
# -gt f:4.9e-10     // global transfer rate
# -lt sl:0.0,1.0,gt // per-family transfer rate
# -lk 0             // transfer recipient does not depend on the evolutionary distance
# -SP f:1           // population size (1 implies not ILS at all)
# -rl f:100         // number of gene family

"""input"""


# -sl f: num_taxa
# -st ln:21.25,0.2 tree height ln(logNormal):location,scale or f:fixed_point
# -su ln:-21.9,0.1  substitution rate ln:location,scale or f:fixed_point
# output: relative path


def tree_generation_simphy(num_taxa: int, height: float, sub_rate: float, output_path: str):
    # command: str = './SimPhy_1.0.2/bin/simphy_lnx64 -sb f:5e-09 -sd f:4.9e-09 -gt f:4.9e-10 -lt sl:0.0,1.0,gt -lk 0 ' \
    #                '-SP f:5 -rl f:10 ' + ' -o ' + output_path + ' -sl f:' + str(num_taxa) + ' -st f:' + str(
    #     height) + ' -su f:' + str(sub_rate)
    # os.system('./SimPhy_1.0.2/bin/simphy_lnx64 -sb f:0.000001 -ld f:0.0000005 -lb f:0.0000005 -lt f:0.0000005 -rs 1 '
    #           '-rs 10 -sp f:10 -sg f:1 -v 2 -od 1 -on 1 -cs 22 ' + ' -o ' + output_path + ' -sl f:'
    # #           + str(num_taxa) + ' -st f:' + str(height) + ' -su f:' + str(sub_rate))
    # os.system('./SimPhy_1.0.2/bin/simphy_lnx64 -sb f:0.000001 -lb f:0.000002 -lt f:0.000005 -rs 10 -rl f:10 -sp f:10 -su f:0.00001 -sg f:1 -sl f:4 -st f:100000 -v 2 -od 1 -op 1 -oc 1 -on 1 '+' -o ' + output_path)
    # os.system('./SimPhy_1.0.2/bin/simphy_lnx64 -rs 2 -rl f:10 -sb ln:-15,1 -st u:200000,20000000 -sl f:5 -so f:1 -sp f:100000 -su f:0.00001 -si f:6 -hh ln:1.2,1 -hl ln:1.4,1 -hg f:200 -v 1 -cs 6656 -om 1 -od 1 -op 1 -oc 1 -on 1 -o' + output_path)
    # subprocess.check_output(command, shell=True)
    # os.system(command)
    # print(command)
    os.system('./SimPhy_1.0.2/bin/simphy_lnx64 -sb f:0.000001 -lb f:0.000002 -lt f:0.000005 -rs 10 -rl f:2 -rg 1 -sp f:10 -su f:0.00001 -sg f:1 -sl f:4 -st f:2 -v 2 -od 1 -op 1 -oc 1 -on 1 '+' -o ' + output_path)


# input_path: path of input tree
# config: path of config file
# output: same as input path
def sequence_generation_indelible(input_path: str, config_path: str):
    # command: str = './SimPhy_1.0.2/scripts/INDELIble_wrapper.pl SimPhy_1.0.2/ILS
    # SimPhy_1.0.2/configuration_files/INDELible_simple.txt 22 1'
    command: str = './SimPhy_1.0.2/scripts/INDELIble_wrapper.pl ' + input_path + ' ' + config_path + ' 22 1 '
    print(command)
    os.system(command)
