# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
from ete3 import Tree
from Bio.Seq import Seq
def print_hi(seq):
    # Use a breakpoint in the code line below to debug your script.
    # print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.
    my_seq = Seq(seq)
    print(my_seq)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi("AGTACACTGGT")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
