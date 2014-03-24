#!/usr/bin/env python
from util import *
# data for phylip run: seed = 1393, jumble 12 times
# read in newick tree from phylip
trees = open("outtree").readlines()

for k in range(1,5001):
strings = read_phyML("all_seq",k)
phylip_tree = newick_tree_to_pygraph(trees[0])
phylip_cost = calc_fitch_cost(phylip_tree,strings,1)
true_tree = read_phyML_tree("all_tree",k)
tt_cost = calc_fitch_cost(true_tree,strings,1)
print tt_cost
print phylip_cost

# score it
