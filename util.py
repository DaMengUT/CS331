
import fileinput
import string
import newick
from union_find import *
from itertools import izip,product
from pygraph.classes.graph import graph
from pygraph.readwrite.dot import write
from pygraph.algorithms.traversal import traversal

def read_warnow (fname):
    lines = fileinput.input(fname)
    strings = {}
    for x, y in grouped(lines, 2):
        strings[x.rstrip().lstrip('>')] = y.rstrip()
    return strings

def read_phyML (fname,data_index):
    lines = fileinput.input(fname)
    res = None
    for _ in range(data_index):
        res = {}
        (n,k) = [int(x) for x in string.split(lines[fileinput.lineno()])]
        for i,l in zip(range(n), lines):
            name,seq = string.split(l)
            res[name] = seq
    fileinput.close()
    return res

def read_phyML_tree(fname,data_index):
    tree_string = ""
    for _,l in zip(range(data_index),fileinput.input(fname)):
        tree_string = l

    fileinput.close()
    tree = newick_tree_to_pygraph(tree_string)
    return tree

def newick_tree_to_pygraph(string):
    g = graph()
    l = [0]
    def sub(tree):
        if isinstance(tree, newick.tree.Leaf):
            # if it's just a label
            g.add_node(tree.identifier)
            return tree.identifier
        else:
            # it's a tree and has edges
            l[0] += 1
            my_l = l[0]
            g.add_node(my_l)
            for i,e in enumerate(tree.get_edges()):
                child = sub(e[0])
                g.add_edge((my_l, child))
            return my_l
    t = newick.parse_tree(string)
    sub(t)
    return g

def print_dot(tree,name="graph.png"):
    dot = write(tree)
    f = open("graph.dot", "w")
    f.write(dot)

def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)

def hamming_distance(s1, s2):
    "Return the Hamming distance between equal-length sequences."
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def kruskals(vertex_set, weights):
    pairs = weights.items()
    uf = UnionFind()
    uf.insert_objects(vertex_set)
    total_weight = 0
    gr = graph()
    gr.add_nodes(vertex_set)

    sorted_pairs = sorted(pairs, key=lambda pair: pair[1])
    for edge,weight in sorted_pairs:
        if uf.find(edge[1]) != uf.find(edge[0]):
            uf.union(edge[0],edge[1])
            gr.add_edge(edge)
            total_weight += weight
    return (total_weight,gr)

def calc_fitch_cost(tree,seqs,start_node=None):
    sites = len(seqs.values()[0])
    if not start_node:
        start_node = seqs.keys()[0]
    #now calculate the cost
    totalcost = 0
    nucleotides = ['A', 'C', 'T', 'G']
    def cost (x,y):
        if x == y:
            return 0
        else:
            return 1

    for i in range(sites):
        #loop every node from tip to root
        seen = set([])
        z = {}
        for v in traversal(tree, start_node, 'post'):
            #add the cost of the site to the total
            n = tree.neighbors(v)
            children = filter(lambda x: x in seen, n)

            if len(children) != 0:
                z[v] = {}
                for t in nucleotides:
                    costs = sum([min([cost(t,t_0) + z[c][t_0] for t_0 in nucleotides]) for c in children])
                    z[v][t] = costs
            else:
                z[v] = {}
                the_sequence = seqs[v]
                the_character = the_sequence[i]
                for t in nucleotides:
                    if (the_character == t):
                        z[v][t] = 0
                    else:
                        z[v][t] = 9E8
            seen.add(v)
        totalcost += min(z[start_node].values())
    return totalcost

