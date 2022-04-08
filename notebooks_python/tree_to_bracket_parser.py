# gml to bracket parser

import networkx as nx
import sys
import subprocess 

#sys.setrecursionlimit(100000) 
jar_path = '/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/notebooks_python/RTED_v1.2.jar'
G1 = '/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_bracket/20/20_1/1.txt'
G2 = '/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_bracket/20/20_3/3.txt'
test = subprocess.call(['java', '-jar', jar_path, '-f', G1, G2, '-c', '1', '1', '1', '-s', 'left', '--switch'])
print(test)

def node_to_bracket(tree, bracket_string, node):
    bracket_string += '{' + str(nx.get_node_attributes(H, 'lbl')[node])
    for child in tree.successors(node):
        print(child)
        bracket_string = node_to_bracket(tree, bracket_string, child)
    bracket_string += '}'
    return bracket_string

def tree_to_bracket(tree):
    bracket_string = ''
    return node_to_bracket(tree, bracket_string, 0)

    # test gml to bracket notation parser

#path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_bracket/20/20_0/0.gml"
#path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/gml_data/PKB4/1.txt"
#H = nx.read_gml(path, label = 'id')
#bracket = tree_to_bracket(H)

#print(bracket)
#f = open("./test.txt","a")
#f = open("/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_bracket/20/20_0/new.txt","a")
#f.write(bracket)