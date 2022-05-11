#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 18:20:28 2021, edited Apr 2022

@author: nboria, jkiederle
"""

import networkx as nx

#bracket to gml parser

def read_tree(path):
    '''helper function for bracket_to_gml to read in trees'''
    with open(path) as fp:
        tree_string = fp.readline().strip()
    return tree_string

def bracket_to_gml_helper(nodes, edges, tree_string, num_nodes):
    '''helper function for bracket_to_gml to parse the data to gml format'''
    node_id = num_nodes
    nodes.append({'id': node_id, 'label': tree_string[1]})
    num_nodes += 1
    tree_string = tree_string[2:-1]
    start = 0
    num_open = 0
    for i in range(len(tree_string)):
        if tree_string[i] == '{':
            num_open += 1
        elif tree_string[i] == '}':
            num_open -= 1
        if num_open == 0:
            edges.append((node_id, num_nodes))
            num_nodes = bracket_to_gml_helper(nodes, edges, tree_string[start:i+1], num_nodes)
            start = i+1
    return num_nodes

def bracket_to_gml(input_path, output_path):
    '''main function for changing data format of trees from bracket to gml format'''
    tree_string = read_tree(input_path)
    nodes = []
    edges = []
    _ = bracket_to_gml_helper(nodes, edges, tree_string, 0)
    with open(output_path, 'w') as fp:
        fp.write('graph [\n directed 1 \n')
        for node in nodes:
            fp.write(f'  node [\n    id {node["id"]}\n    lbl "{node["label"]}"\n  ]\n')
        for edge in edges:
            fp.write(f'  edge [\n    source {edge[0]}\n    target {edge[1]}\n  ]\n')
        fp.write(']\n')



# gml to bracket parser

def gml_to_bracket_helper(tree, bracket_string, visited, node):
    '''helper function for gml_to_bracket to parse the data format to bracket notation'''
    bracket_string += '{' + str(nx.get_node_attributes(tree, 'lbl')[node])
    visited.append(node)
    for child in tree.successors(node):
        if not child in visited:
            #print("now visiting node: ", child)
            bracket_string = gml_to_bracket_helper(tree, bracket_string, visited, child)
        #else:
            #print("already visited node: ", child)
    bracket_string += '}'
    return bracket_string

def gml_to_bracket(input_path):
    '''main function for changing data format of trees from gml to bracket notation'''
    tree = nx.read_gml(input_path)
    bracket_string = ''
    visited = []
    return gml_to_bracket_helper(tree, bracket_string, visited, 0)




