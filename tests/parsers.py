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


def gml_to_bracket_helper(tree, bracket_string, visited, node):
    '''helper function for gml_to_bracket to parse the data format to bracket notation'''
    bracket_string += '{' + str(nx.get_node_attributes(tree, 'lbl')[node])
    visited.append(node)
    for child in tree.successors(node):
        if not child in visited:
            bracket_string = gml_to_bracket_helper(tree, bracket_string, visited, child)
    bracket_string += '}'
    return bracket_string


def gml_to_bracket(input_path):
    '''main function for changing data format of trees from gml to bracket notation'''
    tree = nx.read_gml(input_path)
    bracket_string = ''
    visited = []
    return gml_to_bracket_helper(tree, bracket_string, visited, 0)


def read_ct(input_path):
    with open(input_path) as fp:
        tree_strings = [line.strip() for line in fp.readlines()][1:]
    line = 0
    num_nodes = int(tree_strings[line].split()[0])
    num_edges = num_nodes - 1
    tree = nx.Graph()
    line += 1
    for node in range(num_nodes):
        lbl = tree_strings[line].split()[3]
        tree.add_node(node, lbl=lbl)
        line += 1
    for _ in range(num_edges):
        edge = tree_strings[line].split()
        tree.add_edge(int(edge[0]) - 1, int(edge[1]) - 1)
        line += 1
    return tree


def generated_rooted_trees(tree):
    dists = nx.shortest_path_length(tree)
    n = tree.number_of_nodes()
    permutations = dict()
    for source, dists_from_source in dists:
        nodes_with_dists = list(dists_from_source.items())
        nodes_with_dists.sort(key=lambda t: t[1])
        permutations[source] = {nodes_with_dists[new_id][0]: new_id for new_id in range(len(nodes_with_dists))}
    rooted_trees = []
    for root in range(n):
        permutation = permutations[root]
        rooted_tree = nx.DiGraph()
        for old_id, new_id in permutation.items():
            rooted_tree.add_node(new_id, lbl=tree.nodes[old_id]['lbl'])
        for old_u, old_v in tree.edges():
            new_u = permutation[old_u]
            new_v = permutation[old_v]
            if new_u < new_v:
                rooted_tree.add_edge(new_u, new_v)
            else:
                rooted_tree.add_edge(new_v, new_u)
        rooted_trees.append(rooted_tree)
    return rooted_trees


def nx_to_gml(tree, output_path):
    '''main function for changing data format of trees from bracket to gml format'''
    with open(output_path, 'w') as fp:
        fp.write('graph [\n directed 1 \n')
        for id, attributes in tree.nodes(data=True):
            fp.write(f'  node [\n    id {id}\n    lbl "{attributes["lbl"]}"\n  ]\n')
        for u, v in tree.edges():
            fp.write(f'  edge [\n    source {u}\n    target {v}\n  ]\n')
        fp.write(']\n')

