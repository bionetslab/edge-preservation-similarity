"""
Created on Apr 2022

@authors: jkiederle
"""

import time
from utils import *
import networkx as nx
import pandas as pd
from compute_eps import compute_similarity

import argparse


''' This part is for the CLI only, please refer to the READ ME on git for explanation of usage
    short version:  type in cmd:
    usage: python CLI_eps.py [required arguments] [optional arguments]

    required arguments:
            path for output     string containing the path for the location of the output
            path to trees       string containing one of two options:  
                                - path to .txt file with paths to .gml tree files in each line
                                - paths to all .gml tree files in a row
    optional arguments: 
            --algorithm         possibility to choose version of algorithm, choices: {approx,exact}
                                (default: approx)
            --time_limit=       possibility to set time limit in seconds for exact algorithm (default: 0 meaning no time
                                limit), data type: int
            --normalize         flag to normalize similarity by dividing by max nr. of edges in tree1 and tree2
                                (default: false, meaning no normalization)
            --both_directions   flag to compute similarity between trees in both directions for more    
                                precise output (default: false, meaning just one direction)
            -h, --help          show this help message and exit'''


ALGORITHMS = {
    "approx": "EDGE-PRESERVATION-SIM-APPROX",
    "exact": "EDGE-PRESERVATION-SIM-EXACT"
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CLI for computing the edge preservation similarity of the given trees")
    parser.add_argument("output_path", type=str, help="Path to folder where output should be saved")
    parser.add_argument("graphs", type=str, nargs='+', help="Path to a file with trees or paths to files with trees")
    parser.add_argument("--algorithm", type=str, choices=ALGORITHMS.keys(), help="Version of algorithm (default: approx)")
    parser.add_argument("--time_limit", default=0, dest="limit", type=int, help="Set time limit in seconds for exact algorithm (default: 0 meaning no time limit)")
    parser.add_argument("--normalize", action="store_true", help="Normalize similarity by dividing by max nr. of edges in tree1 and tree2 (default: false meaning no normalization)")
    parser.add_argument("--both_directions", action="store_true", help="Compute similarity between trees in both directions for more precise output (default: false meaning just one direction)")
    parsed_args = parser.parse_args()


    if parsed_args.algorithm != None:
        name_of_algorithm = ALGORITHMS[parsed_args.algorithm]
    else:
        name_of_algorithm = "EDGE-PRESERVATION-SIM-APPROX"

    if len(parsed_args.graphs) == 1:
        #case where we have a file with filepaths to trees as lines
        with open(parsed_args.graphs[0]) as f:
            graph_coll = [line.rstrip() for line in f]
        len_graph_coll = len(graph_coll)
    else:
        #case where we have a list of filepaths to trees
        graph_coll = parsed_args.graphs
        len_graph_coll = len(parsed_args.graphs)

    print("Beginning computation of edge perservation similarity...")
    print("exact or approximated algorithm: " + name_of_algorithm)
    print("normalize similarity: " + str(parsed_args.normalize))

    
    similarity_matrix = np.zeros((len_graph_coll,len_graph_coll))
    duration_matrix = np.zeros((len_graph_coll,len_graph_coll))

    for i in range(len_graph_coll):
        for j in range(len_graph_coll):
            G1=nx.read_gml(graph_coll[i], label="id")
            G2=nx.read_gml(graph_coll[j], label="id")
    
            # check whether you are at a position above the diagonal 
            if i <= j:
                similarity, duration = compute_similarity(name_of_algorithm, G1, G2, parsed_args.limit, parsed_args.normalize)
                similarity_matrix[i,j] = similarity
                duration_matrix[i,j] = duration

                if not parsed_args.both_directions:
                    #if only one direction is computed the matrix is filled with the mirrored values
                    similarity_matrix[j,i] = similarity
                    duration_matrix[j,i] = duration
                


            else:
                #position below diagonal

                if parsed_args.both_directions:
                    #if both directions should be computed the maximum of both similarity values is saved
                    compare_value = similarity_matrix[j,i]
                    similarity, duration = compute_similarity(name_of_algorithm, G1, G2, parsed_args.limit, parsed_args.normalize)

                    if compare_value < similarity:
                        similarity_matrix[i,j] = similarity
                        duration_matrix[i,j] = duration
                        similarity_matrix[j,i] = similarity
                        duration_matrix[j,i] = duration
                    else:
                        similarity_matrix[i,j] = similarity_matrix[j,i]
                        duration_matrix[i,j] = duration_matrix[j,i]

    #parsed_args.output_path
    df_similarity_matrix = pd.DataFrame(data=similarity_matrix)
    df_duration_matrix = pd.DataFrame(data=duration_matrix)
    df_similarity_matrix.to_csv(parsed_args.output_path + '/similarity_' + name_of_algorithm + '.csv')
    df_duration_matrix.to_csv(parsed_args.output_path + '/duration_' + name_of_algorithm + '.csv')

    print("Computation done!")
    print("Results saved to: " + str(parsed_args.output_path))

                

 
    
    