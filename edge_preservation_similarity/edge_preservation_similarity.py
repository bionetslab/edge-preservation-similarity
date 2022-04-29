"""
Created on Apr 2022

@authors: jkiederle
"""

import time
from utils import *
import networkx as nx

import argparse







def compute_similarity(algorithm, G1, G2, time_limit=0, normalize=False):
    '''FUNCTION FOR COMPUTATION OF EDGE PRESERVATION SIMILARITY
        computes the given algorithm and returns the edge preservation similarity

        input:  algorithm: possibilities    'EDGE-PRESERVATION-SIM-APPROX' for approximation
                                            'EDGE-PRESERVATION-SIM-EXACT' for exact measure
                G1:                 first tree as networkx graph object
                G2:                 second tree as networkx graph object
                
                optional:
                time_limit:         time limit in seconds, note: only implemented for 'EDGE-PRESERVATION-SIM-EXACT' as it is NP-hard
                normalize:          flag, if true results are normalized by division by max(#edges in trees G1 or G2)

        output: edge preservation similarity value'''

    E=Evaluator()

    G1 = add_depth(G1)
    G2 = add_depth(G2)

    print("Beginning computation of edge perservation similarity...")
    print("exact or approximated algorithm: " + algorithm)

    if time_limit > 0:
        print("time limit: " + str(time_limit))

    print("normalize similarity: " + str(normalize))

    tic=time.time()

    timeflag = False
    if algorithm == 'EDGE-PRESERVATION-SIM-EXACT':
        GU=Gurobi_solver(0)
        timeflag = GU.compute_duos(G1,G2, time_limit)
        similarity = E.evaluate_sol(G1,G2,GU._sol)

    elif algorithm == 'EDGE-PRESERVATION-SIM-APPROX':
        ALG=Approx_alg()
        ALG.compute_duos(G1,G2)
        similarity = E.evaluate_sol(G1,G2,ALG._sol)

    tac=time.time()
    duration = tac-tic

    print("edge preservation similarity: " + str(similarity) + "\n")
    print("computation done in: " + str(duration) + "s\n")

    if timeflag:
        print("Time limit of " + str(time_limit) + "s exceeded.")
    
            
    if normalize:
        similarity = normalize_similarity(similarity, G1, G2)
    
    print("Computation done!")

    return similarity 



''' This part is for the CLI only, please refer to the READ ME on git for explanation of usage
    short version:  type in cmd:
                    python edge_preservation_similarity.py "string of path of tree1" "string of path of tree 2" 
                    optional:   --time_limit=   your custom timelimit in seconds, data type: int
                                --normalize     flag if you need a normalization'''

ALGORITHMS = {
    "approx": "EDGE-PRESERVATION-SIM-APPROX",
    "exact": "EDGE-PRESERVATION-SIM-EXACT"
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CLI for computing the edge preservation similarity of two trees")
    parser.add_argument("algorithm", type=str, choices=ALGORITHMS.keys(), help="Version of algorithm")
    parser.add_argument("pathG1", type=str, help="Path to tree1")
    parser.add_argument("pathG2", type=str, help="Path to tree2")
    parser.add_argument("--time_limit", default=0, dest="limit", type=int, help="Set time limit in seconds for exact algorithm")
    parser.add_argument("--normalize", action="store_true", help="Normalize similarity by dividing by max nr. of edges in tree1 and tree2")
    parsed_args = parser.parse_args()
    path_G1 = parsed_args.pathG1
    path_G2 = parsed_args.pathG2

    name_of_algorithm = ALGORITHMS[parsed_args.algorithm]

    G1 = nx.read_gml(path_G1, label="id")
    G2 = nx.read_gml(path_G2, label="id")
 
    sim_value= compute_similarity(name_of_algorithm, G1, G2, parsed_args.limit, parsed_args.normalize)
    
    