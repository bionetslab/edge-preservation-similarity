"""
Created on Apr 2022

@authors: jkiederle
"""

import time
from utils import *



def compute_similarity(algorithm, G1, G2, time_limit=0, normalize=False):
    '''FUNCTION FOR COMPUTATION OF EDGE PRESERVATION SIMILARITY
        computes the given algorithm and returns the edge preservation similarity

        input:  algorithm: possibilities    'EDGE-PRESERVATION-SIM-APPROX' for approximation
                                            'EDGE-PRESERVATION-SIM-EXACT' for exact measure
                G1:                 first tree as networkx graph object
                G2:                 second tree as networkx graph object
                time_limit:         time limit in seconds, note: only implemented for 'EDGE-PRESERVATION-SIM-EXACT' as it is NP-hard
                normalize:          flag, if true results are normalized by division by max(#edges in trees G1 or G2)

        output: edge preservation similarity value'''

    E=Evaluator()

    print("Beginning computation of edge perservation similarity...")

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

    print("edge preservation similarity: " + str(similarity) + "\n")
    print("computation done in: " + str(tac) + "s\n")

    if timeflag:
        print("Time limit of " + str(time_limit) + "s exceeded.")
    
            
    if normalize:
        similarity = normalize_similarity(similarity, G1, G2)

    return similarity 