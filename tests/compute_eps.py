import time
from utils import *



''' usage via CLI in CLI_eps.py
    usage: python CLI_eps.py [required arguments] [optional arguments]
    more information in CLI_eps.py or in README on git'''



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

    if time_limit > 0:
        print("time limit: " + str(time_limit))
    time_limit_exceeded = False

    tic=time.time()

    timeflag = False
    if algorithm == 'EDGE-PRESERVATION-SIM-EXACT':
        GU=Gurobi_solver(0)
        time_limit_exceeded = GU.compute_duos(G1,G2, time_limit)
        similarity = E.evaluate_sol(G1,G2,GU._sol)

    elif algorithm == 'EDGE-PRESERVATION-SIM-APPROX':
        ALG=Approx_alg()
        ALG.compute_duos(G1,G2)
        similarity = E.evaluate_sol(G1,G2,ALG._sol)

    tac=time.time()
    duration = tac-tic

    if timeflag:
        print("Time limit of " + str(time_limit) + "s exceeded.")
    
            
    if normalize:
        similarity = normalize_similarity(similarity, G1, G2)

    return similarity, duration, time_limit_exceeded