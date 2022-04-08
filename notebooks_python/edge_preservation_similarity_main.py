import os
import numpy as np
import time
from datetime import date, datetime
import pandas as pd
import subprocess
from edge_preservation_similarity import *

#TODO: also be able to compute tree edit distance results (right now it just prints it but it needs to be saved too)
#TODO: scalability function

# this is computing half of the matrix; intention: scalability
def compute_algorithm_timed_fast(algorithm, graph_coll, duos, times, time_limit, time_limit_gurobi, normalize=False, jar_file_path=''):
    '''computes the given algorithm and returns similarity/distance values, runtime and whether a timelimit was reached
        note:   fast means only one half of the matrix is computed and copied to other half (only important for APPROX)

        input:  algorithm: possibilities:   'EDGE-PRESERVATION-SIM-APPROX' for approximation
                                            'EDGE-PRESERVATION-SIM-EXACT' for exact measure, 
                                            or 'TREE-EDIT-DIST' for tree edit distance
                graph_coll:         graph collection with size (TODO)
                duos:               empty matrix in which similarity or distance measure results are returned
                times:              empty matrix in which runtimes of algorithm are returned
                time_limit:         empty matrix in which true/false are returned whether time limit was reached
                time_limit_gurobi:  time limit in seconds, note: only implemented for gurobi as it is NP-hard
                normalize:          flag, if true results are normalized by division by max(#edges in trees G1 or G2) 
        output: duos, times, time_limit - all matrices defined above'''

    for i in range(len(graph_coll)):
        for j in range(len(graph_coll)):
        
            # check whether you are at a position that needs computation (no computation below diagonal of matrix)

            if i <= j:

                G1=graph_coll[i][0]
                G2=graph_coll[j][0]
                value = 0.0

                # compute value with algorithm
                
                tic=time.time()

                timeflag = False
                if algorithm == 'EDGE-PRESERVATION-SIM-EXACT':            
                    timeflag = GU.compute_duos(G1,G2, time_limit_gurobi)
                    print(timeflag)
                    
                elif algorithm == 'EDGE-PRESERVATION-SIM-APPROX':
                    ALG.compute_duos(G1,G2)

                elif algorithm == 'TREE-EDIT-DIST':
                    test = 0
                    test = subprocess.call(['java', '-jar', jar_file_path, '-f', G1, G2, '-c', '1', '1', '1', '-s', 'left', '--switch'])
                    if test != 0:
                        print("something failed in executing tree edit dist jar file at ", i,j)
                        print("error code: ", test)

                tac=time.time()

                # if time limit is reached put 0.0 and put True in time limit matrix
                # else put false in time limit matrix (is already by default false)

                if timeflag:
                    duos[i][j] = 0.0
                    times[i][j] = 0.0
                    time_limit[i][j] = True
                    time_limit[j][i] = True
                else:
                    times[i][j]=tac-tic
                    times[j][i]=tac-tic
                
                print("Duration ", tac-tic)

                
                if algorithm == 'EDGE-PRESERVATION-SIM-EXACT':
                    value = E.evaluate_sol(G1,G2,GU._sol)
                elif algorithm == 'EDGE-PRESERVATION-SIM-APPROX':
                    value = E.evaluate_sol(G1,G2,ALG._sol)

                #normalization only works for edge preservation similarity because of different data format
                if (algorithm != 'TREE-EDIT-DIST') and normalize:
                    value = normalize_value(value, G1, G2)
                
                duos[i][j] = value
                duos[j][i] = value


    return duos, times, time_limit


def compute_algorithm_timed_exact(algorithm, graph_coll, duos, times, time_limit, time_limit_gurobi=0, normalize=False, jar_file_path=''):
    '''computes the given algorithm and returns similarity/distance values, runtime and whether a timelimit was reached
        note:   exact means the entire matrix is computed. Symmetry is achieved by maximizing over the respective entries
                more information see below (only important for APPROX)

        input:  algorithm: possibilities    'EDGE-PRESERVATION-SIM-APPROX' for approximation
                                            'EDGE-PRESERVATION-SIM-EXACT' for exact measure, 
                                            or 'TREE-EDIT-DIST' for tree edit distance
                graph_coll:         graph collection with size (TODO)
                duos:               empty matrix in which similarity or distance measure results are returned
                times:              empty matrix in which runtimes of algorithm are returned
                time_limit:         empty matrix in which true/false are returned whether time limit was reached
                time_limit_gurobi:  time limit in seconds, note: only implemented for gurobi as it is NP-hard
                normalize:          flag, if true results are normalized by division by max(#edges in trees G1 or G2) 
        output: duos, times, time_limit - all matrices defined above'''

    for i in range(len(graph_coll)):
        for j in range(len(graph_coll)):
        
            # compute all positions

            G1=graph_coll[i][0]
            G2=graph_coll[j][0]

            # compute value with algorithm
            # if time limit is reached put 0.0 and put True in time limit matrix
            # else put false in time limit matrix (is already by default false)
            #note: time limit only possible for 'EDGE-PRESERVATION-SIM-EXACT'

            tic=time.time()

            timeflag = False
            if algorithm == 'EDGE-PRESERVATION-SIM-EXACT':              
                timeflag = GU.compute_duos(G1,G2, time_limit_gurobi)
                print(timeflag)
                
            elif algorithm == 'EDGE-PRESERVATION-SIM-APPROX':
                ALG.compute_duos(G1,G2)

            elif algorithm == 'TREE-EDIT-DIST':
                test = 0
                test = subprocess.call(['java', '-jar', jar_file_path, '-f', G1, G2, '-c', '1', '1', '1', '-s', 'left', '--switch'])
                if test != 0:
                    print("something failed in executing tree edit dist jar file at ", i,j)
                    print("error code: ", test)

            tac=time.time()

            if timeflag:
                duos[i][j] = 0.0
                times[i][j] = 0.0
                time_limit[i][j] = True
            else:
                times[i][j] = tac-tic
            
            print("Duration ", tac-tic)

            #final computation of the edge preservation similarity matrix based on the previously constructed duos matrix
            if algorithm == 'EDGE-PRESERVATION-SIM-EXACT':
                value = E.evaluate_sol(G1,G2,GU._sol)
            elif algorithm == 'EDGE-PRESERVATION-SIM-APPROX':
                value = E.evaluate_sol(G1,G2,ALG._sol)
            
            #normalization only works for edge preservation similarity because of different data format
            if algorithm != 'TREE-EDIT-DIST' and normalize:
                value = normalize_value(value, G1, G2)

            #for getting a symmetric matrix
            if j<i:
                if duos[j,i] > value:
                    duos[j,i] = value
                    duos[i,j] = value
                else:
                    duos[i,j] = duos[j,i]
            else:
                duos[i][j] = value

    return duos, times, time_limit


def sort_matrix(df):
    '''sort matrix by row and col indices for easier comparison, inplace
        input:  df: matrix as pandas dataframe'''
    df.sort_index(key=lambda x: (x.to_series().str[3:].astype(int)), axis = 0, inplace=True)
    df.sort_index(key=lambda x: (x.to_series().str[3:].astype(int)), axis = 1, inplace=True)


def normalize_value(value, G1, G2):
    '''normalize similarity value based on the maximum number of edges of the trees G1 and G2
        input:  value: similarity value
                G1, G2: trees of data type networkx graph'''
    edge_count_G1 = G1.number_of_edges()
    edge_count_G2 = G2.number_of_edges()

    return value / max(edge_count_G1,edge_count_G2)

def get_data_scalability(algorithm, path):
    '''gets the tree data which is in two different formats
        input:  algorithm: 'TREE-EDIT-DIST' -> data in bracket format but graph coll just contains file names
                            else: edge preservation similarity -> data in gml format, saved as networkx graphs'''

    graph_coll = []
    if algorithm == 'TREE-EDIT-DIST':
    
        dir_list=os.listdir(path)
        #count = ''
        for d in dir_list:
            #count = str(d)[-1]
            graph_path = path + "/" + str(d) + "/1.txt"
            graph_coll.append([graph_path])

    else:     

        graph_coll=import_graph_coll(path)   
        graph_coll=graph_coll_edit(graph_coll)


    return graph_coll



def compute_algorithm(algorithm, in_path, out_path, jar_file_path='', time_limit=0, normalize=False):

    if algorithm == 'TREE-EDIT-DIST' and jar_file_path == '':
        print("You need to provide the file path to the RTEF_v1.2.jar file if you want to copmute the tree edit distance")
        return

    E=Evaluator()
    ALG=Approx_alg()
    GU=Gurobi_solver(0)

    date_time = str(date.today().strftime("%b_%d_%Y"))+ "-" + str(datetime.now().strftime("%H_%M"))

    #get data   TODO add possibility to use scalability data and actual data at same time
    graph_coll = get_data_scalability(algorithm, in_path)
    
    graph_names_list = import_graph_names(in_path)
    
    print('in_path ', in_path)   

    #for TODO scalability
    #out_path = out_path + str(max_n) + '_results_' + algorithm + "_" + date_time
    #for everything else
    out_path = out_path + '_results_' + algorithm + "_" + date_time

    os.mkdir(out_path)
    f = open(out_path + "/report.txt","a")

    f.write("Report for " + date_time + "\n")
    f.write("in_path: " + str(in_path) + "\n")
    f.write("outpath: " + str(out_path) + "\n")
    #f.write("max_n: " + str(max_n) + "\n") #TODO only do that when you have scalability data
    f.write("algorithm: " + algorithm + "\n")   
    
    if algorithm == 'EDGE-PRESERVATION-SIM-EXACT':
        f.write("time_limit: " + str(time_limit) + "\n")


    duos= np.zeros((len(graph_names_list),len(graph_names_list)))
    times=np.zeros((len(graph_names_list),len(graph_names_list)))
    time_limit_reached =np.zeros((len(graph_names_list),len(graph_names_list)))

    first_time = time.time()

    if algorithm == 'EDGE-PRESERVATION-SIM-APPROX':
        duos, times, time_limit_reached = compute_algorithm_timed_exact(algorithm, graph_coll, duos, times, time_limit_reached, time_limit, jar_file_path)
    else:
        duos, times, time_limit_reached = compute_algorithm_timed_fast(algorithm, graph_coll, duos, times, time_limit_reached, time_limit, jar_file_path)

    entire_duration = time.time()-first_time


    f.write("entire duration: " + str(entire_duration) + "\n")

    df_duos = pd.DataFrame(data=duos, index=graph_names_list, columns=graph_names_list)
    df_duos_times = pd.DataFrame(data=times, index=graph_names_list, columns=graph_names_list)
    df_duos_limit = pd.DataFrame(data=time_limit_reached, index=graph_names_list, columns=graph_names_list)

    #sort the matrices from small to big according to the name of the files (graph_names_list)
    sort_matrix(df_duos)
    sort_matrix(df_duos_times)
    sort_matrix(df_duos_limit)

    df_duos.to_csv(out_path + '/df_values_' + algorithm + '.csv')
    df_duos_times.to_csv(out_path + '/df_times_' + algorithm + '.csv')
    df_duos_limit.to_csv(out_path + '/df_time_limit_reached_' + algorithm + '.csv')


# choose 'EDGE-PRESERVATION-SIM-APPROX' for approximation or 'EDGE-PRESERVATION-SIM-EXACT' for exact measure, or 'TREE-EDIT-DIST' for tree edit distance
algorithm = 'EDGE-PRESERVATION-SIM-APPROX'
time_limit_gurobi = 600   #10 min but only used for GUROBI
all_max_n = [20, 40]

#path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/gml_data"
#gml_path="/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_gml"
#bracket_path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_bracket"

#jar_file_path = '/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/notebooks_python/RTED_v1.2.jar'

result_path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/final_results_data/results/"

E=Evaluator()
ALG=Approx_alg()
GU=Gurobi_solver(0)

date_time = str(date.today().strftime("%b_%d_%Y"))+ "-" + str(datetime.now().strftime("%H_%M"))

path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/gml_data"

compute_algorithm(algorithm, path, result_path, jar_file_path='', time_limit=0, normalize=False)

'''for i in range(len(all_max_n)):

    max_n = all_max_n[i]


    graph_coll = []
    if algorithm == 'TREE-EDIT-DIST':
        bracket_path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_bracket/" + str(max_n)
        path = bracket_path

    else:
        gml_path="/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_gml/" + str(max_n)
        path = gml_path

    compute_algorithm(algorithm, path, result_path, jar_file_path='', time_limit=0, normalize=False)'''

    