"""
Created on Apr 2022

@authors: jkiederle
"""

import os
import networkx as nx
import numpy as np
import time
from datetime import date, datetime
import pandas as pd
import subprocess
from edge_preservation_similarity import *


def compute_scalability(algorithm, max_n_list, in_path, out_path, jar_file_path= ''):
    '''MAIN FUNCTION FOR USE OF ALGORITHM FOR SCALABILITY data
        computes the given algorithm and returns similarity/distance values, runtime and whether a timelimit was reached
        note:   exact means the entire matrix is computed. Symmetry is achieved by maximizing over the respective entries
                more information see below (only important for APPROX)

        input:  algorithm: possibilities    'EDGE-PRESERVATION-SIM-APPROX' for approximation
                                            'EDGE-PRESERVATION-SIM-EXACT' for exact measure
                time_limit:         time limit in seconds, note: only implemented for 'EDGE-PRESERVATION-SIM-EXACT' as it is NP-hard
                normalize:          flag, if true results are normalized by division by max(#edges in trees G1 or G2)
                fast:               flag, if true only upper half of matrix is computed and then mirrored
                                        -> is less exact for 'EDGE-PRESERVATION-SIM-APPROX'
                jar_file_path:      only necessary if scalability test should be computed with tree edit distance for comparison
        saves:  report file (txt file)
                matrices of computed values, runtime and whether the runtime expired (csv files)'''



    E=Evaluator()
    ALG=Approx_alg()
    GU=Gurobi_solver(0)

    date_time = str(date.today().strftime("%b_%d_%Y"))+ "-" + str(datetime.now().strftime("%H_%M"))


    for i in range(len(max_n_list)):

        max_n = max_n_list[i]

        if max_n == 100:
            print("100")

        graph_coll = []
        if algorithm == 'TREE-EDIT-DIST':
            path = in_path + str(max_n)
        
            dir_list=os.listdir(path)
            count = ''
            for d in dir_list:
                #print(d)
                graph = ''
                count = str(d)[-1]
                graph_path = path + "/" + str(d) + "/1.txt"
                graph_coll.append([graph_path])

        else:
            path = in_path + str(max_n)

            graph_coll=import_graph_coll(path)   
            graph_coll=graph_coll_edit(graph_coll)
        
        graph_names_list = import_graph_names(path)
        
        print('path ', path)   

        outpath = out_path + str(max_n) + '_results_' + algorithm + "_" + date_time

        os.mkdir(outpath)
        f = open(outpath + "/report.txt","a")

        f.write("Report for " + date_time + "\n")
        f.write("inpath: " + str(path) + "\n")
        f.write("outpath: " + str(outpath) + "\n")
        f.write("max_n: " + str(max_n) + "\n")
        f.write("algorithm: " + algorithm + "\n")   
        
        if algorithm == 'EDGE-PRESERVATION-SIM-EXACT':
            f.write("time_limit: " + str(time_limit_gurobi) + "\n")


        duos= np.zeros((len(graph_names_list),len(graph_names_list)))
        times=np.zeros((len(graph_names_list),len(graph_names_list)))
        time_limit_reached =np.zeros((len(graph_names_list),len(graph_names_list)))

        first_time = time.time()

        times, time_limit_reached = compute_scalability_helper(algorithm, graph_coll, times, time_limit_reached, time_limit_gurobi, E, ALG, GU, jar_file_path)
                
        entire_duration = time.time()-first_time


        f.write("entire duration: " + str(entire_duration) + "\n")

        df_duos_times = pd.DataFrame(data=times, index=graph_names_list, columns=graph_names_list)

        df_duos_limit = pd.DataFrame(data=time_limit_reached, index=graph_names_list, columns=graph_names_list)

        df_duos_times.to_csv(outpath + '/df_times_' + algorithm + '.csv')

        df_duos_limit.to_csv(outpath + '/df_time_limit_reached_' + algorithm + '.csv')




def compute_scalability_helper(algorithm, graph_coll, times, time_limit, time_limit_gurobi, E, ALG, GU, jar_file_path=''):
    '''helper funtion for scalability tests'''

    for i in range(len(graph_coll)):
        for j in range(len(graph_coll)):
        
            # check whether you are at a position that needs no computation and fill in a 0.0

            if i >= j:
                times[i][j] = 0.0
                time_limit[i][j] = False

            else:

                G1=graph_coll[i][0]
                G2=graph_coll[j][0]

                # compute value with algorithm
                # if time limit is reached put 0.0 and put True in time limit matrix
                # else put false in time limit matrix (is already by default false)

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
                    times[i][j] = 0.0
                    time_limit[i][j] = True
                else:
                    times[i,j]=tac-tic
                
                print("Duration ", tac-tic)

    return times, time_limit



'''make random trees for scalability tests'''

def make_prufer_trees(max_n, out_path):
    '''MAIN FUNCTION FOR CREATING SCALABILITY TREES AND SAVING AT RIGHT POSITION
        safe the random trees in 10 blocks up to size max_n by using prufer codes
        -> save in gml format '''
    tree_array = get_random_pruefer_trees(max_n)

    for i in range(tree_array[0].size):
        #determine the name of the folder
        max_n = (i+1)*20
        folder_path = out_path + "/" + str(max_n)
        os.mkdir(folder_path)
        #make folder
        for j in range(tree_array[0].size):
            #determine the name of the file
            tree_folder_path = folder_path + "/" + str(max_n) + "_" + str(j)
            os.mkdir(tree_folder_path)
            tree_path = tree_folder_path + '/1.gml'
            #safe file
            nx.write_gml(tree_array[i][j], tree_path)
    

def get_pruefer(min_n=0, max_n=10, len=8):
    '''make pruefercode with max n'''
    return list(np.random.randint(low=min_n, high=(max_n), size=(len)))

def get_pruefer_tree(pruefer_code):
    '''make tree from pruefercode'''
    return nx.from_prufer_sequence(pruefer_code)

def get_random_pruefer_trees(max_n):
    '''make pruefercodes and random trees of different sizes from them
        input:      max_n - maximum height of any generated tree
        output:     array of size (10,10) - 10 buckets of equal size from 5 to max_n (bucket list)
                     - with each: 10 random pruefer trees with height drawn with replacement from bucket list'''

    # split 1 to max_n in 10 buckets
    splitting = np.linspace(1, max_n, num=6).astype(int)
    print(splitting)

    helper_list = []
    # draw 5 times with replacement to get random list of max_n for each bucket
    for i in range(5):
        random_list = get_pruefer(min_n=splitting[i]+1, max_n=splitting[i+1]+1, len=5)
        print(random_list)

        #create pruefer codes and trees for each max_n in random_list
        i_helper_list = []
        for m_n in random_list:
            pruefer_code = get_pruefer(max_n=m_n+1, len=(m_n-1))
            print("max_n", m_n)
            print(len(pruefer_code))
            pruefer_tree = get_pruefer_tree(pruefer_code)
            print(pruefer_tree.size())

            #adjust labels and directedness of the trees
            nx.set_node_attributes(pruefer_tree, "A", "lbl")
            pruefer_tree = pruefer_tree.to_directed()

            #print("test")
            i_helper_list.append(pruefer_tree)
        helper_list.append(np.array(i_helper_list, dtype=object))

    random_pruefer_tree_array = np.array(helper_list, dtype=object)
    return random_pruefer_tree_array


max_n = 100
#tree_array = get_random_pruefer_trees(max_n)

path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees"

make_prufer_trees(max_n, path)


# choose 'EDGE-PRESERVATION-SIM-APPROX' for approximation or 'EDGE-PRESERVATION-SIM-EXACT' for exact measure, or 'TREE-EDIT-DIST' for tree edit distance
algorithm = 'TREE-EDIT-DIST'
time_limit_gurobi = 600   #10 min but only used for GUROBI
all_max_n = [20, 40]

#path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/gml_data"
gml_path="/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_gml/"
bracket_path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_bracket/"

jar_file_path = '/home/jana/Documents/BIONETs/Code/edge-preservation-similarity/notebooks_python/RTED_v1.2.jar'

result_path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/final_results_data/results/"

#compute_scalability(algorithm, all_max_n, bracket_path, result_path, jar_file_path)
