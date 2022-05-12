"""
Created on Apr 2022

@authors: jkiederle
"""

import os
import networkx as nx
import numpy as np
import time
import pandas as pd
import subprocess
import argparse

#import edge_preservation_similarity
#from edge_preservation_similarity.utils import *
#from edge_preservation_similarity import * #TODO does not work
from compute_eps import *
from utils import *
from bracket_gml_parser import *
#TODO comment everything


def compute_scalability(algorithm, max_n_list, input_path, output_path, jar_file_path= ''):
    '''MAIN FUNCTION FOR USE OF ALGORITHM FOR SCALABILITY data
        computes the given algorithm and returns similarity/distance values, runtime and whether a timelimit was reached
        note:   exact means the entire matrix is computed. Symmetry is achieved by maximizing over the respective entries
                more information see below (only important for APPROX)

        input:  algorithm: possibilities    'EDGE-PRESERVATION-SIM-APPROX' for approximation
                                            'EDGE-PRESERVATION-SIM-EXACT' for exact measure
                time_limit:         time limit in seconds, note: only implemented for 'EDGE-PRESERVATION-SIM-EXACT' as it is NP-hard
                jar_file_path:      only necessary if scalability test should be computed with tree edit distance for comparison
        saves:  report file (txt file)
                matrices of duration and whether the runtime was exceeded (csv files)'''

    print("Beginning computation of scalability tests of algorithm" + str(algorithm) + "...")


    for i in range(len(max_n_list)):

        max_n = max_n_list[i]

        graph_coll = []
        if algorithm == 'TREE-EDIT-DIST':
            path = input_path + str(max_n)
        
            dir_list=os.listdir(path)
            for d in dir_list:
                #print(d)
                graph_path = path + "/" + str(d) + "/1.txt"
                graph_coll.append([graph_path])

        else:
            path = input_path + str(max_n)

            graph_coll=import_graph_coll_scal(path)   
            graph_coll=graph_coll_edit(graph_coll)
        
        graph_names_list = import_graph_names(path)
        
        print('path ', path)   

        outpath = output_path + str(max_n) + '_results_' + str(algorithm)

        os.mkdir(outpath)
        f = open(outpath + "/report.txt","a")

        f.write("Report \n")
        f.write("inpath: " + str(path) + "\n")
        f.write("outpath: " + str(outpath) + "\n")
        f.write("max_n: " + str(max_n) + "\n")
        f.write("algorithm: " + str(algorithm) + "\n")   
        
        if algorithm == 'EDGE-PRESERVATION-SIM-EXACT':
            f.write("time_limit: " + str(time_limit_gurobi) + "\n")

        times=np.zeros((len(graph_names_list),len(graph_names_list)))
        time_limit_reached =np.zeros((len(graph_names_list),len(graph_names_list)))

        first_time = time.time()

        duration_matrix, time_limit_exceeded_matrix = compute_scalability_helper(algorithm, graph_coll, times, time_limit_reached, time_limit_gurobi, jar_file_path)
                
        entire_duration = time.time()-first_time


        f.write("entire duration: " + str(entire_duration) + "\n")

        df_duration_matrix = pd.DataFrame(data=duration_matrix, index=graph_names_list, columns=graph_names_list)
        df_time_limit_exceeded_matrix = pd.DataFrame(data=time_limit_exceeded_matrix, index=graph_names_list, columns=graph_names_list)

        df_duration_matrix.to_csv(outpath + '/duration_' + algorithm + '.csv')
        df_time_limit_exceeded_matrix.to_csv(outpath + '/time_limit_reached_' + algorithm + '.csv')

    print("Computation done for " + str(algorithm) + "!")



def compute_scalability_helper(algorithm, graph_coll, times, time_limit, time_limit_gurobi, jar_file_path=''):
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

                time_limit_exceeded = False
                duration = 0.0

                if algorithm == 'TREE-EDIT-DIST':
                    test = 0
                    test = subprocess.call(['java', '-jar', jar_file_path, '-f', G1, G2, '-c', '1', '1', '1', '-s', 'left', '--switch'])
                    if test != 0:
                        print("something failed in executing tree edit dist jar file at ", i,j)
                        print("error code: ", test)
                    duration = time.time() - tic
                else:
                    _, duration, time_limit_exceeded = compute_similarity(algorithm, G1, G2, time_limit_gurobi, False)

                tac=time.time()

                if time_limit_exceeded:
                    times[i][j] = 0.0
                    time_limit[i][j] = True
                else:
                    times[i,j] = duration
                
                print("Duration ", duration)

    return times, time_limit_exceeded



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
            #safe gml file
            nx.write_gml(tree_array[i][j], tree_path)

            #make and safe the bracket version for TREE EDIT DIST usage
            tree_bracket = gml_to_bracket_helper(tree_array[i][j], bracket_string = '', visited = [], node = 0)
            f = open(folder_path + "/" + str(max_n) + "_" + str(j) + "/1.txt","a")
            f.write(tree_bracket)
    

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
    print("Steps: " + str(splitting))

    helper_list = []
    # draw 5 times with replacement to get random list of max_n for each bucket
    for i in range(5):
        random_list = get_pruefer(min_n=splitting[i]+1, max_n=splitting[i+1]+1, len=5)

        #create pruefer codes and trees for each max_n in random_list
        i_helper_list = []
        for m_n in random_list:
            pruefer_code = get_pruefer(max_n=m_n+1, len=(m_n-1))

            pruefer_tree = get_pruefer_tree(pruefer_code)

            #adjust labels and directedness of the trees
            nx.set_node_attributes(pruefer_tree, "A", "lbl")
            pruefer_tree = pruefer_tree.to_directed()

            #print("test")
            i_helper_list.append(pruefer_tree)
        helper_list.append(np.array(i_helper_list, dtype=object))

    random_pruefer_tree_array = np.array(helper_list, dtype=object)
    return random_pruefer_tree_array

def import_graph_coll_scal(path):
    '''reads in gml data and returns a list of lists of graphs''' 
    graph_coll=[]
    listdir=os.listdir(path)
    for dire in listdir:
        print(str(dire))
        graph_coll.append([])
        print(str(path+'/'+dire+'/1.gml'))
        graph_coll[-1].append(nx.readwrite.gml.read_gml(path+'/'+dire+'/1.gml',label='id'))
    return graph_coll




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CLI for the scalability tests of the edge-preservation-similarity")
    parser.add_argument("output_path", type=str, help="Path to folder where output should be saved")
    parser.add_argument("--new_trees", action="store_true", help="flag to compute new scalability trees (default: false, meaning results from paper will be reproduced)")
    parser.add_argument("--time_limit", default=0, dest="limit", type=int, help="Set time limit in seconds for exact eps algorithm (default: 0 meaning no time limit)")
    parsed_args = parser.parse_args()

    # 'EDGE-PRESERVATION-SIM-APPROX' for approximation,'EDGE-PRESERVATION-SIM-EXACT' for exact measure, 'TREE-EDIT-DIST' for tree edit distance
    algorithms = ['EDGE-PRESERVATION-SIM-APPROX', 'EDGE-PRESERVATION-SIM-EXACT', 'TREE-EDIT-DIST']
    time_limit_gurobi = 600   #10 min but only used for GUROBI
    all_max_n = [20, 40, 60, 80, 100]
    out_path = parsed_args.output_path + "/"

    if parsed_args.new_trees:
        print("Computing new trees...")
        max_n = 100
        tree_path = out_path + "/trees"
        os.mkdir(tree_path)
        make_prufer_trees(max_n, tree_path)
        tree_path += "/"
        print("Computing new trees done!")
        print("\n")



    helper_gml_path= os.path.normpath("final_results_data/data/scalability_trees/tree_blocks_gml/")
    gml_path = os.path.abspath(helper_gml_path) + "/"

    #bracket_path = "/home/jana/Documents/BIONETs/Code/test/"
    helper_bracket_path = os.path.normpath("final_results_data/data/scalability_trees/tree_blocks_gml/")
    bracket_path = os.path.abspath(helper_bracket_path) + "/"

    jf_path = 'RTED_v1.2.jar'
    jar_file_path = os.path.abspath(jf_path)


    print("Beginning computation of edge perservation similarity...")

    for alg in algorithms:
        if not parsed_args.new_trees:
            if alg == 'TREE-EDIT-DIST':
                tree_path = bracket_path
            else:
                tree_path = gml_path
        compute_scalability(alg, all_max_n, tree_path, out_path, jar_file_path)
        print("\n")

    print("Scalability computation done!")
    print("Results saved to: " + str(parsed_args.output_path))

                

 
    
    