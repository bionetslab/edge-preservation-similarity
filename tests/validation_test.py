"""
Created on Apr 2022

@authors: jkiederle
"""

import functools
import os
import networkx as nx
import numpy as np
import time
import pandas as pd
import subprocess
import argparse
import itertools
from sklearn.metrics.cluster import normalized_mutual_info_score, mutual_info_score, adjusted_mutual_info_score
from numpy import linalg as LA

#import edge_preservation_similarity
#from edge_preservation_similarity.utils import *
#from edge_preservation_similarity import * #TODO does not work
from compute_eps import *
from utils import *
from bracket_gml_parser import *
#TODO comment everything


def compute_similarities_for_validation(algorithm, input_path, output_path, jar_file_path= ''):
    '''MAIN FUNCTION FOR USE OF ALGORITHM FOR VALIDATION
        computes the given algorithm and returns similarity values
        note:   exact means the entire matrix is computed. Symmetry is achieved by maximizing over the respective entries
                more information see below (only important for APPROX)

        input:  algorithm: possibilities    'EDGE-PRESERVATION-SIM-APPROX' for approximation
                                            'EDGE-PRESERVATION-SIM-EXACT' for exact measure
                                            'TREE-EDIT-DIST'
                jar_file_path:      path to jar file RTED, necessary if scalability test should be computed with tree edit distance
        saves:  report file (txt file)
                similarity matrix (csv file)'''

    print("Beginning computation of similarity with algorithm" + str(algorithm) + "...") 

    #reading in all trees into the graph collection
    graph_coll = []
    if algorithm == 'TREE-EDIT-DIST':
        path = input_path
    
        dir_list=os.listdir(path)
        for d in dir_list:
            #print(d)
            graph_path = path + "/" + str(d) + "/1.txt"
            graph_coll.append([graph_path])

    else:
        path = input_path

        graph_coll=import_graph_coll(path)   
        graph_coll=graph_coll_edit(graph_coll)
    
    graph_names_list = import_graph_names(path)
    
    print('path ', path)   

    outpath = output_path + '_results_' + str(algorithm)

    os.mkdir(outpath)
    f = open(outpath + "/report.txt","a")

    f.write("Report \n")
    f.write("inpath: " + str(path) + "\n")
    f.write("outpath: " + str(outpath) + "\n")
    f.write("algorithm: " + str(algorithm) + "\n")   
    

    similarity = np.zeros((len(graph_names_list),len(graph_names_list)))

    first_time = time.time()

    similarity_matrix = compute_scalability_helper(algorithm, graph_coll, similarity, jar_file_path)
            
    entire_duration = time.time()-first_time


    f.write("entire duration: " + str(entire_duration) + "\n")

    df_similarity_matrix = pd.DataFrame(data=similarity_matrix, index=graph_names_list, columns=graph_names_list)
 
    sort_matrix(df_similarity_matrix)

    df_similarity_matrix.to_csv(outpath + '/similarity_' + algorithm + '.csv')

    print("Computation done for " + str(algorithm) + "!")



def compute_scalability_helper(algorithm, graph_coll, similarity, jar_file_path=''):
    '''helper funtion for scalability tests'''

    for i in range(len(graph_coll)):
        for j in range(len(graph_coll)):
            G1=graph_coll[i,0]
            G2=graph_coll[j,0]

            if i <= j:
                # above diagonal of matrix              

                tic=time.time()
                duration = 0.0

                if algorithm == 'TREE-EDIT-DIST':
                    #TODO: problem with getting values for TED
                    test = 0
                    test = subprocess.call(['java', '-jar', jar_file_path, '-f', G1, G2, '-c', '1', '1', '1', '-s', 'left', '--switch'])
                    if test != 0:
                        print("something failed in executing tree edit dist jar file at ", i,j)
                        print("error code: ", test)
                    duration = time.time() - tic
                else:
                    sim, _, _ = compute_similarity(algorithm, G1, G2, True)


                similarity[i,j] = sim
                if algorithm != 'EDGE-PRESERVATION-SIM-APPROX':
                    similarity[j,i] = sim

                
                print("Duration ", duration)
            else:
                #position below diagonal, similarity of both directions should be computed, the maximum of both similarity values is saved
                #only necessary for approximation of edge-preservation-similarity as 'EDGE-PRESERVATION-SIM-EXACT' and 'TREE-EDIT-DIST' are exact measures

                if algorithm == 'EDGE-PRESERVATION-SIM-APPROX':
                    compare_value = similarity[j,i]
                    sim, _, _ = compute_similarity(algorithm, G1, G2, True)

                    if compare_value < similarity:
                        similarity[i,j] = similarity
                        similarity[j,i] = similarity
                    else:
                        similarity[i,j] = similarity[j,i]

        if algorithm == 'TREE-EDIT-DIST':
            #transform the distance matrix to a similarity matrix for easy comparison
            similarity = dist_to_similarity_matrix(similarity)

    return similarity

def sort_matrix(df):
    '''sort matrix by index'''
    df.sort_index(key=lambda x: (x.to_series().str[3:].astype(int)), axis = 0, inplace=True)
    df.sort_index(key=lambda x: (x.to_series().str[3:].astype(int)), axis = 1, inplace=True)


def dist_to_similarity_matrix(matrix_df):
    '''gets distance matrix and returns similarity matrix'''
    matrix_np = matrix_df.to_numpy()
    max_df = matrix_np.max()
    sim_np = max_df - matrix_np
    matrix_df.index.name = None
    index_m = matrix_df.index
    sim_df = pd.DataFrame(data=sim_np, index=index_m, columns=index_m)
    sort_matrix(sim_df) 
    return sim_df


#helper functions for computing functional similarity matrix based on GO terms

def compute_jaccard_matrix(df):
    # Iterate through columns and compute jaccard index

    sim_df = pd.DataFrame(columns=df.columns, index=df.columns)
    for col_pair in itertools.combinations(df.columns, 2):
        u1= col_pair[0]
        u2 = col_pair[1]
        sim_df.loc[col_pair] = compute_jaccard(set(df[u1].dropna()), set(df[u2].dropna()))
    
    for i in sim_df.index:
        sim_df[i].loc[i] = 1.0
    
    return sim_df


def compute_jaccard(user1_vals, user2_vals):
    intersection = user1_vals.intersection(user2_vals)
    union = user1_vals.union(user2_vals)
    if float(len(union)) == 0 or len(user1_vals) == 0 or len(user2_vals) == 0:
        return 0.0
    jaccard = len(intersection)/float(len(union))
    return jaccard

#main function to compute functional similarity matrix based on GO terms
def compute_functional_matrix(excel_path):
    excel_df = pd.read_excel (excel_path)
    
    #get all GO data
    GO_data = excel_df.loc[:,"Gene Ontology": "Unnamed: 22"]

    functional_sim_GO = compute_jaccard_matrix(GO_data.T)

    #set row and col indices for easier understanding
    PKB = excel_df['PKBno.']
    functional_sim_GO.set_index(PKB, inplace=True)
    functional_sim_GO.columns = PKB

    #mirror matrix on diagonal for easier comparison
    zeros_functional_sim_GO = functional_sim_GO.fillna(0)
    np_functional_sim_GO = zeros_functional_sim_GO.to_numpy()
    np_functional_sim_GO = np_functional_sim_GO + np_functional_sim_GO.T - np.diag(np.diag(np_functional_sim_GO))

    #make dataframe from functional matrix
    functional_matrix = pd.DataFrame(data=np_functional_sim_GO, index=functional_sim_GO.index, columns=functional_sim_GO.index)

    sort_matrix(functional_matrix)

    return functional_matrix

def evaluate_similarity_results(df_approx, df_gurobi, df_tree_edit, functional_similarity_matrix):
    #cropping is necessary as some PKB values are not used in the reference paper of Wang et al. (2020)
    cropped_functional = functional_similarity_matrix[df_approx.index].loc[df_approx.index]

    np_gurobi = df_gurobi.to_numpy().flatten()
    np_approx = df_approx.to_numpy().flatten()
    np_tree_edit_flat = df_tree_edit.to_numpy().flatten()
    np_crop_func_flat = cropped_functional.to_numpy().flatten()

    dist_approx = df_approx.to_numpy() - cropped_functional.to_numpy()
    dist_gurobi = df_gurobi.to_numpy() - cropped_functional.to_numpy()
    dist_tree_edit = df_tree_edit.to_numpy() - cropped_functional.to_numpy()

    norm_mi_approx = normalized_mutual_info_score(np_approx, np_crop_func_flat)
    adjusted_mi_approx = adjusted_mutual_info_score(np_approx, np_crop_func_flat)
    mi_approx = mutual_info_score(np_approx, np_crop_func_flat)

    norm_mi_gurobi = normalized_mutual_info_score(np_gurobi, np_crop_func_flat)
    adjusted_mi_gurobi = adjusted_mutual_info_score(np_gurobi, np_crop_func_flat)
    mi_gurobi = mutual_info_score(np_gurobi, np_crop_func_flat)

    norm_mi_tree_edit = normalized_mutual_info_score(np_tree_edit_flat, np_crop_func_flat)
    adjusted_mi_tree_edit = adjusted_mutual_info_score(np_tree_edit_flat, np_crop_func_flat)
    mi_tree_edit = mutual_info_score(np_tree_edit_flat, np_crop_func_flat)

    corr_approx = np.corrcoef(np_approx,np_crop_func_flat)[0,1]
    corr_gurobi = np.corrcoef(np_gurobi,np_crop_func_flat)[0,1]
    corr_tree_edit = np.corrcoef(np_tree_edit_flat,np_crop_func_flat)[0,1]

    frobenius_approx = LA.norm(dist_approx)
    frobenius_gurobi = LA.norm(dist_gurobi)
    frobenius_tree_edit = LA.norm(dist_tree_edit)

    summary_df = pd.DataFrame(np.array([[mi_approx, mi_gurobi, mi_tree_edit], [norm_mi_approx, norm_mi_gurobi, norm_mi_tree_edit], [adjusted_mi_approx, adjusted_mi_gurobi, adjusted_mi_tree_edit], [corr_approx, corr_gurobi, corr_tree_edit], [frobenius_approx, frobenius_gurobi, frobenius_tree_edit]]), columns=['approx', 'gurobi', 'tree edit'], index=['mi', 'normalized mi', 'adjusted mi', 'correlation coeff', 'frobenius of dist'])

    summary_df.to_csv(out_path+ 'final_validation_results.csv')
    return summary_df





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CLI for the validation tests of the edge-preservation-similarity")
    parser.add_argument("output_path", type=str, help="Path to folder where output should be saved")
    parser.add_argument("--time_limit", default=0, dest="limit", type=int, help="Set time limit in seconds for exact eps algorithm, only for single computations (default: 0 meaning no time limit)")
    parsed_args = parser.parse_args()

    # 'EDGE-PRESERVATION-SIM-APPROX' for approximation or 'EDGE-PRESERVATION-SIM-EXACT' for exact measure, or 'TREE-EDIT-DIST' for tree edit distance
    algorithms = ['EDGE-PRESERVATION-SIM-APPROX', 'TREE-EDIT-DIST']
    time_limit_gurobi = 600   #10 min but only used for GUROBI

    out_path = parsed_args.output_path + "/"



    helper_gml_path= os.path.normpath("final_results_data/data/RNA_trees/RNA_trees_gml/")
    gml_path = os.path.abspath(helper_gml_path) + "/"

    #bracket_path = "/home/jana/Documents/BIONETs/Code/test/"
    helper_bracket_path = os.path.normpath("final_results_data/data/RNA_trees/RNA_trees_bracket/")
    bracket_path = os.path.abspath(helper_bracket_path) + "/"

    helper_excel_path = os.path.normpath("final_results_data/data/PKBdatasetGeneOntologyInformationCollection.xlsx")
    excel_path = os.path.abspath(helper_excel_path)


    jf_path = 'RTED_v1.2.jar'
    jar_file_path = os.path.abspath(jf_path)

    print("Beginning computation of edge perservation similarity...")

    for alg in algorithms:
        if alg == 'TREE-EDIT-DIST':
            tree_path = bracket_path
        else:
            tree_path = gml_path
        #compute_similarities_for_validation(alg, tree_path, out_path, jar_file_path)
        print("\n")

    print("Similarity computation done!")
    print("Results saved to: " + str(parsed_args.output_path))


    print("Starting evaluation...")

    df_approx = pd.read_csv (out_path + '/similarity_EDGE-PRESERVATION-SIM-APPROX.csv', index_col=0)
    df_gurobi = pd.read_csv (out_path + '/similarity_EDGE-PRESERVATION-SIM-EXACT.csv', index_col=0)
    df_tree_edit = pd.read_csv(out_path + '/similarity_TREE-EDIT-DIST.csv', index_col=0)         
    #get jaccard GO matrix (functional similarity matrix)
    functional_similarity_matrix = compute_functional_matrix(excel_path)

    summary_df = evaluate_similarity_results(df_approx, df_gurobi, df_tree_edit, functional_similarity_matrix)
    
    print("Evaluation done!\n")
    print(summary_df)

                

 
    
    