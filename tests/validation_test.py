"""
Created on Apr 2022

@authors: jkiederle
"""

''' This test should be run using the CLI, please refer to the READ ME on git for detailed explanation of usage
    short version:  type in cmd:
    usage: python validation_test.py [required arguments] [optional arguments]

    required arguments:
        output_path  Path to folder where output should be saved

    optional arguments:
        -h, --help   show this help message and exit'''


import inspect
import os
import sys
import numpy as np
import time
import pandas as pd
import subprocess
import argparse
import itertools
from sklearn.metrics.cluster import normalized_mutual_info_score, mutual_info_score, adjusted_mutual_info_score
from numpy import linalg as LA

UTILS_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
BASE_DIR = os.path.dirname(UTILS_DIR)
sys.path.append(BASE_DIR)
from edge_preservation_similarity.utils import *
from edge_preservation_similarity.compute_eps import *
from parsers import *

def compute_similarities_for_validation(algorithm, input_path, output_path, jar_file_path= '', dataset="rna"):
    '''MAIN FUNCTION FOR USE OF ALGORITHM FOR VALIDATION
        computes the given algorithm and returns similarity values
        note:   exact means the entire matrix is computed. Symmetry is achieved by maximizing over the respective entries
                more information see below (only important for APPROX)

        input:  algorithm: possibilities    'EDGE-PRESERVATION-SIM-APPROX' for approximation
                                            'EDGE-PRESERVATION-SIM-EXACT' for exact measure
                                            'TREE-EDIT-DIST'
                jar_file_path:      path to jar file RTED, necessary if scalability test should be computed with tree edit distance
        saves:  report file (txt file)
                similarity matrix (csv file)
        output: sim_list        list of all flattened similarity values for approximation ratio'''

    print("Beginning computation of similarity with algorithm " + str(algorithm) + "...") 

    #reading in all trees into the graph collection
    graph_coll = []
    if algorithm == 'TREE-EDIT-DIST':
        path = input_path
    
        listdir=os.listdir(path)
        for dire in listdir:
            if not "._" in str(dire):
                graph_coll.append([])
                graphdir=os.listdir(path+'/'+dire)
                for graph_file in graphdir:
                    if not "._" in str(graph_file):
                        graph_coll[-1].append(path + "/" + dire + "/" + graph_file)

    else:
        path = input_path

        graph_coll=import_graph_coll(path)   
        graph_coll=graph_coll_edit(graph_coll)
    
    graph_names_list = import_graph_names(path)
    
    print('path ', path)
    print('dataset: ', dataset)  

    print("graph_names_list: ", graph_names_list)
    

    similarity_matrix = np.zeros((len(graph_names_list),len(graph_names_list)))

    first_time = time.time()

    similarity_matrix, sim_list = compute_similarity_helper(algorithm, graph_coll, similarity_matrix, jar_file_path)
            
    entire_duration = time.time()-first_time

    df_similarity_matrix = pd.DataFrame(data=similarity_matrix, index=graph_names_list, columns=graph_names_list)
 
    sort_matrix(df_similarity_matrix, dataset)

    if algorithm == 'TREE-EDIT-DIST' and dataset == "rna":
        #transform the distance matrix to a similarity matrix for easy comparison
        df_similarity_matrix = dist_similarity_matrix_parser(df_similarity_matrix)
        sort_matrix(df_similarity_matrix, dataset)

    if not algorithm == 'TREE-EDIT-DIST' and not dataset == "rna":
        #transform the similarity matrix to a distance matrix for easy comparison
        df_similarity_matrix = dist_similarity_matrix_parser(df_similarity_matrix)
        sort_matrix(df_similarity_matrix, dataset)

    df_similarity_matrix.to_csv(output_path + '/similarity_' + algorithm + '.csv')

    print("Computation done for " + str(algorithm) + "!")

    return sim_list



def compute_similarity_helper(algorithm, graph_coll, similarity_matrix, jar_file_path=''):
    '''helper funtion for scalability tests
        output:     similarity_matrix   all best pairwise similarities
                    sim_list            list of all flattened similarity values for approximation ratio'''
    
    sim_list = []
    for i in range(len(graph_coll)):
        for j in range(len(graph_coll)):

            if i <= j:
                # above diagonal of matrix    

                sim_values = np.zeros((len(graph_coll[i]),len(graph_coll[j])))
                for ii in range(len(graph_coll[i])):
                    for jj in range(len(graph_coll[j])):
                        
                        G1 = graph_coll[i][ii]
                        G2 = graph_coll[j][jj]

                        if algorithm == 'TREE-EDIT-DIST':
                            sim_values[ii,jj] = 0
                            sim_values[ii,jj] = str(subprocess.check_output(['java', '-jar', 'RTED_v1.2.jar', '-f', G1, G2, '-c', '1', '1', '1', '-s', 'left', '--switch']))[2:-3]
                            similarity = np.min(sim_values)
                        else:
                            sim_values[ii,jj], _, _ = compute_similarity(algorithm, G1, G2, normalize=True)
                            similarity = np.max(sim_values)

                print("similarity: ",similarity)
                sim_list.extend(list(sim_values.flatten()))


                similarity_matrix[i][j] = similarity

                if algorithm == 'TREE-EDIT-DIST':
                    similarity_matrix[j][i] = similarity
       

            else:
                #position below diagonal, similarity of both directions should be computed, the maximum of both similarity values is saved
                #not necessary for 'TREE-EDIT-DIST'

                if algorithm != 'TREE-EDIT-DIST':
                    
                    sim_values = np.zeros((len(graph_coll[i]),len(graph_coll[j])))
                    for ii in range(len(graph_coll[i])):
                        for jj in range(len(graph_coll[j])):
                            
                            G1 = graph_coll[i][ii]
                            G2 = graph_coll[j][jj]

                            sim_values[ii,jj], _, _ = compute_similarity(algorithm, G1, G2, normalize=True)
                    similarity = np.max(sim_values)
                    sim_list.extend(list(sim_values.flatten()))

                    compare_value = similarity_matrix[j][i]

                    if compare_value < similarity:
                        similarity_matrix[i][j] = similarity
                        similarity_matrix[j][i] = similarity
                    else:
                        similarity_matrix[i][j] = similarity_matrix[j][i]

    return similarity_matrix, sim_list

def sort_matrix(df, dataset):
    '''sort matrix by index'''
    if dataset == "rna":
        df.sort_index(key=lambda x: (x.to_series().str[3:].astype(int)), axis = 0, inplace=True)
        df.sort_index(key=lambda x: (x.to_series().str[3:].astype(int)), axis = 1, inplace=True)
    else:
        df.sort_index(axis = 0, inplace=True)
        df.sort_index(axis = 1, inplace=True)


def dist_similarity_matrix_parser(matrix_df):
    '''gets distance matrix and returns similarity matrix and vice versa'''
    matrix_np = matrix_df.to_numpy()
    max_df = matrix_np.max()
    new_np = max_df - matrix_np
    matrix_df.index.name = None
    index_m = matrix_df.index
    new_df = pd.DataFrame(data=new_np, index=index_m, columns=index_m)
    return new_df


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

    summary_df = pd.DataFrame(np.array([[mi_approx, mi_gurobi, mi_tree_edit], [norm_mi_approx, norm_mi_gurobi, norm_mi_tree_edit], [adjusted_mi_approx, adjusted_mi_gurobi, adjusted_mi_tree_edit], [corr_approx, corr_gurobi, corr_tree_edit], [frobenius_approx, frobenius_gurobi, frobenius_tree_edit]]), columns=['EPS_APPROX', 'EPS_EXACT', 'TED'], index=['mi', 'normalized mi', 'adjusted mi', 'correlation coeff', 'frobenius of dist'])

    summary_df.to_csv(out_path+ 'final_validation_results.csv')
    return summary_df


DATASETS = {
    "rna": "RNA",
    "acyclic": "Acyclic",
    "alkane": "Alkane"
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CLI for the validation tests of the edge-preservation-similarity")
    parser.add_argument("output_path", type=str, help="Path to folder where output should be saved")
    parser.add_argument("dataset", type=str, choices=DATASETS.keys(), help="you can choose rna, acyclic or alkane")
    parsed_args = parser.parse_args()

    # 'EDGE-PRESERVATION-SIM-APPROX' for approximation or 'EDGE-PRESERVATION-SIM-EXACT' for exact measure, or 'TREE-EDIT-DIST' for tree edit distance
    algorithms = ['EDGE-PRESERVATION-SIM-APPROX', 'EDGE-PRESERVATION-SIM-EXACT', 'TREE-EDIT-DIST']

    out_path = parsed_args.output_path + "/"
    gml_path = ""
    bracket_path = ""

    helper_gml_path= os.path.normpath("final_results_data/data/" + str(DATASETS[parsed_args.dataset]) + "_trees/" + str(DATASETS[parsed_args.dataset]) + "_trees_gml/")
    gml_path = os.path.abspath(helper_gml_path) + "/"

    #bracket_path = "/home/jana/Documents/BIONETs/Code/test/"
    helper_bracket_path = os.path.normpath("final_results_data/data/" + str(DATASETS[parsed_args.dataset]) + "_trees/" + str(DATASETS[parsed_args.dataset]) + "_trees_bracket/")
    bracket_path = os.path.abspath(helper_bracket_path) + "/"


    if parsed_args.dataset == "rna":
        helper_excel_path = os.path.normpath("final_results_data/data/PKBdatasetGeneOntologyInformationCollection.xlsx")
        excel_path = os.path.abspath(helper_excel_path)
        #get jaccard GO matrix (functional similarity matrix)
        functional_matrix = compute_functional_matrix(excel_path)
    elif parsed_args.dataset == "acyclic":
        print("This is acyclic")
        path_to_functional_similarity_matrix = os.path.abspath(os.path.normpath("final_results_data/data/" + str(DATASETS[parsed_args.dataset]) + "_trees/bp_distances.csv"))
        functional_matrix = pd.read_csv(path_to_functional_similarity_matrix, index_col=0)
    elif parsed_args.dataset == "alkane":
        print("This is alkane")
        path_to_functional_matrix = os.path.abspath(os.path.normpath("final_results_data/data/" + str(DATASETS[parsed_args.dataset]) + "_trees/bp_distances.csv"))
        functional_matrix = pd.read_csv(path_to_functional_matrix, index_col=0)
    else:
        print("You did not provide the name of a dataset.")


    jf_path = 'RTED_v1.2.jar'
    jar_file_path = os.path.abspath(jf_path)

    print("Beginning computation of edge perservation similarity...")
    approx_sim_list = []
    exact_sim_list = []

    for alg in algorithms:
        if alg == 'TREE-EDIT-DIST':
            tree_path = bracket_path
            _ = compute_similarities_for_validation(alg, tree_path, out_path, jar_file_path, parsed_args.dataset)
        elif alg == 'EDGE-PRESERVATION-SIM-APPROX':
            tree_path = gml_path
            approx_sim_list.extend(compute_similarities_for_validation(alg, tree_path, out_path, jar_file_path, parsed_args.dataset))
        else:
            tree_path = gml_path
            exact_sim_list.extend(compute_similarities_for_validation(alg, tree_path, out_path, jar_file_path, parsed_args.dataset))
        
        print("\n")

    approximation_ratios = np.array(exact_sim_list)/np.array(approx_sim_list)
    np.savetxt(out_path + '/approx_ratios.csv', approximation_ratios, delimiter=',', fmt='%f')

    print("Similarity computation done!")
    print("Results saved to: " + str(parsed_args.output_path))


    print("Starting evaluation...")

    df_approx = pd.read_csv(out_path + '/similarity_EDGE-PRESERVATION-SIM-APPROX.csv', index_col=0)
    df_gurobi = pd.read_csv(out_path + '/similarity_EDGE-PRESERVATION-SIM-EXACT.csv', index_col=0)
    df_tree_edit = pd.read_csv(out_path + '/similarity_TREE-EDIT-DIST.csv', index_col=0)         
    
    
    summary_df = evaluate_similarity_results(df_approx, df_gurobi, df_tree_edit, functional_matrix)
    
    print("Evaluation done!\n")
    print(summary_df)
