import import_trees as im
import gurobipy as gu
import os
import networkx as nx
import numpy as np
import time
from datetime import date, datetime
import pandas as pd
import signal
import subprocess
from compute_distances import *


def compute_with_max_n(algorithm, max_n, graph_coll, duos, times, time_limit, time_limit_gurobi):

    for i in range(len(graph_coll)):
        for j in range(len(graph_coll)):
        
            # check whether you are at a position that needs no computation and fill in a 0.0

            if i >= j:
                duos[i][j] = 0.0
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
                if algorithm == 'GUROBI':              
                    timeflag = GU.compute_duos(G1,G2, time_limit_gurobi)
                    print(timeflag)
                    
                elif algorithm == '4-APPROX':
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
                    times[i,j]=tac-tic
                
                print("Duration ", tac-tic)


                if algorithm == 'GUROBI':
                    duos[i,j]=E.evaluate_sol(G1,G2,GU._sol)
                elif algorithm == '4-APPROX':
                    duos[i,j]=E.evaluate_sol(G1,G2,ALG._sol)

    return duos, times, time_limit



# choose '4-APPROX' for approximation or 'GUROBI' for exact measure, or 'TREE-EDIT-DIST' for tree edit distance
algorithm = '4-APPROX'
time_limit_gurobi = 600   #10 min but only used for GUROBI
all_max_n = [20, 40, 60, 80, 100]

#path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/gml_data"
gml_path="/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_gml"
bracket_path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_bracket"

jar_file_path = '/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/notebooks_python/RTED_v1.2.jar'

result_path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/final_results_data/results/"

E=Evaluator()
ALG=Approx_alg()
GU=Gurobi_solver(0)

date_time = str(date.today().strftime("%b_%d_%Y"))+ "-" + str(datetime.now().strftime("%H_%M"))


for i in range(len(all_max_n)):

    max_n = all_max_n[i]

    if max_n == 100:
        print("100")

    graph_coll = []
    if algorithm == 'TREE-EDIT-DIST':
        bracket_path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_bracket/" + str(max_n)
        path = bracket_path
    
        dir_list=os.listdir(path)
        count = ''
        for d in dir_list:
            #print(d)
            graph = ''
            count = str(d)[-1]
            graph_path = path + "/" + str(d) + "/" + str(count) + ".txt"
            graph_coll.append([graph_path])

    else:
        gml_path="/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/tree_blocks_gml/" + str(max_n)
        graph_coll=graph_coll_edit(graph_coll)
        path = gml_path

        graph_coll=import_graph_coll(path)   
        graph_coll=graph_coll_edit(graph_coll)
    
    graph_names_list = import_graph_names(path)
      
    print('path ', path)   

    outpath = result_path + str(max_n) + '_results_' + algorithm + "_" + date_time

    os.mkdir(outpath)
    f = open(outpath + "/report.txt","a")

    f.write("Report for " + date_time + "\n")
    f.write("inpath: " + str(path) + "\n")
    f.write("outpath: " + str(outpath) + "\n")
    f.write("max_n: " + str(max_n) + "\n")
    f.write("algorithm: " + algorithm + "\n")   
    
    if algorithm == 'GUROBI':
        f.write("time_limit: " + str(time_limit_gurobi) + "\n")


    duos= np.zeros((len(graph_names_list),len(graph_names_list)))
    times=np.zeros((len(graph_names_list),len(graph_names_list)))
    time_limit_reached =np.zeros((len(graph_names_list),len(graph_names_list)))

    first_time = time.time()

    duos, times, time_limit_reached = compute_with_max_n(algorithm, max_n, graph_coll, duos, times, time_limit_reached, time_limit_gurobi)
            
    entire_duration = time.time()-first_time


    f.write("entire duration: " + str(entire_duration) + "\n")

    df_duos = pd.DataFrame(data=duos, index=graph_names_list, columns=graph_names_list)

    df_duos_times = pd.DataFrame(data=times, index=graph_names_list, columns=graph_names_list)

    df_duos_limit = pd.DataFrame(data=time_limit_reached, index=graph_names_list, columns=graph_names_list)

    df_duos.to_csv(outpath + '/df_duos_' + algorithm + '.csv')

    df_duos_times.to_csv(outpath + '/df_times_' + algorithm + '.csv')

    df_duos_limit.to_csv(outpath + '/df_time_limit_reached_' + algorithm + '.csv')
    #break