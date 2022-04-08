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

max_n = 50

#path = "/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/gml_data"
path="/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/all_other_data/scalability_trees/max_n_"+ str(max_n) +"_trees"
print('path ', path)

outpath = path +'_results_' + str(date.today().strftime("%b_%d_%Y")+ "-" + str(datetime.now().strftime("%H_%M")))

os.mkdir(outpath)
f = open(outpath + "/report.txt","a")

f.write("Report for " + str(date.today().strftime("%b_%d_%Y")+ "-" + str(datetime.now().strftime("%H_%M")) + "\n"))
f.write("inpath: " + str(path) + "\n")
f.write("outpath: " + str(outpath) + "\n")
f.write("max_n: " + str(max_n) + "\n")

     
graph_coll=import_graph_coll(path)    
graph_coll=graph_coll_edit(graph_coll)
graph_names_list = import_graph_names(path)

#print(Sols)
E=Evaluator()
ALG=Approx_alg()
GU=Gurobi_solver(0)

#graph_coll=[graph_coll[i] for i in range(10)]
#graph_names_list=[graph_names_list[i] for i in range(10)]

duos_matrix=np.zeros((2,len(graph_coll),len(graph_coll)))
times_matrix=np.zeros((2,len(graph_coll),len(graph_coll)))

timeout = 600000 #seconds
first_time = time.time()
#dist_alg._name
stop_time_gurobi = False    #if true then gurobi takes too long to compute
stop = False    # if true then stop entire loop

#while time.time() < first_time + timeout:

for i in range(len(graph_coll)):
    for j in range(len(graph_coll)):
        duos= np.zeros((len(graph_coll[i]),len(graph_coll[j])))
        times=np.zeros((len(graph_coll[i]),len(graph_coll[j])))
        for alg_index,dist_alg in enumerate([ALG,GU]):
            if not stop:
                if not (dist_alg._name == 'GUROBI' and stop_time_gurobi):
                    for ii in range(len(graph_coll[i])):
                        for jj in range(len(graph_coll[j])):
                            G1=graph_coll[i][ii]
                            G2=graph_coll[j][jj]
                            tic=time.time()
                            dist_alg.compute_duos(G1,G2)
                            tac=time.time()

                            if (tac - tic) > timeout:
                                if dist_alg._name == 'GUROBI':
                                    stop_time_gurobi = True
                                else:
                                    stop = True

                            duos[ii,jj]=E.evaluate_sol(G1,G2,dist_alg._sol)
                            times[ii,jj]=tac-tic
                            print("Duration ", tac-tic)
                    max_duo=np.max(duos)
                    min_time=times[np.unravel_index(duos.argmax(),duos.shape)]
                    duos_matrix[alg_index,i,j]=max_duo
                    times_matrix[alg_index,i,j]=min_time
                    #print('duo_ratio=',duos_matrix[0,i,j]/duos_matrix[1,i,j],'time ratio=',times_matrix[0,i,j]/times_matrix[1,i,j])
            print("dist_alg name ", dist_alg._name)
print('avg duo ratio=',np.mean(duos_matrix[0]/duos_matrix[1]), 'avg time ratio=',np.mean(times_matrix[0]/times_matrix[1]))

entire_duration = time.time()-first_time

f.write("entire duration: " + str(entire_duration) + "\n")
f.write('avg duo ratio=' + str(np.mean(duos_matrix[0]/duos_matrix[1])))
f.write('avg time ratio=' + str(np.mean(times_matrix[0]/times_matrix[1])))

#np.savetxt("duos_matrix_approx.csv", duos_matrix[0], delimiter=",")
#np.savetxt("duos_matrix_gurobi.csv", duos_matrix[1], delimiter=",")

df_duos_approx = pd.DataFrame(data=duos_matrix[0], index=graph_names_list, columns=graph_names_list)
df_duos_gurobi = pd.DataFrame(data=duos_matrix[1], index=graph_names_list, columns=graph_names_list)

df_duos_times_approx = pd.DataFrame(data=times_matrix[0], index=graph_names_list, columns=graph_names_list)
df_duos_times_gurobi = pd.DataFrame(data=times_matrix[1], index=graph_names_list, columns=graph_names_list)

df_duos_approx.to_csv(outpath + '/df_duos_approx.csv')
df_duos_gurobi.to_csv(outpath + '/df_duos_gurobi.csv')

df_duos_times_approx.to_csv(outpath + '/df_times_approx.csv')
df_duos_times_gurobi.to_csv(outpath + '/df_times_gurobi.csv')
