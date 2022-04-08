#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 18:20:28 2021

@author: nboria
"""

import import_trees as im
import os

path="/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/data_wang_jcb2020/data_wang_jcb2020/tree_result_MutiRoot"

outpath="/home/jana/Documents/BIONETs/Code/tree_match_approx_validator/data_test"

list_dir=os.listdir(path)
in_dir=[]
out_dir=[]

for dire in list_dir:
     if not dire.startswith('.'):
        os.mkdir(outpath+'/'+dire)

        list_trees=os.listdir(path+'/'+dire)
        for tree_file in list_trees:

            if not tree_file.startswith('.'):
                im.bracket_to_gml(path+'/'+dire+'/'+tree_file,outpath+'/'+dire+'/'+tree_file)

    
#graph_coll=[]
#for dire in lis_dir