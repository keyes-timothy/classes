#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:57:57 2019

@author: tkeyes
"""
import time

#tanimoto

self = Tanimoto(drugs_path = './drugs.csv', targets_path = './targets.csv', out_path = './result.csv')

self.load_data()

start = time.time()
self.calc_all_tanimotos() #currently takes 16 seconds. 
end = time.time()
end-start

start = time.time()
self.find_all_shared_targets() #currently takes 33 seconds. 
end = time.time()
end-start


#pvalue
num_iterations = 1000
seed = 214
drug_path = 'drugs.csv'
target_path = 'targets.csv'

protein_a = 'P21918'
protein_b = 'P18089'

protein_a = 'P00734'
protein_b = 'P18089'

protein_a = 'P00734'
protein_b = 'P22888'


my_p = pvalue(drug_path = drug_path, target_path = target_path, num_iterations = 50, seed = 214)
self = my_p

start = time.time()
self.precompute_all_tanimotos() #takes about 18 seconds
end = time.time()
end-start


start = time.time()
T_summaries = self.run_full_analysis(protein_a, protein_b) #takes about 0.01 seconds.
end = time.time()
end-start

args.uniprot_1 = protein_a
args.uniprot_2 = protein_b


#networkgen.py
node_path = 'protein_nodes.csv'
start = time.time()
self = networkgen(drug_path, target_path, node_path) #takes about 16 seconds
end = time.time()
end-start


###################
"""
Select 100 random seeds

For iteration_number in 100, 500, 1000:
        For seed in random seeds:   
            python pvalue.py -r ${seed} -n ${iteration_number} drugs.csv targets.csv P54577 Q7RTX0
        Save results for ${iteration_number} to file
"""

my_random_seeds = range(0, 10000, 100)

#create some data structure

final_list = []
for iteration_number in [100, 500, 1000]:
    iteration_list = []
    for seed in my_random_seeds:
        my_p = pvalue(drug_path = drug_path, target_path = target_path, num_iterations = iteration_number, seed = seed)
        my_p.precompute_all_tanimotos()
        p = my_p.run_full_analysis('P54577', 'Q7RTX0')
        iteration_list.append(p)
    final_list.append(iteration_list)


import statistics
mean_list = []
sd_list = []
for my_list in final_list: 
    my_mean = statistics.mean(my_list)
    my_sd = statistics.stdev(my_list)
    mean_list.append(my_mean)
    sd_list.append(my_sd)
    

