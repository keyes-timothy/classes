#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 09:59:44 2019

@author: tkeyes
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 11:34:52 2019

@author: tkeyes
"""

import sys
import argparse
import numpy as np
import os
import pandas as pd
import random
from tanimoto import Tanimoto 
from pvalue import pvalue


class networkgen: 
    
    def __init__(self, drug_path, targets_path, nodes_path): 
        """
        Reads in relevant files and initializes the datastructures for the networkgen object. 
        
        Inputs: 
            drug_path = path to the drugs.csv file
            targets_path = path to the targets.csv file
            nodes_path = path to the protein_nodes.csv file
        Returns:
            Nothing
        Side-effects: 
            loads self with a pvalue object and a pandas dataFrame storing the information in nodes.csv
        """
        #create a pvalue object that will do most of the heavy lifting
        self.my_pvalue = pvalue(num_iterations = 500, seed = 214, drug_path = drug_path, target_path = targets_path)
        self.my_pvalue.precompute_all_tanimotos()
        
        #save paths
        self.nodes_df = pd.read_csv(nodes_path)
        self.pair_list = self.find_unique_pairs(self.nodes_df.uniprot_accession.values)
    
    def find_unique_pairs(self, protein_array): 
        """
        A function that calculates all unique pairs of proteins in a list of protein names (protein_array). 
        
        
        Inputs: 
            protein_array = array of strings in which each string represents a protein that will serve as one of the nodes of our network. 
            
        Returns: 
            pair_list = a list of tuples in which each tuple represents a unique pairing of two proteins from our input protein_nodes.csv file. 
        """
        pair_set = set()
        for protein_1 in protein_array: 
            for protein_2 in protein_array: 
                if protein_1 != protein_2: 
                    new_set = frozenset([protein_1, protein_2])
                    pair_set.add(new_set)
                    
        pair_list = [tuple(pair) for pair in pair_set]
        return pair_list
        
    def compute_pvalues(self): 
        """
        A funtion that computes the p-values for the ligand set associations between all unique pairs of proteins in protein_nodes.csv
        
        Inputs: 
            None
            
        Returns: 
            pvalues = a list containing the p-values associated with each pair of proteins in the input self.pair_list. 
        """
        pvalues = []
        for pair in self.pair_list: 
            protein_a = pair[0]
            protein_b = pair[1]
            p_value = self.my_pvalue.run_full_analysis(protein_a, protein_b)
            pvalues.append(p_value)
        return pvalues
    
    def write_output(self, pvalues): 
        """
        A function that writes the requested output for networkgen.py. Specifically, a .csv file is written that contains 
        the edgelist for all protein pairs with an associated p-value less than 0.05.
        
        Inputs: 
            pvalues = a list containing the p-values associated with each pair of proteins in the input self.pair_list. 
            
        Returns: 
            Nothing
        
        Side-effects: 
            Writes a .csv file named "network_edgelist.txt" that contains the edgelist as specified in the project page. 
            
        """
        #separate tuples from one another
        protein_1_list, protein_2_list = zip(*self.pair_list)
        
        #create pandas dataFrame
        output_frame = pd.DataFrame(list(zip(protein_1_list, protein_2_list, pvalues)))
        
        #select rows with pvalue <0.05
        final_frame = output_frame.loc[output_frame[2] <= 0.05].iloc[:,0:2]
        
        #write out with write_csv. 
        final_frame.to_csv('network_edgelist.txt', header = False, index = False, sep = '\t')

    
        

def main():
    
    #handle input arguments
    drug_path = sys.argv[1] 
    targets_path = sys.argv[2]
    nodes_path = sys.argv[3]
    
    #call functions that do necessary calculations and write our output files. 
    my_networkgen = networkgen(drug_path, targets_path, nodes_path)
    pvalues = my_networkgen.compute_pvalues()
    my_networkgen.write_output(pvalues)

if __name__ == "__main__": 
    main()