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
from tanimoto import Tanimoto #may need to change

class pvalue: 
    
    def __init__(self, num_iterations, seed, drug_path, target_path, protein_a, protein_b): 
        #store paths 
        self.drug_path = drug_path
        self.target_path = target_path
        
        #store names
        self.protein_a = protein_a
        self.protein_b = protein_b
        
        #store information about the bootstrapping procedure
        self.num_iterations = num_iterations
        self.seed = seed
        
        #initialize variables that will be used later
        self.tanimoto = Tanimoto(drug_path, target_path, None)
    
    def precompute_all_tanimotos(self): 
        """
        Precompute all tanimoto coefficients for all drugs in drugs.csv
        
        Inputs: 
            None
        Returns: 
            None
        Side-effects: 
            FILL IN
        """
        self.tanimoto.load_data()
        self.tanimoto.calc_all_tanimotos()
        
    def subset_tanimotos(self, rows, cols):
        """
        After all tanimoto values in the dataset are computed and stored in a symmetric matrix, this function will return a 
        subset of rows and columns of that matrix. These rows/columns represent the tanimoto values between a subset of the 
        drugs in drugs.csv. 
        
        Inputs: 
            rows = indices of the rows to be included in the subset
            cols = indices of the columns to be included in the subset
        Returns: 
            tanimoto_subset = 2D numpy array representing the subsetted data
        """
        tanimoto_subset = self.tanimoto.tanimoto_matrix[np.ix_(rows, cols)] #check that this works
        return tanimoto_subset
        
        
    def calc_T_summary(self, my_tanimotos): 
        """
        Calculates the T-summary statistic for an input matrix of tanimoto scores in which the i,jth entry in the matrix 
        represents the tanimoto score between the ith drug in the first protein's ligand set and the jth drug in the second 
        protein's ligand set. 
        """
        T_summary = 0
        for i in range(my_tanimotos.shape[0]): 
            for j in range(my_tanimotos.shape[1]): 
                current_score = my_tanimotos[i, j]
                if current_score >= 0.5: 
                    T_summary = T_summary + current_score
        return T_summary
        
    #calculate a single bootstrap
    def single_boot(self, ligand_size_a, ligand_size_b): 
        """
        calculates a single bootstrap sample and returns the T-summary value for that step
        
        Inputs: 
            ligand_size_a = an integer indicating the number of ligands in the first protein's ligand set
            ligand_size_b = an integer indicating the number of ligands in the second protein's ligand set
            
        Returns: 
            the T-summary value for the current bootstrapped step
        """
        #select rows and columns 
        rows = np.array([random.randint(0, len(self.tanimoto.drug_data)-1) for i in range(ligand_size_a)])
        cols = np.array([random.randint(0, len(self.tanimoto.drug_data)-1) for i in range(ligand_size_b)])
        
        #compute T-summary statistic
        current_T = self.calc_T_summary(my_tanimotos = self.subset_tanimotos(rows, cols))
        
        #return the T-summary statistic
        return current_T
        
    
    #calculate all bootstraps within the bootstrapping procedure
    
    def all_bootstrapping(self, ligand_size_a, ligand_size_b):
        """
        Calculates all bootstrapping steps and storing the T_summaries obtains with each. 
        
        Inputs: 
            FILL IN 
        Returns: 
            FILL IN 
        """
        bootstrapped_T_summaries = []
        for i in range(self.num_iterations): 
            new_T_summary = self.single_boot(ligand_size_a, ligand_size_b)
            bootstrapped_T_summaries.append(new_T_summary)
        
        return bootstrapped_T_summaries
    
    
    def calc_p_value(self, observed_T_summary, bootstrapped_T_summaries): 
        """
        Documentation
        
        Inputs: 
            FILL IN 
        Returns: 
            FILL IN 
        """
        p = sum(observed_T_summary<pd.Series(bootstrapped_T_summaries))/self.num_iterations
        return p
        

def main():
    
    #handle input arguments
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-n', type = int, dest = 'num_iterations', help = 'number of iterations', default = 500)
    parser.add_argument('-r', type = int, dest = 'seed', help = 'parameter that sets the state of the pseudo-random number generator', default = 214) 
    parser.add_argument('drug_path', help = 'path to drugs.csv file')
    parser.add_argument('target_path', help = 'path to target.csv file')
    parser.add_argument('uniprot_1', help = 'first uniprot accession ID')
    parser.add_argument('uniprot_2', help = 'second uniprot accession ID')

    args = parser.parse_args()
    
    #set up pvalue object that will do the calculations
    my_pvalue = pvalue(args.num_iterations, args.seed, args.drug_path, args.target_path, args.uniprot_1, args.uniprot_2)
    my_pvalue.precompute_all_tanimotos()
    
    #find observed T-summary value for our two proteins of interest
    protein_a_drug_names = my_pvalue.tanimoto.protein_data[my_pvalue.protein_a]
    protein_b_drug_names = my_pvalue.tanimoto.protein_data[my_pvalue.protein_b]
    
    #find the subsetted matrix for protein_a and protein_b of interest and calculate observed T-summary
    protein_a_indices = [i for i, drug in enumerate(list(my_pvalue.tanimoto.drug_data.keys())) if drug in protein_a_drug_names]
    protein_b_indices = [i for i, drug in enumerate(list(my_pvalue.tanimoto.drug_data.keys())) if drug in protein_b_drug_names]
    
    my_protein_tanimotos = my_pvalue.subset_tanimotos(rows = protein_a_indices, cols = protein_b_indices)
    observed_T_summary = my_pvalue.calc_T_summary(my_protein_tanimotos)
    
    #perform bootstrapping and find p-value
    bootstrapped_Ts = my_pvalue.all_bootstrapping(len(protein_a_indices), len(protein_b_indices))
    p_value = my_pvalue.calc_p_value(observed_T_summary, bootstrapped_Ts)
    print(p_value)
    

if __name__ == "__main__": 
    main()