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
    
    def __init__(self, num_iterations, seed, drug_path, target_path): 
        """
        Initialization function that saves relevant file paths and initializes instance variables for the pvalue object.
        
        Inputs: 
            num_interations = number of iterations to use for bootstrapping.
            seed = integer representing the random seed to be used for the bootstrapping.
            drugs_path = path to drugs.csv
            targets_path = path to targets.csv
        
        Returns: 
            Nothing
        
        Side-effects: 
            initializes several instance variables that store the paths to each input file (as well as the the eventual output), 
            instance variables that store the random seed and number of iterations to use for the boostrapping, and sets up a 
            Tanimoto object to perform the tanimoto calculations.
        """
        #set seed
        np.random.seed(seed)
        
        #store paths
        self.drug_path = drug_path
        self.target_path = target_path
        
        #store information about the bootstrapping procedure
        self.num_iterations = num_iterations
        self.seed = seed
        
        #initialize variables that will be used later
        self.tanimoto = Tanimoto(drug_path, target_path, None)
    
    def precompute_all_tanimotos(self): 
        """
        A function that precomputes all tanimoto coefficients for all drugs in drugs.csv and loads them into the pvalue object.
        
        Inputs: 
            None
            
        Returns: 
            None
            
        Side-effects: 
            Calls the load_data() and calc_all_tanimotos methods of self.tanimoto (a Tanimoto object). See tanimoto.py for details.
        """
        self.tanimoto.load_data()
        self.tanimoto.calc_all_tanimotos()
        
    def subset_tanimotos(self, protein_a, protein_b):
        """
        After all tanimoto values in the dataset are computed and stored in a symmetric matrix, this function will return a 
        subset of rows and columns of that matrix. These rows/columns represent the tanimoto values between a subset of the 
        drugs in drugs.csv (only the drugs that bind to protein_a and the drugs that bind to protein_b). 
        
        Inputs: 
            protein_a = character string representing the ID of the first protein
            protein_b = character string representing the ID of the second protein
            
        Returns: 
            tanimoto_subset = 2D numpy array representing the subsetted data
        """
        
        #look up the names of the drugs that bind to protein_a and protein_b
        protein_a_drug_names = self.tanimoto.protein_data[protein_a]
        protein_b_drug_names = self.tanimoto.protein_data[protein_b]
    
        #find the subsetted matrix for protein_a and protein_b of interest and calculate observed T-summary
        protein_a_indices = [i for i, drug in enumerate(list(self.tanimoto.drug_data.keys())) if drug in protein_a_drug_names]
        protein_b_indices = [i for i, drug in enumerate(list(self.tanimoto.drug_data.keys())) if drug in protein_b_drug_names]
        
        #return the subset
        tanimoto_subset = self.tanimoto.tanimoto_matrix[np.ix_(protein_a_indices, protein_b_indices)]
        return tanimoto_subset
        
        
    def calc_T_summary(self, my_tanimotos): 
        """
        A function that calculates the T-summary statistic for an input matrix of tanimoto scores in which the i,jth entry in the matrix 
        represents the tanimoto score between the ith drug in the first protein's ligand set and the jth drug in the second 
        protein's ligand set. 
        
        Inputs: 
            my_tanimotos = a 2D NumPy array containing the tanimoto scores for drugs that bind to either of our proteins of interest. 
                           This subset will be generated by a pvalue object's method self.subset_tanimotos. 
                           
        Returns: 
            T_summary = a T_summary statistic for the proteins whose ligands' pairwise Tanimoto scores are represented in my_tanimotos. 
        """
        T_summary = 0
        for i in range(my_tanimotos.shape[0]): 
            for j in range(my_tanimotos.shape[1]): 
                current_score = my_tanimotos[i, j]
                if current_score > 0.5: 
                    T_summary = T_summary + current_score
        return T_summary
        
    def single_boot(self, ligand_size_a, ligand_size_b): 
        """
        A function that calculates a single bootstrap step and returns the T-summary value calculated for that step. Sampling of 
        random drugs is performed by picking random row and column indices from self.tanimoto.tanimoto_matrix (the matrix in which 
        all pairwise tanimoto scores for all drugs in drugs.csv are stored).
        
        Inputs: 
            ligand_size_a = an integer indicating the number of ligands in the first protein's ligand set
            ligand_size_b = an integer indicating the number of ligands in the second protein's ligand set
            
        Returns: 
            the T-summary value for the current bootstrapped step
        """
        #select rows and columns to sample from self.tanimoto.tanimoto_matrix by generating random integers.
        rows = np.array([random.randint(0, len(self.tanimoto.drug_data)-1) for i in range(ligand_size_a)])
        cols = np.array([random.randint(0, len(self.tanimoto.drug_data)-1) for i in range(ligand_size_b)])
        
        #compute and return T-summary statistic
        current_T = self.calc_T_summary(my_tanimotos = self.tanimoto.tanimoto_matrix[np.ix_(rows, cols)])
        return current_T
            
    def all_bootstrapping(self, ligand_size_a, ligand_size_b):
        """
        A function that performs all bootstrapping steps and stores the T_summaries obtained in each step. 
        
        Inputs:
            ligand_size_a = an integer indicating the number of ligands in the first protein's ligand set
            ligand_size_b = an integer indicating the number of ligands in the second protein's ligand set
            
        Returns: 
            bootstrapped_T_summaries = a list of floating decimal values. Each entry represents the bootstrapped T_summary 
                                       statistic obtained from a single iteration of the bootstrapping. 
        """
        bootstrapped_T_summaries = []
        for i in range(self.num_iterations): 
            new_T_summary = self.single_boot(ligand_size_a, ligand_size_b)
            bootstrapped_T_summaries.append(new_T_summary)
        
        return bootstrapped_T_summaries
    
    
    def calc_p_value(self, observed_T_summary, bootstrapped_T_summaries): 
        """
        A function that calculates the p-value associated with a pair of proteins after bootstrapping is performed. 
        
        Inputs: 
            observed_T_summary = a floating decimal value indicating the actual, observed T_summary statistic for a pair of proteins.
            bootstrapped_T_summaries = a list of bootstrapped T_summary statistics obtained via random sampling. This list can be 
                                       obtained by running the self.all_bootstrapping method. 
        Returns: 
            p = the p-value associated with a pair or proteins.
        """
        p = sum(observed_T_summary<pd.Series(bootstrapped_T_summaries))/self.num_iterations
        return p
    
    def run_full_analysis(self,protein_a, protein_b):
        """
        A function that runs the entire analysis to associate the binding patterns of two proteins. Performs all necessary calculations and 
        outputs a p-value for the association (calculated via our bootstrapping procedure).
        
        Inputs: 
            protein_a = ID for the first protein
            protein_b = ID for the second protein
        Returns: 
            p_value = p_value for the association between the two proteins. 
        """
        #find observed T-summary value for our two proteins of interest
        my_protein_tanimotos = self.subset_tanimotos(protein_a, protein_b)
        observed_T_summary = self.calc_T_summary(my_protein_tanimotos)
    
        #perform bootstrapping and find p-value
        bootstrapped_Ts = self.all_bootstrapping(my_protein_tanimotos.shape[0], my_protein_tanimotos.shape[1])
        p_value = self.calc_p_value(observed_T_summary, bootstrapped_Ts)
        return p_value
        

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
    my_pvalue = pvalue(args.num_iterations, args.seed, args.drug_path, args.target_path)
    my_pvalue.precompute_all_tanimotos()
    
    #run full analysis and print result
    print(my_pvalue.run_full_analysis(args.uniprot_1, args.uniprot_2))
    

if __name__ == "__main__": 
    main()