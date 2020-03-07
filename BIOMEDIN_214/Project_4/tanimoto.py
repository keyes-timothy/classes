#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 13:28:01 2019

@author: tkeyes
"""
import sys
import argparse
import numpy as np
import os
import pandas as pd

class Tanimoto: 
    """
    Class that performs pairwise tanimoto score calculations for all pairs of drugs in an input file. Will also find if each drug pair
    shares a common target and write an output file storing this information. 
    """
    
    def __init__(self, drugs_path, targets_path, out_path): 
        """
        Initialization function that saves the relevant file paths to be used by later functions. 
        
        Inputs: 
            drugs_path = path to drugs.csv
            targets_path = path to targets.csv
            out_path = path to output file
        
        Returns: 
            Nothing
        
        Side-effects: 
            initializes several instance variables that store the paths to each input file (as well as the the eventual output), 
            dictionaries for fast lookup of which drugs bind to which proteins (and which proteins bind to which drugs), and a 2D array for 
            storing tanimoto scores between all drugs in the input drugs.csv file. 
        """
        #store paths
        self.drugs_path = drugs_path
        self.targets_path = targets_path
        self.out_path = out_path
        
        #initialize input data structures
        self.drug_data = {}
        self.target_data = {}
        
        #initialize way of storing tanimoto scores
        self.tanimoto_matrix = np.empty(dtype = float, shape = (len(self.drug_data), len(self.drug_data)))
        self.target_matrix = np.empty(dtype = int, shape = (len(self.drug_data), len(self.drug_data))) 
    
    def load_data(self):
        """
        Loads data from .csvs containing information about each drug's chemical features and their target proteins. 
        
        Inputs: 
            None.
        
        Outputs: 
            self.drug_data = a dictionary in which keys are drug IDs and values are sets of the keys' chemical features. 
                             Saved as a global variable within each class instance. 
            self.target_data = a dictionary in which keys are drug IDs and values are sets of the keys' target proteins. 
                             Saved as a global variable within each class instance. 
            self.protein_data = a dictionary in which keys are protein IDs and values are sets of the proteins' ligands. 
                                Saved as a global variable within each class instance.
        """
        #set up data structure for data in drugs.csv
        drug_dataFrame = pd.read_csv(self.drugs_path)
        featureSets = []
        for i in range(len(drug_dataFrame.index)): 
            featureSet = set(drug_dataFrame.maccs.values[i].strip().split())
            featureSets.append(featureSet)
        
        self.drug_data = pd.Series(featureSets, index = drug_dataFrame.db_id.values).to_dict()  #this is a dictionary in which the drug id is the key and the set of features is the value.
        
        #set up data stucture for data in targets.csv
        targets_dataFrame = pd.read_csv(self.targets_path)
        
        drug_set = set(drug_dataFrame.db_id.values)
        protein_set = set(targets_dataFrame.uniprot_accession.values)
        
        target_data = pd.Series([set() for index in range(len(drug_set))], index = list(drug_set))
        protein_data = pd.Series([set() for index in range(len(protein_set))], index = list(protein_set))
        
        for i in range(len(targets_dataFrame.db_id.values)): 
            drug = targets_dataFrame.db_id.values[i]
            drug_target = targets_dataFrame.uniprot_accession.values[i]
            
            #add drug target to set corresponding to that drug in target_data
            target_data[drug].add(drug_target)
            
            #add drug to ligand set corresponding to its protein in protein_data
            protein_data[drug_target].add(drug)
            
        self.target_data = target_data  #this is a dictionary in which the drug id is the key and the set of target proteins is the value.
        self.protein_data = protein_data #this is a dictionary in which the protein id is the key and the set of ligand drugs is the value.
    
    def calc_tanimoto(self, drug_a, drug_b): 
        """
        Calculates the tanimoto score between two drugs. 
        
        Inputs: 
            drug_a = string encoding the database id for the first drug. 
            drug_b = string encoding the database id for the second drug. 
            
        Returns: 
            Tanimoto score for those two proteins
        """
        set_a = self.drug_data[drug_a]
        set_b = self.drug_data[drug_b] 
        
        my_tanimoto = len(set_a.intersection(set_b))/len(set_a.union(set_b))
        
        return my_tanimoto
    
    def has_common_target(self, drug_a, drug_b): 
        """
        Figures out if two drugs have a common target. 
        
        Inputs: 
            drug_a = the database name for the first drug
            drug_b = the database name for the second drug
            
        Returns: 
            0 if the drugs have no common targets 
            1 if the drugs have 1 or more common targets. 
        """
        set_a = self.target_data[drug_a]
        set_b = self.target_data[drug_b] 
        
        num_common_targets = len(set_a.intersection(set_b))
        
        if num_common_targets == 0: 
            return 0
        else: 
            return 1
        
    def calc_all_tanimotos(self): 
        """
        Calculates all pairwise tanimoto coefficient values between drugs in the original drugs.csv dataset.
        
        Inputs: 
            None
            
        Returns: 
            None
            
        Side-effects: 
            Updates self.tanimoto_matrix = 2D numpy array that stores the pairwise tanimoto scores of each pair of drugs by 
                calling self.calc_tanimoto on each possible pair. 
            
        """
        #reshape storage structures for tanimoto and common target calculation results
        self.tanimoto_matrix = np.empty(dtype = float, shape = (len(self.drug_data), len(self.drug_data)))
        #self.target_matrix = np.empty(dtype = int, shape = (len(self.drug_data), len(self.drug_data))) 
        
        for i in range(len(self.drug_data)): #for each row in tanimoto matrix
            for j in range(len(self.drug_data)): #for each column in tanimoto matrix
                drug_a = list(self.drug_data.keys())[i]
                drug_b = list(self.drug_data.keys())[j]
                
                #calculate tanimoto score
                self.tanimoto_matrix[i,j] = self.calc_tanimoto(drug_a, drug_b)
                
        
    def find_all_shared_targets(self):
        """
        A function that calculates whether the pairs of drugs in drugs.csv share protein binding partners. Calls self.has_common_target
        for each possible drug pair. 
        
        Inputs: 
            None
            
        Returns: 
            None
            
        Side-effects: 
            Updates self.target_matrix = 2D numpy array that stores a 0 at its i,jth location if the ith and jth drugs do not 
                have any common targets and a 1 otherwise. 
        """
        self.target_matrix = np.empty(dtype = int, shape = (len(self.drug_data), len(self.drug_data))) 
        
        for i in range(len(self.drug_data)): #for each row in tanimoto matrix
            for j in range(len(self.drug_data)): #for each column in tanimoto matrix
                drug_a = list(self.drug_data.keys())[i]
                drug_b = list(self.drug_data.keys())[j]
                self.target_matrix[i,j] = self.has_common_target(drug_a, drug_b)
                    
    
    def write_output(self): 
        """
        A function that writes the requested output file (a .csv file containing the tanimoto scores for all pairs of drugs in drugs.csv, 
        as well as a 0 or 1 indicating if each pair shares a protein binding target) for tanimoto.py. 
        
        Inputs: 
            None
            
        Returns: 
            None
            
        Side-effects: 
            writes a .csv file to self.out_path with 4 columns: 
                Column 1 = First drug's name
                Column 2 = Second drug's name
                Column 3 = Tanimoto score for the pair of drugs
                Column 4 = 0 if the drugs do not share a common protein target, 1 if the drugs do share a common protein target.
        """
        #obtain only unique pairs
        unique_pairs = set()
        for drug_a in self.drug_data.keys():
            for drug_b in self.drug_data.keys():
                if drug_a != drug_b: 
                    unique_pairs.add(frozenset([drug_a, drug_b]))
        unique_pairs = list(unique_pairs)
        
        #switch to pandas dataFrames for indexing convenience, which are very slow but that should not matter much in this final step.
        tanimoto_frame = pd.DataFrame(self.tanimoto_matrix, columns = list(self.drug_data.keys()), index = self.drug_data.keys())
        target_frame = pd.DataFrame(self.target_matrix, columns = list(self.drug_data.keys()), index = self.drug_data.keys())
        
        #collect necessary values for each unique pair to be put in the output csv 
        drug_a_list = []
        drug_b_list = []
        tanimoto_list = []
        target_list = []
        for pair in unique_pairs: 
            drug_a = list(pair)[0]
            drug_b = list(pair)[1]
            
            drug_a_list.append(drug_a)
            drug_b_list.append(drug_b)
            tanimoto_list.append(round(tanimoto_frame.loc[drug_a, drug_b], ndigits = 6))
            target_list.append(target_frame.loc[drug_a, drug_b])
            
        #read out final csv
        output_frame = pd.DataFrame(list(zip(drug_a_list, drug_b_list, tanimoto_list, target_list)))
        output_frame.to_csv(self.out_path, header = False, index = False)
        
        

def main():
    
    drugs_path = sys.argv[1] 
    targets_path = sys.argv[2]
    out_path = sys.argv[3]
    
    my_tanimoto = Tanimoto(drugs_path = drugs_path, targets_path = targets_path, out_path = out_path)
    my_tanimoto.load_data()
    my_tanimoto.calc_all_tanimotos()
    my_tanimoto.find_all_shared_targets()
    my_tanimoto.write_output()

if __name__ == "__main__": 
    main()