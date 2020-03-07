#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File name: knn.py

Author: tkeyes

Last Updated: 10/23/19

Description: Creates a class that runs the KNN classifier algorithm. 
"""

#import modules
import sys
import pandas as pd
import numpy as np

class GSEA: 
    """
    Documentation
    """
    def __init__(self):
        """
        Documentation
        """      
        #initialize data structures to be used later - decide how to handle these later
        self.expression_matrix = pd.DataFrame()
        self.metadata = pd.DataFrame()
        self.genesets = {}
        self.gene_list = []
        
    def load_data(self, expfile, sampfile, genesets): 
        """
        Documentation
        """
        #read in expression data 
        self.expression_matrix = pd.read_csv(expfile, sep = '\t')
        
        #read in patient metadata
        self.metadata = pd.read_csv(sampfile, sep = '\t', header = None, names = ['patient', 'condition'])
        
        #read in gene sets
        my_file = open(genesets) 
        for line in my_file: 
            line = line.strip().split()
            del(line[1])
            self.genesets[line[0]] = line[1:]
        my_file.close()
    
    def get_fold_changes(self, patient_names, control_names): 
        """
        Function to calculate the FC between mean patient and mean control expression values. Called by get_gene_rank_order(). 
        """
        
        #set up two dataFrames for patient and control samples 
        patient_frame = self.expression_matrix[patient_names]
        patient_frame.index = self.expression_matrix['SYMBOL']
        healthy_frame = self.expression_matrix[control_names]
        healthy_frame.index = self.expression_matrix['SYMBOL']
        
        #calculate means for patients and controls across individuals
        patient_means = patient_frame.mean(axis = 1)
        healthy_means = healthy_frame.mean(axis = 1)
        
        FCs = patient_means - healthy_means
        
        return FCs.sort_values(ascending = False)
    
    def get_gene_rank_order_helper(self, patient_names, control_names): 
        """
        Documentation
        
        Can be used with permutations.
        """
        fold_changes = self.get_fold_changes(patient_names, control_names)
        gene_list = list(fold_changes.index)
        return gene_list
    
    def get_gene_rank_order(self):
        """
        Documentation
        
        Will return the "true" gene rank order by calling get_gene_rank_order_helper(). 
        """
        patient_names = self.metadata.loc[self.metadata['condition']== 1]['patient'].values
        control_names = self.metadata.loc[self.metadata['condition']== 0]['patient'].values
        
        return self.get_gene_rank_order_helper(patient_names, control_names)
        
    
    def get_enrichment_score_helper(self, geneset, ranked_gene_list): #LIKELY WHERE THE PROBLEM IS
        """
        Documentation
        
        Can be used with permutations (in which the ranked_gene_list will differ)
        """
        set_genes = self.genesets[geneset]
        set_genes_in_dataset = list(set(set_genes) - (set(set_genes) - set(self.expression_matrix.SYMBOL)))
        
        num_total_genes = len(self.expression_matrix.index)
        num_set_genes = len(set_genes_in_dataset)
        
        hit_points = pow((num_total_genes - num_set_genes)/num_set_genes, 0.5)
        miss_points = pow((num_set_genes/(num_total_genes - num_set_genes)), 0.5)
                
        #iterate through every gene in ranked_gene_list
        scores = []
        current_score = 0
        for gene in ranked_gene_list: #where ranked_gene_list is the sorted list of all genes according to log2FC
            if gene in set_genes_in_dataset: 
                new_score = current_score + hit_points
            else: 
                new_score = current_score - miss_points
            scores.append(new_score)
            current_score = new_score
        
        return round(max(scores),2)
    
    
    def get_enrichment_score(self, geneset):
        """
        Function that returns the ES of a geneset given a ranked list of differentially expressed genes between patients and controls. 
        
        Inputs: 
            FILL IN
        Returns: 
            FILL IN
            
        Will return the "true" observed values.
        """
        ranked_gene_list = self.get_gene_rank_order()
        return self.get_enrichment_score_helper(geneset, ranked_gene_list)
    

    def get_all_enrichment_scores(self, patient_names, control_names):
        """
        Function that returns a pandas Series of all genesets' enrichment scores using get_enrichment_score_helper()
        
        Input: 
            patient_names = array of strings indicating which columns in self.expression_matrix are patients 
            control_names = array of strings indicating which columns in self.expression_matrix are controls 
        """
        #define ranked_gene_list
        ranked_gene_list = self.get_gene_rank_order_helper(patient_names, control_names)
        
        enrichment_scores = {}
        for geneset in self.genesets: 
            current_score = self.get_enrichment_score_helper(geneset, ranked_gene_list)
            enrichment_scores[geneset] = current_score
        enrichment_scores = pd.Series(enrichment_scores)
        enrichment_scores.sort_values(ascending = False, inplace = True)
        return enrichment_scores
    
    def write_output_files(self, observed_enrichment_scores): 
        """
        Documentation
        """
        observed_enrichment_scores.to_csv('kegg_enrichment_scores.txt', sep = '\t', header = False)
        
    def get_sig_sets(self, p): 
        """
        Documentation
        """
        real_patients = self.metadata.loc[self.metadata['condition']== 1]['patient'].values
        real_controls = self.metadata.loc[self.metadata['condition']== 0]['patient'].values
        observed_enrichment_scores = self.get_all_enrichment_scores(real_patients, real_controls)
        self.write_output_files(observed_enrichment_scores)
        
        #document some information to use for the permutations
        all_sample_names = list(real_patients) + list(real_controls)
        num_patients = len(real_patients)
        num_controls = len(real_controls)
        
        #initialize a dataframe that will hold the results of each iteration
        simulated_ESs = pd.DataFrame(index = observed_enrichment_scores.index, columns = range(100), data = 0)
        
        #do 100 interations
        for i in range(100): 
            #randomly assign samles to patients and controls 
            permuted_sample_names = np.random.permutation(all_sample_names)
            permuted_patients = permuted_sample_names[0:num_patients]
            permuted_controls = permuted_sample_names[num_patients:]
            
            #call self.get_all_enrichment_scores()
            permuted_ESs = self.get_all_enrichment_scores(patient_names = permuted_patients, control_names = permuted_controls)
            
            #save that as the ith column of simulated_ESs
            simulated_ESs[i] = permuted_ESs
        
        p_values = pd.Series(index = observed_enrichment_scores.index, data = -1.0)
        for geneset in p_values.index: 
            p_values[geneset] = sum(simulated_ESs.loc[geneset] >= observed_enrichment_scores.loc[geneset])/len(simulated_ESs.loc[geneset])
        
        self.p_values = p_values.sort_values(ascending = True)
        
        #correct p-values for multiple comparisons
        threshold = p/len(self.genesets)
        return(list(p_values.loc[p_values<threshold].index.values))
        

            
    

def main(): 
    """
    Documentation
    """
    
    #check that the file is being properly used
    None
    
    #input variables
    expfile = sys.argv[1]
    sampfile = sys.argv[2]
    keggfile = sys.argv[3]
    
    #create a GSEA object and run
    my_gsea = GSEA(expfile, sampfile, keggfile)
    my_gsea.load_data(expfile, sampfile, keggfile)
    my_gsea.get_sig_sets(0.5)
    
if __name__ == "__main__": 
    main()
