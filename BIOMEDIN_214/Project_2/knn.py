# -*- coding: utf-8 -*-
"""
File name: knn.py

Author: tkeyes

Last Updated: 10/23/19

Description: Creates a class that runs the KNN classifier algorithm. 
"""

import sys
import pandas as pd
import numpy as np

##global variables 
k = 3
fn = 0.5


#Define my class

class KNN: 
    """
    Documentation
    """
    def __init__(self):
        self.expression_matrix = pd.DataFrame()
        self.metadata = pd.DataFrame()
        self.my_assignments = {}
    
    def load_data(self, expfile, sampfile): 
        """
        A function that stores the input expression and sample data in the KNN object. 
        
        Inputs: 
            expfile = a string indicating the path to the expression file 
            sampfile = a string indicating the path to the sample file
            
        Returns: 
            Nothing
        
        Side-effects: 
            Stores the information in expfile and sampfile within the KNN object. 
        """
        self.expression_matrix = pd.read_csv(expfile, sep = '\t')
        self.metadata = pd.read_csv(sampfile, sep = '\t', header = None, names = ['patient', 'condition'])
        
        #make expression_matrix and metadata slightly easier to work with
        self.expression_matrix.index = self.expression_matrix.SYMBOL
        self.expression_matrix.drop('SYMBOL', axis = 1, inplace = True)
        self.metadata.index = self.metadata.patient
        self.metadata.drop('patient', axis = 1, inplace = True)
    
    def get_assignments(self, k, fn): 
#        self.expression_matrix.index = self.expression_matrix.SYMBOL
#        self.expression_matrix.drop('SYMBOL', axis = 1, inplace = True)
#        self.metadata.index = self.metadata.patient
#        self.metadata.drop('patient', axis = 1, inplace = True)
        
        #for each patient
        for i in range(0, len(self.expression_matrix.columns)): 
            patient_name = self.expression_matrix.columns[i]
            
            #remove it from the dataset
            my_patient = self.expression_matrix.iloc[:,i]
            my_data = self.expression_matrix.drop(patient_name, axis = 1)
            
            #find its k closest neighbors
            my_diff = my_data.sub(my_patient, axis = 'index') #subtract my_patient data from each of the other patients
            my_diff_2 = pow(my_diff, 2)
            my_distances = pow(my_diff_2.sum(axis = 'index'), 0.5)
            nearest_neighbor_names = my_distances.sort_values(ascending = True).index.values[0:k]
            
            #query what those neighbors' conditions are in self.metadata
            average_class = self.metadata.loc[nearest_neighbor_names].values.mean()
        
            #return 1 if average class rating is > fn, if not return 0.
            if average_class > fn: 
                self.my_assignments[patient_name] = 1
            else: 
                self.my_assignments[patient_name] = 0
        
        my_predictions = []
        for patient in self.metadata.index.values: 
            print(patient)
            my_predictions.append(self.my_assignments[patient])
        
        return my_predictions
        
    def calc_metrics(self, k, fn):
        true_values = self.metadata.condition.values
        my_predictions = self.get_assignments(k, fn)
        TPs = 0
        FPs = 0
        TNs = 0
        FNs = 0
        for i in range(len(true_values)): 
            if true_values[i] == 0 and my_predictions[i] == 0: 
                TNs +=1
            elif true_values[i] == 0 and my_predictions[i] == 1: 
                FPs +=1
            elif true_values[i] == 1 and my_predictions[i] == 1: 
                TPs +=1
            elif true_values[i] == 1 and my_predictions[i] == 0: 
                FNs +=1
        
        sensitivity = TPs/(TPs + FNs)
        specificity = TNs/(TNs + FPs)
        
        return sensitivity, specificity
    

#run the main method

def main(): 
    
    #check that the file is being properly used
    None
    
    #input variables
    expfile = sys.argv[1]
    sampfile = sys.argv[2]
    
    #create a KNN object and run
    my_knn = KNN(expfile, sampfile)
    my_knn.load_data(expfile, sampfile)
    my_knn.get_assignments(k, fn)
    my_knn.calc_metrics(k, fn)

if __name__ == "__main__": 
    main()

    