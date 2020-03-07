#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 12:06:31 2019

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
from networkgen import networkgen
import networkx
import matplotlib.pyplot as plt


class plot_graph: 
    
    def __init__(self, edgelist_path, nodes_path, out_path): 
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
        self.graph = networkx.read_edgelist(edgelist_path)
        self.nodes_df = pd.read_csv(nodes_path)
        self.out_path = out_path
        

    def plot_network(self): 
        """
        A fuction that plots our protein network according to the specifications of the project page. 
        
        Inputs: 
            None
            
        Returns: 
            Nothing
            
        Side-effects:
            Saves a .png file of the network plotted using networkx.draw_networkx(). This image represents each protein as a node, 
            colored by indication, and with nodes connected if the associated p-value between them is less than 0.05.
            
        """
        #relabel
        name_map = {}
        for i in range(len(self.nodes_df.uniprot_accession)): 
            code = self.nodes_df.uniprot_accession[i]
            name = self.nodes_df.uniprot_id.values[i]
            name_map[code] = name
        self.graph = networkx.relabel_nodes(G = self.graph, mapping = name_map, copy = True)
        
        #color
        colors = []
        color_map = {"bp":"red", "bp;cholesterol":"green", "bp;cholesterol;diabetes":"blue", "bp;diabetes":"purple"}
        for node in self.graph:
            node_type = self.nodes_df.loc[self.nodes_df.uniprot_id == node].iloc[:,2].values[0]
            node_color = color_map[node_type]
            colors.append(node_color)
        
        plt.figure(figsize=(8,8)) 
        networkx.draw_networkx(G = self.graph, node_color = colors)
        plt.savefig(self.out_path, format = 'PNG', dpi = 150)
        

        

def main():
    
    #handle input arguments
    edgelist_path = sys.argv[1] 
    nodes_path = sys.argv[2]
    out_path = sys.argv[3]
    
    #plot the network
    my_plot = plot_graph(edgelist_path, nodes_path, out_path)
    my_plot.plot_network()

if __name__ == "__main__": 
    main()