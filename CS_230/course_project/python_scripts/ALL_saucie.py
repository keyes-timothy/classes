#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 18:58:17 2020

@author: tkeyes
"""

# Parameters 


# ============================================================================

from model import SAUCIE
from loader import Loader
import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt

tf.reset_default_graph()


# read in my data
x = pd.read_csv(
    "/Users/tkeyes/GitHub/classes/CS_230/course_project/data/saucie_data.csv", 
    sep = ",", 
    header = 0
    )
x = x.to_numpy(dtype = 'float64')


lambda_cs = [i/10 for i in range(0, 1)]
lambda_ds = [i/10 for i in range(0, 10)]

for lambda_c in lambda_cs: 
    print("lambda_c = " + str(lambda_c))
    for lambda_d in lambda_ds:
        print("lambda_d = " + str(lambda_d))
        tf.reset_default_graph()
        # Construct and train SAUCIE model
        my_loader = Loader(x, shuffle=False)
        my_saucie = SAUCIE(x.shape[1], lambda_c=0.01, lambda_d=lambda_d)
        my_saucie.train(load = my_loader, steps = 100, batch_size = 3000)
        
        #extract features from SAUCIE
        embedding = my_saucie.get_embedding(my_loader)
        num_clusters, clusters = my_saucie.get_clusters(my_loader)
        reconstruction = my_saucie.get_reconstruction(my_loader)
        
        # save files
        output_frame = pd.DataFrame({"clusters": clusters, "embedding_1": embedding[:,0], "embedding_2": embedding[:,1]})
        output_frame.to_csv("/Users/tkeyes/GitHub/classes/CS_230/course_project/data/saucie_output_" + str(lambda_c) + "_" + str(lambda_d) + ".csv")
        
        reconstruction_frame = pd.DataFrame(reconstruction)
        reconstruction_frame.to_csv("/Users/tkeyes/GitHub/classes/CS_230/course_project/data/saucie_reconstruction_" + str(lambda_c) + "_" + str(lambda_d) + ".csv")

# cluster_counts = []

# for i in range(num_clusters): 
#     cluster_counts.append(sum(clusters == i))
# cluster_counts


# # save files
# output_frame = pd.DataFrame({"clusters": clusters, "embedding_1": embedding[:,0], "embedding_2": embedding[:,1]})
# output_frame.to_csv("/Users/tkeyes/GitHub/classes/CS_230/course_project/data/saucie_output_new.csv")

# reconstruction_frame = pd.DataFrame(reconstruction)
# reconstruction_frame.to_csv("/Users/tkeyes/GitHub/classes/CS_230/course_project/data/saucie_reconstruction_new.csv")


