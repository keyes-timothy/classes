#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 12:41:08 2019

@author: tkeyes
"""

args = parser.parse_args()


######reading in fragmentSet

fragfile = 'helix_9mers.frag'
rmsdfile = 'helix_9mers.rmsd'

my_file = open(fragfile)
my_line = my_file.readline().strip().split()


self = FragmentSet(fragfile, rmsdfile)
self.get_lowRMS_fragments(pos = 9, N = 1)


#####testing FragmentSampler class

my_protein = Protein(sequence = my_sequence)
fragment_set = FragmentSet(fragfile, rmsdfile)

my_frag_sampler = MCMCSampler(my_protein, fragment_set, nmers = 9, start_temp = 100, end_temp = 1, nfrags = 3)

result = my_frag_sampler.simulate(start_temp = 100, end_temp = 1) #doesn't seem to go down that much after the first step...

my_9mer_helix_sampler = MCMCSampler(my_protein, fragment_set, nmers = 9, start_temp = 100, end_temp = 1, nfrags = 3)
result = my_9mer_helix_sampler.simulate(start_temp = 100, end_temp = 1)


my_9mer_helix_sampler = MCMCSampler(my_protein, fragment_set_9mer, nmers = 9, start_temp = 100, end_temp = 1, nfrags = 3)
result = my_9mer_helix_sampler.simulate(start_temp = 100, end_temp = 1)










