"""
This file contains the main fragment sampling class, which performs a Monte Carlo simulated annealing procedure to fold a protein.
"""

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.scoring import *
init(extra_options='-mute all  -constant_seed')

import numpy as np
import utils
from Protein import Protein
from FragmentSet import FragmentSet

import os
import random
import math
import pandas as pd


class MCMCSampler(object):
    def __init__(self, protein, fragment_set, nmers, start_temp, end_temp, nfrags, anneal_rate):
        """
        Initialize necessary variables for an MCMCSampler object.
        
        Inputs: 
            protein = a Protein object
            fragment_set = a FragmentSet object
            nmers = an int representing the number of residues in each fragment of the fragment_set
            start_temp = starting temperature for the simulated annealing
            end_temp = ending temperature for the simulated annealing
            nfrags = the number of lowest RMSD fragments that the MCMCSampler will use during sampling
        """
        
        #store variables used to initiate the object
        self.scorefxn = create_score_function('score3')
        self.current_protein = protein
        self.best_protein = protein
        self.my_fragment_set = fragment_set
        self.nmers = nmers
        self.nfrags = nfrags
        self.end_temp = end_temp
        self.anneal_rate = anneal_rate
        
        #initialize some starter values (for tracking later...)
        self.current_energy = self.compute_energy(self.current_protein)
        self.current_iteration = 0
        self.current_T = start_temp

        
        #initialize a data structure to keep track of which fragments (at a given position) have been sampled already during each sampling step
        self.sampled_fragments = {}
        
        #create a dictionary of nfrag candidate fragments for each position in the sequence
        self.candidate_frag_list = {}
        for position in range(1, protein.length-self.nmers+1):
            self.candidate_frag_list[position] = self.my_fragment_set.get_lowRMS_fragments(position, nfrags)
            
        #for reporting information to the log file later
        self.temperature = [self.current_T]
        self.iteration = [self.current_iteration]
        self.energy = [self.current_energy]
        
    def compute_energy(self, protein):
        """
        Compute energy of a protein.
        --------
        Params:
            - protein (Protein object): protein to score
        Return:
            - energy of conformation (float)
        """        
        return utils.score_pose(protein.pose, self.scorefxn)

    def perturb_fragment(self, protein, position): # you may want to add more arguments
        """
        Sample from possible fragments for a position, and replace torsion angles of that fragment in the protein.
        ---------
        Params:
            - protein = an input protein object that you want to copy and perturb
            - position = a position (1-indexed) in the input protein at which you would like to start the perturbation
        Returns:
            - perturbed fragment (a protein object)
        """
        
        #check which fragments have already been sampled from the current position during the current step
        chosen_indices = self.sampled_fragments[position]
        fragments_to_sample_from = set(range(self.nfrags)) - chosen_indices
        
        #choose a candidate fragment at this position, then add that candidate to the list of previously chosen fragments (during this step)
        chosen_candidate = random.choice(list(fragments_to_sample_from))
        self.sampled_fragments[position].add(chosen_candidate)
        
        #after candidate fragment is chosen, make perturbed fragment and return
        new_positions = self.candidate_frag_list[position][chosen_candidate]
        perturbed_fragment = Protein(pose = protein.pose)
        for i in range(len(new_positions)): 
            perturbed_fragment.set_torsion((position + i), new_positions[i][0], new_positions[i][1])
            
        return(perturbed_fragment)
        

    def metropolis_accept(self, old_protein, new_protein, T):
        """
        Calculate probability of accepting or rejecting move based on Metropolis criterion.
        --------
        Params:
            old_protein = Protein object from previous step of the sampling procedure
            new_protein = Protein object from the current step of the sampling procedure (i.e. a candidate protein for the current move)
            T = current annealing temperature at which the metropolis criterion is being calculated
        Returns:
            - p - the probability (float between 0 and 1) of accepting the move based on Metropolis criteria. 
        """
        delta_E = self.compute_energy(new_protein) - self.compute_energy(old_protein)
        
        if delta_E <= 0: 
            p = 1
        else: 
            p = math.e**(-delta_E/T)
        return(p)
            

    def anneal_temp(self, T):
        """
        Anneal temperature using exponential annealing schedule. Consider kT to be a single variable (i.e. ignore Boltzmann constant)
        --------
        Params:
            - T: the temperature at the current step
        Returns:
            - new_T: the new temperature after annealing
        """
        new_T = self.anneal_rate * T
        return(new_T)

    def step(self, T):
        """
        Take a single MCMC step. Each step should do the following:
        1. sample position in chain
            - Note: think about positions you can sample a k-mer fragment from. 
              For example, you cannot sample from position 1 because there is no phi angle
        2. sample fragment at that position and replace torsions in a *copied version* of the protein
        3. measure energy after replacing fragment
        4. accept or reject based on Metropolis criterion
            - if accept: incorporate proposed insertion and anneal temperature
            - if reject: sample new fragment (go to step 3)
            
        Params: 
            T = temperature at which the step is taking place
        Returns: 
            Nothing
        Side-effects: 
            Updates self.current_T, self.best_protein, self.current_energy, and self.current_iteration in the MCMCSampler Object.
        
            
        """
        #reset structure that keeps track of all the fragments you have tried during this step
        #print(self.current_iteration)
        self.sampled_fragments = {}
        self.current_tries = 0
        for position in range(1, self.current_protein.length-self.nmers+1):
            self.sampled_fragments[position] = set()
            
        #sample position in the protein sequence
        sequence_position_range = range(2, len(self.current_protein.sequence)-self.nmers+1)
        
        while True: #keep looping until you accept a move
            
            #chooose a position in the sequence
            position = random.choice(sequence_position_range)
            
            #if you have already tried all the fragments at a particular position, choose a different position
            while len(self.sampled_fragments[position]) == self.nfrags:
                position = random.choice(sequence_position_range)                
                
            #sample a single fragment at the position and replace torsions in a copied version of the protein
            new_protein = self.perturb_fragment(self.current_protein, position)
            
            #measure energy after replacing the fragment
            new_energy = self.compute_energy(new_protein)
            best_energy = self.compute_energy(self.best_protein)
            
            #update self.best_protein if necessary
            if new_energy < best_energy: 
                self.best_protein = new_protein
            
            #accept or reject move using Metropolis criterion and update T as necessary
            accept_probability = self.metropolis_accept(self.current_protein, new_protein, self.current_T)
            accept = np.random.binomial(size = 1, n = 1, p = accept_probability)[0] 
            
            if accept == 1: 
                self.current_protein = new_protein
                self.current_T = self.anneal_temp(self.current_T)
                self.current_energy = new_energy
                self.current_iteration = self.current_iteration + 1
                break
        

    def simulate(self):
        """
        Run ONE full MCMC simulation from start_temp to end_temp. 
        Be sure to save the best (lowest-energy) structure, so you can access it after.
        It is also a good idea to track certain variables during the simulation (temp, energy, and more).
        -------- 
        Params:
            start_temp = starting temperature for the annealing
            end_temp = ending temperature for the annealing
        Returns:
            log_table = a pandas dataFrame that holds the iteration number, energy, and temperature at each step of sampling
        """
        #loop to perform additional steps until the current temperature is no longer greater than the ending_temperature
        while self.current_T >= self.end_temp: 
            self.step(self.current_T)
            
            #log various parameters that changed in the MCMCSampler object after a single step
            self.temperature.append(self.current_T)
            self.iteration.append(self.current_iteration)
            self.energy.append(self.current_energy)
        #return a pandas dataframe that will hold all of the information requested above
        log_table = pd.DataFrame(list(zip(self.iteration, self.energy, self.temperature)), columns =['iteration', 'energy', 'temperature']) 
        return(log_table)





    



