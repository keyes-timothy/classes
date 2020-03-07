"""
This is the master file, which you should use to set up and run the simulations.
You may define functions or classes as necessary

For an input sequence, do the following (see project page for details):
	1. load sequence from fasta file and initialize protein into extended configuration
	2. Run a series of simulations (n=5 or 10):
		- Perform MCMC sampling with 9-mer fragments from kT=100 to kT=1 (assembly stage)
		- Perform MCMC sampling with 3-mer fragments from kT=1 to kT=0.1 (refinement stage)
		- Take best (lowest-energy) structure after refinement and perform energy minimization (see utils.relax)
		- Log energy and RMSD to native structure after minimization
	3. Visualize lowest-RMSD structure in PyMol

"""

from pyrosetta import *
from pyrosetta.rosetta import *
init(extra_options='-mute all -constant_seed')

from Protein import Protein
from FragmentSampler import MCMCSampler
from FragmentSet import FragmentSet

import utils
import argparse
import numpy as np
import time
import os
import pandas as pd



def main():
    #Read in arguments from the command line
    parser = argparse.ArgumentParser(description = 'Performs ab initio protein structure prediction')
    parser.add_argument('--fasta', help = 'name of the fasta file containing the sequence')
    parser.add_argument('--logdir', help = 'directory to save all log files', default = '.') #changed so does not have the / afterwards
    parser.add_argument('--nsims', type = int, help = 'number of simulations', default = 1)
    parser.add_argument('--nfrags', type = int, help = 'number of fragments to sample at each iteration', default = 3)
    parser.add_argument('--anneal_rate', type = float, help = 'temperature annealing parameter', default = 0.999)
    
    args = parser.parse_args()
        
    #load sequence from fasta file
    my_file = open(args.fasta)
    my_sequence = my_file.readlines()
    my_sequence = my_sequence[1].strip()
    my_file.close()
    
    #initialize protein into extended conformation
    my_protein = Protein(sequence = my_sequence)
    
    #set up FragmentSets for later use 
    protein_name = args.fasta[:-6]
    fragfile_9 = protein_name + '_9mers.frag'
    rmsdfile_9 = protein_name + '_9mers.rmsd'
    fragment_set_9mer = FragmentSet(fragfile_9, rmsdfile_9)
    
    fragfile_3 = protein_name + '_3mers.frag'
    rmsdfile_3 =  protein_name + '_3mers.rmsd'
    fragment_set_3mer = FragmentSet(fragfile_3, rmsdfile_3)
    
    
    #create lists for saving certain information about each simulation 
    sim_numbers = []
    sim_rmsds = []
    sim_energies = []
    sim_proteins = []
    
    #loop through 10 sims
    for sim_number in range(1,args.nsims+1): 
        
        #create FragmentSampler object and run coarse simulation with 9mers
        coarse_sampler = MCMCSampler(my_protein, fragment_set_9mer, nmers = 9, start_temp = 100, end_temp = 1, nfrags = 3, anneal_rate = args.anneal_rate)
        coarse_result = coarse_sampler.simulate()
        
        #create FragmentSampler object and run refinement simulation with 3mers
        fine_sampler = MCMCSampler(coarse_sampler.best_protein, fragment_set_3mer, nmers = 3, start_temp = 1, end_temp = 0.1, nfrags = 3, anneal_rate = args.anneal_rate)
        fine_result = fine_sampler.simulate()
        
        #save the final structure from the simulation
        fine_sampler.best_protein.save_pdb(filename = 'best.pdb')
        
        #save a log_file from the simulation
        total_result = pd.concat([coarse_result, fine_result]) #might need some fixing depending on what the quiz asks for 
        total_result.to_csv((args.logdir + '/sim_' + str(sim_number) + '_log.txt'), sep = '\t')
        
        #relax the simulated structure
        relaxed_protein, final_rmsd, final_energy = utils.relax(pdb = 'best.pdb', native = args.fasta.split(sep = '.')[0] + '.pdb')

        #log some parameters...sim_number, energy, rmsd
        sim_numbers.append(sim_number)
        sim_rmsds.append(final_rmsd)
        sim_energies.append(final_energy)
        sim_proteins.append(relaxed_protein)
        
    #return sim log
    sim_log = pd.DataFrame(list(zip(sim_numbers, sim_energies, sim_rmsds)), columns =['sim_number', 'energy', 'rmsd']) 
    sim_log.to_csv((args.logdir + '/simulation_summary.txt'), sep = '\t', index = False) #may need to be tweaked, i.e. removing indices in first column
    


if __name__=='__main__':
    main()

