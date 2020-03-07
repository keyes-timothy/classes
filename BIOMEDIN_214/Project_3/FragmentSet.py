import utils
import numpy as np
import pandas as pd
import os

class FragmentSet(object):
    def __init__(self, fragfile, rmsdfile):
        """
        This class contains the fragment library for the input protein. It must do the following:
        - Read in fragment file and parse fragments at each position. Fragment files are of the form <protein>_<frag_len>mers.frag
        - Read in RMSD file containing pre-calculated RMSD to native structure for each fragment at each position.
        - Based on fragments and their corresponding RMSDs, rank fragments at each position by RMSD
        """
        #ultimately gives a a fragment library in which each position in the library (a list) contains a dictionary of each fragment 
        #to a list of tuples encoding the phi and psi angles of its 1st, 2nd...nth position. 
        my_fragment_library = []
        current_position_frags = {}
        with open(fragfile) as my_file:
            while True: 
                my_line = my_file.readline()
                if my_line == '': 
                    my_file.close()
                    break
                else: 
                    while True: #loop through all lines until you find a "position:" - single position loop (must be put into a larger loop)
                        my_line = my_line.strip().split()
                        if len(my_line) == 0: 
                            my_fragment_library.append(current_position_frags)
                            current_position_frags = {}
                            break
                        elif my_line[0] == "position:": 
                            #store old position's dictionary
                            my_fragment_library.append(current_position_frags)
                            current_position_frags = {}
                            junk = my_file.readline()
                            break
                        else: 
                            #set up new fragment
                            current_frag_name = my_line[len(my_line)-1]
                            if current_frag_name[0] == 'F': 
                                current_frag_name = current_frag_name[1:len(current_frag_name)]
                            current_frag_name = int(current_frag_name)
                            current_phi = float(my_line[5])
                            current_psi = float(my_line[6])
                            current_frag_torsions = [(current_phi, current_psi)]
                            
                            #loop through the rest of the lines containing information about the current fragment
                            while True: 
                                my_line = my_file.readline().strip().split()
                                if my_line == []:
                                    #once the fragment is done, add it to the dictionary containing all fragments for that section
                                    current_position_frags[current_frag_name] = current_frag_torsions
                                    my_line = my_file.readline()
                                    break
                                else: 
                                    current_phi = float(my_line[5])
                                    current_psi = float(my_line[6])
                                    current_frag_torsions.append((current_phi, current_psi))
                                    
            self.my_fragment_library = my_fragment_library
            
            #read in RMSD file (using pandas) 
            RMSDs = pd.read_csv(rmsdfile, sep = '\t', header = None)
            RMSDs.columns = ['position', 'fragment', 'RMSD']
            self.RMSDs = RMSDs.sort_values('RMSD', ascending = True).groupby('position')


    def get_lowRMS_fragments(self, pos, N):
        """
        Returns the top-ranked fragments by RMSD at a defined position in the chain
        --------
        Params
            - pos (int): fragment position in chain (1-indexed)
			- N (int): number of fragments to return
		Returns
			- lowRMS_fragments (list): top N fragments at pos by RMSD. This should be a list of lists of (phi, psi) tuples. 
			  For example, a 3-mer fragment could be represented as the following: [(-60.892, 142.456), (-72.281, 128.933), (-132.337, -175.477)]
		"""
        low_fragment_ids = self.RMSDs.head(N).loc[self.RMSDs.head(N)['position'] == pos].fragment.values
        
        lowRMS_fragments = []
        for low_fragment in low_fragment_ids: 
            lowRMS_fragments.append(self.my_fragment_library[pos][low_fragment])
        return(lowRMS_fragments)





