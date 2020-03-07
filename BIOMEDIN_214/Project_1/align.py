
"""

This file provides skeleton code for align.py. 

Locations with "FILL IN" in comments are where you need to add code.

Note - you do not need to follow this set up! It is just a suggestion, and may help for program design and testing.

Usage: python align.py input_file output_file

"""

import sys
import pandas as pd
import numpy as np

#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a, b):
    """
    Checks if two floating point numbers are equivalent.
    """
    epsilon = 10**(-6) 
    return (abs(a - b) < epsilon)
    

#### ------- CLASSES ------- ####

class MatchMatrix:
    """
    Match matrix class stores the scores of matches in a data structure
    """
    def __init__(self, alphabet_a, alphabet_b, len_alphabet_a, len_alphabet_b): #might have to change arguments here
        #make a pandas dataFrame with each element initialized as NaN.
        self.MM = pd.DataFrame(index=range(len_alphabet_a),columns=range(len_alphabet_b))
        self.MM.columns = list(alphabet_b)
        self.MM.index = list(alphabet_a) 
    
    def set_score(self, a, b, score):
        """
        Updates or adds a score for a specified match

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
           score = the score to set it for
        """
        self.MM.loc[[a], [b]] = score
    
    def get_score(self, a, b):
        """
        Returns the score for a particular match, where a is the
        character from sequence a and b is from sequence b.

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
        Returns:
           the score of that match
        """
        return self.MM.at[a, b]

class ScoreMatrix:
    """
    Object to store a score matrix to be generated during the alignment process. The score matrix consists of a pandas DataFrame
    called "score_matrix" that holds the scores (floats) and a pandas DataFrame called "pointer_matrix" that holds 
    a list of tuples formatted as ("matrix_name", row, col) that represent pointers from each entry in score_matrix to their
    predecessors as discussed in class. 
    """
    def __init__(self, name, nrow, ncol):
        self.name = name # identifier for the score matrix - Ix, Iy, or M
        self.nrow = nrow
        self.ncol = ncol
        self.score_matrix = pd.DataFrame(index=range(nrow),columns=range(ncol))
        self.pointer_matrix = pd.DataFrame(index=range(nrow),columns=range(ncol)) 
        
    def get_score(self, row, col):
        """
        Returns the score at a particular entry in the score_matrix. 
        
        Inputs: 
            row = index of the row (0-indexed) of the score being requested. 
            col = index of the column (0-indexed) of the column being requested. 
        
        Outputs:
            returns a float representing the score at entry [row, col]. 
        """
        return self.score_matrix.iloc[row, col]
    
    def set_score(self, row, col, score):    
        """
        Sets the score for a particular entry in the score_matrix. 
        
        Inputs: 
            row = index of the row (0-indexed) of the score being set. 
            col = index of the column (0-indexed) of the column being set. 
            score = value of the score to be placed at entry [row, col]. 
        
        Outputs:
            None
        """
        self.score_matrix.iloc[row, col] = score
        
    def get_pointers(self, row, col):
        """
        Returns the indices of the entries that are pointed to
        This should be formatted as a list of tuples:
         ex. [(1,1), (1,0)]
        """
        return self.pointer_matrix.iloc[row, col]
    
    def set_pointers(self, row, col, pointers):
        """
        Sets the pointers for a particular entry in the pointer_matrix. 
        
        Inputs: 
            row = index of the row (0-indexed) of the score being set. 
            col = index of the column (0-indexed) of the column being set. 
            pointers = tuple formatted as ("matrix_name", row, col) to be placed at entry [row, col]. 
            
        Outputs:
            None
        """
        self.pointer_matrix.iloc[row, col] = pointers
    
    def print_scores(self):
        """
        Returns a nicely formatted string containing the scores in the score matrix. Use this for debugging!
        
        Example:
        M=
            0.0, 0.0, 0.0, 0.0, 0.0
            0.0, 1.0, 0.0, 0.0, 0.0
            0.0, 1.0, 1.0, 1.0, 1.0
            0.0, 0.0, 1.0, 1.0, 1.0
            0.0, 0.0, 2.0, 2.0, 1.0
            0.0, 0.0, 1.0, 2.0, 3.0
            
        """
        return str(self.score_matrix) 
    
    def print_pointers(self):
        """
        Returns a nicely formatted string containing the pointers for each entry in the score matrix. Use this for debugging!
        """
        return str(self.pointer_matrix)

class AlignmentParameters:
    """
    Object to hold a set of alignment parameters from an input file.
    """
    
    def __init__(self):
        # default values for variables that are filled in by reading
        # the input alignment file
        self.seq_a = ""
        self.seq_b = ""
        self.global_alignment = False 
        self.dx = 0
        self.ex = 0
        self.dy = 0
        self.ey = 0
        self.alphabet_a = "" 
        self.alphabet_b = ""
        self.len_alphabet_a = 0
        self.len_alphabet_b = 0
        self.match_matrix = MatchMatrix(self.alphabet_a, self.alphabet_b, self.len_alphabet_a, self.len_alphabet_b)
    
    def load_params_from_file(self, input_file): 
        """
        Reads the parameters from an input file and stores in the object

        Input:
           input_file = specially formatted alignment input file
        """
        #open input file
        with open(input_file, 'r') as my_file:
            self.seq_a = my_file.readline().strip()
            self.seq_b = my_file.readline().strip()
            
            #extract alignment type information
            if my_file.readline().strip() == '0': 
                self.global_alignment = True
                
            #extract information about gap penalties
            self.dx, self.ex, self.dy, self.ey = [float(i) for i in my_file.readline().strip().split()]
            
            #extract information about alphabets
            self.len_alphabet_a = int(my_file.readline().strip())
            self.alphabet_a = my_file.readline().strip()
            self.len_alphabet_b = int(my_file.readline().strip())
            self.alphabet_b = my_file.readline().strip()
            my_file.close()
        
        #extract information about the match matrix    
        self.match_matrix = MatchMatrix(self.alphabet_a, self.alphabet_b, self.len_alphabet_a, self.len_alphabet_b)
        mm_entries = pd.read_csv(input_file, skiprows = [i for i in range(8)], header = None, delim_whitespace=True)
        for i in range(len(mm_entries.index)): 
            my_a = mm_entries.iloc[i, 2]
            my_b = mm_entries.iloc[i, 3]
            my_score = mm_entries.iloc[i, 4]
            self.match_matrix.set_score(my_a, my_b, my_score)


class Align:
    """
    Object to hold and run an alignment; running is accomplished by using "align()"
    """
    def __init__(self, input_file, output_file):
        """
        Input:
            input_file = file with the input for running an alignment
            output_file = file to write the output alignments to
        """
        self.input_file = input_file
        self.output_file = output_file
        self.align_params = AlignmentParameters() 
        
        ### FILL IN - note: be careful about how you initialize these! ### Might still have to change them...
        self.m_matrix = ScoreMatrix("M", 0, 0)
        self.ix_matrix = ScoreMatrix("Ix", 0, 0)
        self.iy_matrix = ScoreMatrix("Iy", 0, 0)
        
    def align(self):
        """
        Main method for running alignment.
        """
        
        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)
        
        # populate the score matrices based on the input parameters
        self.populate_score_matrices() 
        
        # perform a traceback and write the output to an output file
        max_val, max_loc = self.find_traceback_start()
        
        self.traceback(max_loc)
        
        self.write_output(max_val)
        
    def populate_score_matrices(self):
        """
        Method to populate the score matrices based on the data in align_params.
        Should call update(i,j) for each entry in the score matrices
        """
        
        #re-initialize the score matrices
        self.m_matrix = ScoreMatrix('M', len(self.align_params.seq_a) + 1, len(self.align_params.seq_b) + 1)
        self.m_matrix.score_matrix.iloc[0,:] = 0
        self.m_matrix.score_matrix.iloc[:,0] = 0
        self.m_matrix.pointer_matrix.iloc[0,:] = None
        self.m_matrix.pointer_matrix.iloc[:,0] = None
        
        self.ix_matrix = ScoreMatrix('M', len(self.align_params.seq_a) + 1, len(self.align_params.seq_b) + 1)
        self.ix_matrix.score_matrix.iloc[0,:] = 0
        self.ix_matrix.score_matrix.iloc[:,0] = 0
        self.ix_matrix.pointer_matrix.iloc[0,:] = None
        self.ix_matrix.pointer_matrix.iloc[:,0] = None
        
        self.iy_matrix = ScoreMatrix('M', len(self.align_params.seq_a) + 1, len(self.align_params.seq_b) + 1)
        self.iy_matrix.score_matrix.iloc[0,:] = 0
        self.iy_matrix.score_matrix.iloc[:,0] = 0
        self.iy_matrix.pointer_matrix.iloc[0,:] = None
        self.iy_matrix.pointer_matrix.iloc[:,0] = None
        
        #call update(i, j) on each entry in the score matrices 
        for i in range(1, len(self.align_params.seq_a)+1): 
            for j in range(1, len(self.align_params.seq_b)+1): 
                self.update(i, j) 
        
    def update(self, row, col):
        """
        Method to update the matrices at a given row and column index.
        Input:
           row = the row index to update
           col = the column index to update
        """
        self.update_m(row, col)
        self.update_ix(row, col)
        self.update_iy(row, col)
    
    def update_m(self, row, col):
        """
        Helper method to update the m_matrix at a given row and column index using the recursion equations discussed in class. 
            Input: 
                row = the row index to update
                col = the column index to update
        """
        #extract score from current match from match matrix
        match_value = self.align_params.match_matrix.get_score(a = self.align_params.seq_a[row-1], b = self.align_params.seq_b[col-1])
        
        #first recursion
        current_max_score = self.m_matrix.get_score(row-1, col-1) + match_value
        current_pointers = [('M', row-1, col-1)]
        
        #second recursion
        new_value = self.ix_matrix.get_score(row-1, col-1) + match_value
        if fuzzy_equals(new_value, current_max_score): 
            current_pointers.append(('Ix', row-1, col-1))
        elif new_value > current_max_score: 
            current_max_score = new_value
            current_pointers = [('Ix', row-1, col-1)]
            
        #third recursion
        new_value = self.iy_matrix.get_score(row-1, col-1) + match_value 
        if fuzzy_equals(new_value, current_max_score):
            current_pointers.append(('Iy', row-1, col-1))
        elif new_value > current_max_score:
            current_max_score = new_value
            current_pointers = [('Iy', row-1, col-1)]
        
        #so that local alignment will not have negative values (will reset local alignment)
        if not self.align_params.global_alignment: 
            if 0 > current_max_score: 
                current_max_score = 0
                current_pointers = None
        
        #update m_matrix value and pointers accordingly
        self.m_matrix.set_score(row, col, current_max_score)
        self.m_matrix.set_pointers(row, col, current_pointers)
    
    def update_ix(self, row, col):
        """
        Helper method to update the ix_matrix at a given row and column index using the recursion equations discussed in class. 
            Input: 
                row = the row index to update
                col = the column index to update
        """
        current_max_score = self.m_matrix.get_score(row-1, col) - self.align_params.dy
        current_pointers = [('M', row-1, col)]
        
        new_value = self.ix_matrix.get_score(row-1, col) - self.align_params.ey
        if fuzzy_equals(new_value, current_max_score): 
            current_pointers.append(('Ix', row-1, col))
        elif new_value > current_max_score: 
            current_max_score = new_value
            current_pointers = [('Ix', row-1, col)]
            
        if not self.align_params.global_alignment: 
            if 0 > current_max_score: 
                current_max_score = 0
                current_pointers = None
        
        #update ix_matrix value and pointers accordingly
        self.ix_matrix.set_score(row = row, col = col, score = current_max_score)
        self.ix_matrix.set_pointers(row = row, col = col, pointers = current_pointers)
        
    def update_iy(self, row, col):
        """
        Helper method to update the iy_matrix at a given row and column index using the recursion equations discussed in class. 
            Input: 
                row = the row index to update
                col = the column index to update
        """
        #first recursion
        current_max_score = self.m_matrix.get_score(row, col-1) - self.align_params.dx 
        current_pointers = [('M', row, col-1)]
        
        #second recursion
        new_value = self.iy_matrix.get_score(row, col-1) - self.align_params.ex
        if fuzzy_equals(new_value, current_max_score): 
            current_pointers.append(('Iy', row, col-1))
        elif new_value > current_max_score: 
            current_max_score = new_value
            current_pointers = [('Iy', row, col-1)]
        
        #take care of local alignment case
        if not self.align_params.global_alignment: 
            if 0 > current_max_score: 
                current_max_score = 0
                current_pointers = None
        
        #update iy_matrix value and pointers accordingly
        self.iy_matrix.set_score(row = row, col = col, score = current_max_score)
        self.iy_matrix.set_pointers(row = row, col = col, pointers = current_pointers)
        
    def find_traceback_start(self):
        """
        Finds the location to start the traceback. For global, it finds the maximum value in the rightmost column and the 
        lowermost row in any of the ScoreMatrices. For local, it finds the maximum value at any entry in the ScoreMatrices. 
        
        Returns:
            (max_val, max_loc) where max_val is the best score
            max_loc is a list [] containing tuples with the (i,j) location(s) to start the traceback
             (ex. [(1,2), (3,4)])
        """
        if self.align_params.global_alignment: 
            #look for max value in any of the last rows or columns of our 3 matrices
            m_col_max = self.m_matrix.score_matrix.max().iloc[len(self.align_params.seq_b)]
            m_row_max = self.m_matrix.score_matrix.max(axis = 1).iloc[len(self.align_params.seq_a)]
            ix_col_max = self.ix_matrix.score_matrix.max().iloc[len(self.align_params.seq_b)]
            ix_row_max = self.ix_matrix.score_matrix.max(axis = 1).iloc[len(self.align_params.seq_a)]
            iy_col_max = self.iy_matrix.score_matrix.max().iloc[len(self.align_params.seq_b)]
            iy_row_max = self.iy_matrix.score_matrix.max(axis = 1).iloc[len(self.align_params.seq_a)]
            
            max_val = max([m_col_max, m_row_max, ix_col_max, ix_row_max, iy_col_max, iy_row_max])
            
            #save all the pointers for positions that match the max value
            max_loc = []
            final_row = len(self.m_matrix.pointer_matrix.index)-1
            final_col = len(self.m_matrix.pointer_matrix.columns)-1
            
            for col in self.m_matrix.pointer_matrix.columns: 
                #iterate across all columns of final row
                if fuzzy_equals(self.m_matrix.get_score(final_row, col), max_val): 
                    max_loc.append(("M", final_row, col))
                if fuzzy_equals(self.ix_matrix.get_score(final_row, col), max_val): 
                    max_loc.append(("Ix", final_row, col))
                if fuzzy_equals(self.iy_matrix.get_score(final_row, col), max_val): 
                    max_loc.append(("Iy", final_row, col))
                
            
            for row in self.m_matrix.pointer_matrix.index: 
                #iterate across all rows of final column
                if fuzzy_equals(self.m_matrix.get_score(row, final_col), max_val):
                    max_loc.append(("M", row, final_col))
                if fuzzy_equals(self.ix_matrix.get_score(row, final_col), max_val): 
                    max_loc.append(("Ix", row, final_col))
                if fuzzy_equals(self.iy_matrix.get_score(row, final_col), max_val): 
                    max_loc.append(("Iy", row, final_col))
                
            #to remove any duplicates (the bottom right entry is checked twice)
            max_loc = list(set(max_loc))
            
            return (max_val, max_loc)
            
        else: 
            #look for max value anywhere in the 3 matrices for local alignment
            m_max = self.m_matrix.score_matrix.values.max()
            ix_max = self.ix_matrix.score_matrix.values.max()
            iy_max = self.iy_matrix.score_matrix.values.max()
            max_val = max([m_max, ix_max, iy_max])
            
            max_loc = []
            for row in self.m_matrix.pointer_matrix.index: 
                for col in self.m_matrix.pointer_matrix.columns: 
                    if fuzzy_equals(self.m_matrix.get_score(row, col), max_val): 
                        max_loc.append(("M", row, col))
                    if fuzzy_equals(self.ix_matrix.get_score(row, col), max_val): 
                        max_loc.append(("Ix", row, col))
                    if fuzzy_equals(self.iy_matrix.get_score(row, col), max_val): 
                        max_loc.append(("Iy", row, col))
            
            return (max_val, max_loc)
    
    def end_traceback(self, my_row, my_col, my_score): 
        """
        A helper function that returns a boolean value to indicate if it is time to end the traceback procedure at a 
        particular entry in the ScoreMatrices. 
        
        Inputs: 
            my_row = row (0-indexed) of the ScoreMatrix entry being considered
            my_col = column (0-indexed) of the ScoreMatrix entry being considered
            my_score = score of the ScoreMatrix entry being considered 
        """
        #for global alignment, end ithe traceback f either sequence has been fully traced back to its first letter
        if self.align_params.global_alignment: 
            if my_row == 1 or my_col == 1: 
                return True
            else: 
                return False
       #for local alignment, end the traceback if any scores are equal to 0
        else: 
            if my_score == 0: 
                return True
            else: 
                return False
    
    
    def find_paths(self, start_loc):
        """
        A function that returns a list of lists of tuples denoting the order in which entries were encountered during the traceback 
        procedure. Each list represents a distinct path, and each tuple represents a ScoreMatrix entry within that path. This function
        is similar to a depth-first search.
        
        These pathways are used to assemble the sequences that were aligned in  the self.traceback() function.
        
        Inputs: 
            start_loc = a tuple representing the starting location of the traceback.
        
        Returns: 
            A list of lists of tuples (pointers) denoting the order in which entries were encountered during the traceback. 
        """
        to_visit = [(start_loc, [start_loc])]
        while len(to_visit) > 0: 
            current_loc = to_visit.pop()
            current_pointer, current_path = current_loc
            current_matrix, current_row, current_col = current_pointer
            if (current_matrix == 'M'): 
                current_score = self.m_matrix.get_score(current_row, current_col)
                next_pointers = self.m_matrix.get_pointers(current_row, current_col)
            elif (current_matrix == "Ix"): 
                current_score = self.ix_matrix.get_score(current_row, current_col)
                next_pointers = self.ix_matrix.get_pointers(current_row, current_col)
            else: 
                current_score = self.iy_matrix.get_score(current_row, current_col)
                next_pointers = self.iy_matrix.get_pointers(current_row, current_col)
                
            if(self.end_traceback(current_row, current_col, current_score)): 
                if(self.align_params.global_alignment == True): 
                    yield (current_path)
                else: 
                    if (current_matrix == 'M'): #prevent reading out 0s from Ix and Iy encountered during local alignment
                        yield (current_path[:-1]) #have to cut the last pointer out because, 
                                                  #if we're getting a 0, then this pair isn't included in the answer. 
            else: 
                for new_pointer in next_pointers:
                    to_visit.append((new_pointer, current_path + [new_pointer]))
    
    def find_new_letters(self, my_matrix, my_row, my_col): 
        """
        A helper function that returns the alignment pair associated with a given entry of a given ScoreMatrix. 
        
        Inputs: 
            my_matrix = a string ("M", "Ix", or "Iy") that indicates the name of the ScoreMatrix. 
            my_row = an integer indicating the (0-indexed) row of the ScoreMatrix
            my_col = an integer indicating the (0-indexed) column of the ScoreMatrix
            
        Returns: 
            two strings called "my_a" and "my_b" that represent the character to input into the alignment in sequence A 
            and sequence B, respectively.
        """
        if(my_matrix == "M"):
                my_letter_a = str('_' + self.align_params.seq_a)[my_row]
                my_letter_b = str('_' + self.align_params.seq_b)[my_col]
        elif(my_matrix == "Ix"): 
            my_letter_a = str('_' + self.align_params.seq_a)[my_row]
            my_letter_b = '_'
        else:
            my_letter_a =  '_'
            my_letter_b = str('_' + self.align_params.seq_b)[my_col]
        return my_letter_a, my_letter_b
                
    
    def traceback(self, max_loc):
        """
        Performs the traceback procedure by calling find_paths on each pointer in max_loc, then converting each of those paths 
        into an alignment string. Onces the alignment strings for Sequence A and Sequence B are computed, they are stored in the 
        class variables self.aligned_a and self.aligned_b as lists of strings (corresponding items in self.aligned_a and 
        self.aligned_b are paired with one another). 
        
        Inputs: 
            max_loc =  a list [] containing tuple pointers with the (row,col) location(s) to start the traceback
        """
        my_paths = []
        for loc in max_loc: 
            my_paths.extend(list(self.find_paths(loc)))
        
        #follow each path and return actual alignments for a and b
        self.aligned_a = []
        self.aligned_b = []
        for path in my_paths: #for each path in the list of paths
            my_a = ""
            my_b = ""
            for location in path: #for each tuple in a given path
                my_matrix, my_row, my_col = location
                new_a, new_b = self.find_new_letters(my_matrix, my_row, my_col)
                my_a = new_a + my_a 
                my_b = new_b + my_b
            self.aligned_a.append(my_a)
            self.aligned_b.append(my_b)
    
    
    #write output
    def write_output(self, max_val):
        """
        A function that writes the output file as described on the project page. 
        
        Inputs: 
            max_val = a float representing the alignments' maximum score. 
            
        Outputs: 
            my_output = a .output file formatted as requested. 
        """
        my_output = open(self.output_file,'w+')
        my_output.write(str(round(max_val, 1)))
        my_output.write('\n')
        
        for i in range(len(self.aligned_a)): 
            my_output.write('\n')
            my_output.write(self.aligned_a[i] + '\n')
            my_output.write(self.aligned_b[i] + '\n')
        
        my_output.close()
        

def main():

    # check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return
        
    # input variables
    input_file = sys.argv[1]
    print(input_file) #for testing
    output_file = sys.argv[2]
    print(output_file) #for testing
    
    # create an align object and run
    align = Align(input_file, output_file)
    align.align()


if __name__=="__main__":
    main()
