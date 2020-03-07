#my_unit_tests.py

#set up files 
input_file = 'alignment_example6.input'
output_file = 'test_results.output'

#testing

my_Align_test = Align(input_file, output_file)

self = my_Align_test
self.align_params.load_params_from_file(self.input_file)

print(self.align_params.seq_a)
print(self.align_params.seq_b)
print(self.align_params.match_matrix.MM)


self.m_matrix.print_scores()
self.m_matrix.print_pointers()
self.ix_matrix.print_scores()
self.ix_matrix.print_pointers()
self.iy_matrix.print_scores()
self.iy_matrix.print_pointers()

#test dynamic programming of the matrices
self.populate_score_matrices()

self.m_matrix.print_scores()
self.m_matrix.print_pointers()
self.ix_matrix.print_scores()
self.ix_matrix.print_pointers()
self.iy_matrix.print_scores()
self.iy_matrix.print_pointers()

#test finding the starting location of the traceback 
max_val, max_loc = self.find_traceback_start()

my_paths = list(self.find_paths(max_loc[0]))

self.traceback(max_loc)

self.write_output(max_val)

#test alignment function

start_point = ('M', 2, 3)




