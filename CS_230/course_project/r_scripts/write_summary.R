# Making summary tables (excel) and printing information about run to a text file
write_summary <- function(classifications, pop_names, output_folder, start_time, stop_time, pops_dir, 
                          sample_folder, threshold, classifier_names, give_summary) {
  
  # Making summary tables
  table1 <- NULL; table2 <- NULL
  for (i in 1:length(classifications)) {
    classes <- table(factor(classifications[[i]]$MahID, levels = 0:length(pop_names)))
    table1 <- rbind(table1, classes)
    table2 <- rbind(table2, classes/sum(classes)*100)
  }
  rownames(table1) <- rownames(table2) <- names(classifications); colnames(table1) <- colnames(table2) <- c('Unclassified', pop_names)
  if(give_summary == TRUE){
    write.xlsx(x = list('Counts' = table1, 'Percent' = table2), file = file.path(output_folder, 'Classification_Summary.xlsx'), row.names = T)
    message('Finished writing excel output file.')
  }
  my.table <- table1[,1:length(pop_names)+1]
  
  # Writing info about run to a file
  write(paste0('The Developmental Classifier started parsing the input cells at ', start_time, ', and finished at ', stop_time, '.\n\n',
               'A total of ', sum(table1), ' cells from ', length(classifications), ' files were parsed and ', round(sum(table1[,1])/sum(table1)*100, 2), ' % of these were left unclassified.\n\n',
               'The reference file directory was ', pops_dir, ' and the number of reference populations was ', length(pop_names), '.\n\n',
               'The input directory was ', sample_folder, ' and the output directory was ', output_folder, '.\n\n',
               'The threshold(s) for unclassified cells was ', threshold, ' and the number of markers used was ', length(classifier_names), '.\n\n',
               'List of reference files:\n', paste(list.files(pops_dir, pattern = paste0('^(', paste(sprintf('%02d', 1:length(pop_names)), collapse = '|'), ').+.fcs*')), collapse = ', '), '\n\n',
               'List of input files:\n', paste(names(classifications), collapse = ', '), '\n\n',
               'List of applied markers:\n', paste(classifier_names, collapse = ', ')), file = file.path(output_folder, 'run_info.txt'))
  return(my.table)
}