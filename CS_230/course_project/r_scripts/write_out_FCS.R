# Making a single FCS file per input sample containing the raw data and columns containing information about classification and make one file per population
write_out_FCS <- function(file_to_classify, fcs, class_stats, output_folder, pop_names) {
  
  # Making main FCS with all events
  fcs_edited <- cbind(fcs@exprs, class_stats)
  
  output_file <- paste0(output_folder, '/Classified_', basename(file_to_classify))
  if (file.exists(output_file)) {
    file.remove(output_file)
  }
  write_modified_FCS(as.matrix(fcs_edited), output_file, channel_descriptions = fcs@parameters$desc,
                     reference_description = flowCore::description(fcs))
  
  cat('Main output FCS files written for', basename(file_to_classify), '\n')
  
  
  # Making population FCSs
  population_folder <- paste0(output_folder, '/', file_path_sans_ext(basename(file_to_classify)), '_pops')
  if (file.exists(population_folder)) {
    unlink(population_folder, recursive = T)
  }
  dir.create(population_folder, showWarnings = F)
  
  pop_names <- c('Unclassified', pop_names)
  
  for (pop in unique(class_stats$MahID)) {
    output_file <- paste0(population_folder, '/', sprintf('%02d', pop), '_', gsub(' ', '_', pop_names[pop+1]), '_', basename(file_to_classify))
    
    fcs_pop <- fcs_edited[fcs_edited[,'MahID']==pop,]
    write_modified_FCS(as.matrix(fcs_pop), output_file, channel_descriptions = fcs@parameters$desc,
                       reference_description = flowCore::description(fcs))
  }
  cat('Wrote', length(unique(class_stats$MahID)), 'population FCS files for', basename(file_to_classify), '\n')
}