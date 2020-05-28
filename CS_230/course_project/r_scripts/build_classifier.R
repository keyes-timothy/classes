################      function for building classifier       #################

build_classifier <- function(
  pops_path = NULL, 
  pop_names = NULL, 
  classifier_names = NULL, 
  other_marker_names = NULL, 
  metadata = NULL
) { 
  message("Building classifier from healthy cells...")
  # A vector for saving the distance from the most distant cell in each reference population
  # Not fully implemented yet
  farthest <- c()
  
  #reading in population data
  
  message("Reading in healthy populations...")
  healthy_populations <- 
    read_rds(pops_path) %>% 
    ungroup() %>% 
    select(one_of(metadata, classifier_names, other_marker_names)) %>% 
    mutate_at(vars(-one_of(metadata)), ~ asinh(. / 5)) %>% 
    nest(
      classifier_markers = one_of(classifier_names), 
      remaining_markers = one_of(other_marker_names)
    )
  
  # Calculate mean and covariance matrix for each population
  message("Calculating mean and covariance matrix for all healthy populations")
  healthy_populations <- 
    healthy_populations %>% 
    mutate(
      population, 
      means = 
        map(
          classifier_markers, 
          ~ (.) %>% 
            summarize_all(mean) %>% 
            pivot_longer(everything()) %>% 
            deframe()
        ), 
      cov_mat = map(classifier_markers, get_cov)
    )
  message("Done! Returning classifier_fit object")
  healthy_populations
}
