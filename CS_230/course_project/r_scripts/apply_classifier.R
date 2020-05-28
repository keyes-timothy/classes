##############      function for applying classifier       #################


apply_classifier <- function(
  classifier_fit = NULL, 
  nCores = NULL, 
  cancer_path = NULL, 
  metadata = "file_name"
) {
  
  classifier_names <- colnames(classifier_fit$cov_mat[[1]])
  
  #read in cancer data
  cancer_data <- 
    read_rds(cancer_path) %>% 
    ungroup() %>% #should take this out eventually
    mutate_if(is.numeric, ~ asinh(. / 5)) %>% 
    nest(
      classifier_markers = one_of(classifier_names), 
      remaining_markers = -one_of(classifier_names, metadata)
    )
  
  my_cluster <- makeCluster(nCores, outfile = "")
  registerDoParallel(my_cluster)
  
  cancer_data <- 
    cancer_data %>% 
    mutate(
      classification_data = 
        foreach(
          expression_matrix = cancer_data$classifier_markers,
          .combine = list, 
          .packages = c("dplyr", "purrr"), 
          .export = c("classify_cell", "classifier_fit"), 
          .multicombine = TRUE, 
          .maxcombine = nrow(cancer_data)
        ) %dopar%
        classify_cell(classifier_fit, expression_matrix)
    )
  stopCluster(my_cluster)
  
  cancer_data
}
