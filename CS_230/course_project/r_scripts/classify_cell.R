# Function to assign cells to populations
classify_cell <- 
  function(classifier_fit, expression_matrix){
    # Calculate Mahalanobis distance to each population
    classifications <- 
      classifier_fit %>%
      mutate(
        classification_data = 
          map2(
            .x = means,
            .y = cov_mat, 
            .f = ~ mahalanobis(expression_matrix, center = .x, cov = .y)
          )
      ) %>% 
      pull(classification_data) %>% 
      bind_cols()
    
    colnames(classifications) <- classifier_fit$population
    
    
    #This has to be optimized
    classifications <- 
      classifications %>% 
      mutate(
        row_min = pmap_dbl(classifications, min)
      )
    
    #this part has to be optimized
    classifications <- 
      classifications %>% 
      mutate(
        row_index = 
          map_int(
            1:nrow(classifications), 
            ~ which(classifications[.,] == classifications$row_min[[.]])[[1]]),
        classified_population = 
          map_chr(
            row_index,
            ~ colnames(classifications)[[.]]
          )
      ) %>% 
      dplyr::select(-row_index, -row_min)
    
    classifications

  }

