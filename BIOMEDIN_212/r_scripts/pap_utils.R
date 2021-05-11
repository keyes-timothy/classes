#pap_read_tif
# description: 
#
# input:
#
# output:
pap_read_tif <- function(filepath) {
  
  read_tif <- quietly(tiff::readTIFF)
  
  result <- 
    filepath %>% 
    read_tif(source = .) %>% 
    pluck("result") %>% 
    as_tibble() %>% 
    mutate(y = n():1) %>% 
    pivot_longer(
      cols = -y, 
      names_to = "x", 
      values_to = "values", 
      names_prefix = "V", 
    ) %>% 
    mutate(
      x = as.integer(x), 
      values = if_else(values == 1, "surfactant", "no surfactant") 
    ) %>% 
    select(x, y, values) %>% 
    arrange(x, y)
  
  return(result)
}
