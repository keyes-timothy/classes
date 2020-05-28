#sample_clustered_data.R

# Description

# Author: Timothy Keyes
# Version: 2020-06-20

# Libraries
library(flowCore)
library(tidyverse)

# Parameters

healthy_path <- here::here("CS_230", "course_project", "data", "healthy.rds")
cancer_path <- here::here("CS_230", "course_project", "data", "population_data.rds")
out_path <- here::here("CS_230", "course_project", "data")
marker_path <- here::here("CS_230", "course_project", "docs", "ALL_panel.csv")

CLASSIFIER_POPULATIONS <-
  c(
    "HSC", 
    "Progenitor_1", 
    "Progenitor_2", 
    "Progenitor_3", 
    "Pre_Pro_B", 
    "Pro_B1", 
    "Pro_B2", 
    "Pre_B1", 
    "Pre_B2", 
    "Early_Progenitors", 
    "Late_Progenitors", 
    "Immature_B1", 
    "Immature_B2", 
    "Mature_B", 
    "Mature_Non_B"
  )

#===============================================================================
source(
  file = 
    file.path("~", "GitHub", "classes", "CS_230", "course_project", "r_scripts", "build_classifier.R"
    )
)

source(
  file = 
    file.path("~", "GitHub", "classes", "CS_230", "course_project", "r_scripts", "apply_classifier.R"
    )
)
#===============================================================================

my_classifier <- 
  build_classifier(pops_path = )





my_data <- 
  input_path %>% 
  read_rds() %>% 
  dplyr::filter(stimulation == "Basal", !str_detect(patient, "Relapse")) %>% 
  group_by(patient) %>% 
  dplyr::filter(n() >= 10000) %>% 
  ungroup() %>% 
  sample_n(size = 10000)

write_rds(x = my_data, path = file.path(out_data, "sampled_cluster_data.rds"))
