#sample_clustered_data.R

# Description

# Author: Timothy Keyes
# Version: 2020-06-20

# Libraries
library(flowCore)
library(tidyverse)
library(doParallel)

#source functions
source(here::here("CS_230", "course_project", "r_scripts", "get_cov.R"))
source(here::here("CS_230", "course_project", "r_scripts", "apply_classifier.R"))
source(here::here("CS_230", "course_project", "r_scripts", "classify_cell.R"))


# Parameters

healthy_path <- here::here("CS_230", "course_project", "data", "population_data.rds")
cancer_path <- here::here("CS_230", "course_project", "data", "sampled_DDPR_data.rds")
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

CLASSIFIER_MARKERS <- 
  c(
    'CD19', 
    'CD20', 
    'CD34', 
    'CD38', 
    'Tdt',
    'CD45', 
    'IgMs', 
    'IgMi',
    'CD179a', 
    'CD179b', 
    'CD127'
  )

MARKERS <- 
  c(
    'CD19', 
    'CD20', 
    'CD34', 
    'CD38', 
    'IgMi', 
    'IgMs', 
    'CD179a', 
    'CD179b',
    'CD127', 
    'Tdt',
    'CD45', 
    'PLCg2', 
    'CD22', 
    'p4EBP1', 
    'Ikaros', 
    'CD79b', 
    'pSTAT5',
    'CD123',
    'Ki67', 
    'IgL kappa',
    'IKAROS_i', 
    'CD10', 
    'pAkt',
    'CD24', 
    'CRLF2', 
    'RAG1', 
    'Pax5',
    'pSyk',
    'CD43',
    'CD58',
    'HIT3a',
    'CD16',
    'pS6',
    'pErk',
    'HLADR',
    'pCreb'  
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
  build_classifier(
    pops_path = healthy_path,
    pop_names = CLASSIFIER_POPULATIONS,
    classifier_names = CLASSIFIER_MARKERS,
    other_marker_names = setdiff(MARKERS, CLASSIFIER_MARKERS), 
    metadata = c("population")
  )

my_result <- 
  apply_classifier(
    classifier_fit = my_classifier,
    nCores = 10, 
    cancer_path = cancer_path, 
    metadata = "file_name"
  ) %>%
  unnest(cols = c(classifier_markers, remaining_markers, classification_data))

write_rds(x = my_result, path = file.path(out_path, "classifier_saucie_result.rds"))
