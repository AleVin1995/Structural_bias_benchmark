library(tidyverse)
library(magrittr)
library(ggsci)

source("src/exec/LDO_and_GAM_func.r")

# Load data
data <- readRDS("data/lfc_exp_cn.rds")

######################################
########### LDO CORRECTION ###########
######################################

# run Local Drop Out correction
data <- LDO( data,
             params_weigthed_mean = list(
               omega = 1e6,
               decay_func = (function(x)exp(x)),
               in_subunit_weight = 0,
               side = "both"
             ),
             params_regression = list(
               minbucket = 2, 
               minsubunits = 3, 
               cp_init = -3, 
               cp_iteration = 0.1),
             verbose = FALSE)



######################################
########### GAM CORRECTION ###########
######################################

# Smooth out multi-alignment and CNA effect: using generalized additive models
cat("  Smooth out CNA effect using generalized additive models \n")

data2 <- GAM(data, formula = "DEPENDENCY_SCORE ~ CNA", subunit = "SAMPLE_NAME", verbose = FALSE)