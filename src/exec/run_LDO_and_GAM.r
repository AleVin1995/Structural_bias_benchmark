library(tidyverse)
library(magrittr)
library(ggsci)

source("src/exec/LDO_and_GAM_func.r")

format_omics <- function(df, values_to){
  df <- df %>% 
    pivot_longer(-colnames(.)[1], names_to = 'GENE_NAME', values_to = values_to) %>%
    rename(SAMPLE_NAME = colnames(.)[1]) %>%
    separate(col = "GENE_NAME", sep = " \\(", into = c("GENE_NAME", "code")) %>%
    select(-code)
}

# Load data
data <- readRDS("data/lfc_exp_cn.rds")
# data <- read_csv('data/AvanaGuideMap.csv') %>%
#   filter(UsedByChronos == TRUE) %>%
#   separate(col = "Gene", sep = " \\(", into = c("GENE_NAME", "code")) %>%
#   select(-c(code, nAlignments, DropReason, UsedByChronos)) %>%
#   separate(col = "GenomeAlignment", sep = '_', into = c("CHROMOSOME", "POSITION", "STRAND")) %>%
#   select(-STRAND) %>%
#   mutate(POSITION = as.numeric(POSITION)) %>%
#   inner_join(read_csv('data/raw/Avana_sgrna_raw_LFC.csv') %>%
#     pivot_longer(-sgRNA, names_to = 'SAMPLE_NAME', values_to = 'DEPENDENCY_SCORE'), 
#     by = c("sgRNA")) %>%
#   left_join(read_csv('data/OmicsCNGene.csv') %>%
#     format_omics(., values_to = 'CNA'),
#     by = c("GENE_NAME", "SAMPLE_NAME")) %>%
#   left_join(read_csv('data/OmicsExpressionProteinCodingGenesTPMLogp1.csv') %>%
#     format_omics(., values_to = 'EXP'),
#     by = c("GENE_NAME", "SAMPLE_NAME"))

# Pan-cancer essential genes
essential_genes <- read_csv("data/AchillesCommonEssentialControls.csv") %>% 
  separate(col = "Gene", sep = " \\(", into = c("Gene", "code")) %>%
  pull(Gene)

######################################
########### LDO CORRECTION ###########
######################################

# rename columns
data <- data %>%
  mutate(ESSENTIAL = GENE_NAME %in% essential_genes)

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