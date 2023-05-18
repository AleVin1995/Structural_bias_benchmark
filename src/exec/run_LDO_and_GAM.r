library(tidyverse)
library(magrittr)
library(ggsci)

source("src/exec/LDO_and_GAM_func.r")

######################################
########### LDO CORRECTION ###########
######################################

# rename columns
data <- data %>%
  mutate( ESSENTIAL = GENESYMBOL %in% essential_genes) %>% 
  dplyr::rename(DEPENDENCY_SCORE = LogFC_org, GENE_NAME = GENESYMBOL, SAMPLE_NAME = CLEANNAME, CHROMOSOME = CHR)

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

# Run the main function
## For instance:
## Rscript src/exec/run_CCR.r data/raw/Avana_sgrna_raw_LFC.csv data/AvanaGuideMap.csv data/corrected/ Avana
if (!interactive()){
    sys.args <- commandArgs(trailingOnly = TRUE)

    if (length(sys.args) < 4){
        stop(cat("Please provide these args in the following order: 
            \n- raw_LFC_path: path to the raw LFC file
            \n- GuideMap_path: path to the GuideMap file
            \n- output_dir: path to the output directory
            \n- lib: library name of the output files"))
    }

    raw_LFC_path <- sys.args[1]
    GuideMap_path <- sys.args[2]
    output_dir <- sys.args[3]
    lib <- sys.args[4]

    # Load data
    data <- read_csv('data/raw/Avana_sgrna_raw_LFC.csv') %>%
      column_to_rownames('sgRNA')
    guide_map <- read_csv('data/AvanaGuideMap.csv')


    # Pan-cancer essential genes
    essential_genes <- read_csv("data/AchillesCommonEssentialControls.csv") %>% 
      separate(col = "Gene", sep = " \\(", into = c("Gene", "code")) %>%
      pull(Gene)

    main(raw_LFC_path, GuideMap_path, output_dir, lib)
}