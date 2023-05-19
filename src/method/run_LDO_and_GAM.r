library(tidyverse)
library(magrittr)
library(ggsci)

source("src/method/LDO_and_GAM_func.r")


run_method <- function(data, method){
  cell_lines <- unique(data$SAMPLE_NAME)

  # Apply method to all cell lines
  for (i in 1:length(cell_lines)){
    print(paste0("Processing ", cell_lines[i], " ", i, "/", length(cell_lines)))
    cell_line <- cell_lines[i]

    data_sub <- data %>% 
      filter(SAMPLE_NAME == cell_line) %>%
      na.omit()
    
    if (nrow(data_sub) == 0){
      print(paste0("No data for ", cell_line, ", skipping..."))
      next
    }

    if (method == "LDO"){
      ## LDO correction
      data_sub <- LDO(data_sub,
                      params_weigthed_mean = list(
                          omega = 1e6,
                          decay_func = (function(x)exp(x)),
                          in_subunit_weight = 0,
                          side = "both"),
                      params_regression = list(
                          minbucket = 2, 
                          minsubunits = 3, 
                          cp_init = -3, 
                          cp_iteration = 0.1),
                      verbose = FALSE)

      data_sub <- data_sub %>% 
        select(GENE_NAME, DEPENDENCY_SCORE_LDO) %>%
        group_by(GENE_NAME) %>%
        summarise(DEPENDENCY_SCORE_LDO = median(DEPENDENCY_SCORE_LDO, na.rm = TRUE)) %>%
        ungroup() %>%
        rename(Gene = GENE_NAME, !!cell_line := DEPENDENCY_SCORE_LDO)
    } else if (method == "GAM"){
      ## GAM correction
      data_sub <- GAM(data_sub, formula = "DEPENDENCY_SCORE ~ CNA", subunit = "SAMPLE_NAME", verbose = FALSE)
      data_sub <- data_sub %>% 
        select(GENE_NAME, DEPENDENCY_SCORE_GAM) %>%
        group_by(GENE_NAME) %>%
        summarise(DEPENDENCY_SCORE_GAM = median(DEPENDENCY_SCORE_GAM, na.rm = TRUE)) %>%
        ungroup() %>%
        rename(Gene = GENE_NAME, !!cell_line := DEPENDENCY_SCORE_GAM)
    } else {
      stop("Please provide a valid method")
    }
    
    if (i == 1){
      df <- data_sub
    } else {
      ## concatenate results
      df <- full_join(df, data_sub, by = "Gene")
    }
  }

  return(df)
}


# main function
main <- function(combined_data_path, method, output_dir, lib){
  data <- readRDS(combined_data_path)

  ## Run method
  df <- run_method(data = data, method = method)
  write_csv(df, paste0(output_dir, "/", lib, "_gene_", method, ".csv"))
}


# Run the main function
## For instance:
## Rscript src/method/run_LDO_and_GAM.r data/Avana_lfc_exp_cn.rds LDO data/corrected/ Avana
if (!interactive()){
    sys.args <- commandArgs(trailingOnly = TRUE)

    if (length(sys.args) < 4){
        stop(cat("Please provide these args in the following order: 
            \n- combined_data_path: LFC + GuideMap + CN + Exp data
            \n- method: LDO or GAM
            \n- output_dir: path to the output directory
            \n- lib: library name of the output files"))
    }

    combined_data_path <- sys.args[1]
    method <- sys.args[2]
    output_dir <- sys.args[3]
    lib <- sys.args[4]

    main(combined_data_path, method, output_dir, lib)
}