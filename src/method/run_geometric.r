library(tidyverse)


# add arm information
arm_info <- function(df, cytoband_path){
    genes <- unique(df$GENE_NAME)

    ## get cytoband information
    cyto_band <- read_csv(cytoband_path) %>% 
        filter(symbol %in% genes) %>% 
        distinct() %>% 
        as_tibble() %>% 
        mutate(ARM = ifelse(grepl("p", map), "p", ifelse(grepl("q", map), "q", NA))) %>%
        rename("GENE_NAME" = "symbol") %>%
        group_by(GENE_NAME) %>%
        mutate(COUNT = n()) %>% ## remove gene with contradictory arm information
        ungroup() %>%
        filter(COUNT == 1) %>%
        filter(!is.na(ARM)) %>%
        dplyr::select(GENE_NAME, ARM) %>%
        distinct()
    
    df <- inner_join(df, cyto_band, by = "GENE_NAME")

    return(df)
}


# compute average LFC for unexpressed genes in each arm
avg_lfc_unexp <- function(df, thr = 1){
    arm_base_lfc <- df %>%
        dplyr::select(CHROMOSOME, ARM, EXP, DEPENDENCY_SCORE) %>%
        filter(EXP < thr) %>%
        group_by(CHROMOSOME, ARM) %>%
        summarise(MEAN_LFC_UNEXP = mean(DEPENDENCY_SCORE, na.rm = TRUE)) %>%
        ungroup()

    return(arm_base_lfc)
}


# Run geometric correction for proximity bias
run_geometric <- function(data){
    cell_lines <- unique(data$SAMPLE_NAME)

    # Apply method to all cell lines
    for (i in 1:length(cell_lines)){
        print(paste0("Processing ", cell_lines[i], " ", i, "/", length(cell_lines)))
        cell_line <- cell_lines[i]

        data_sub <- data %>% 
            filter(SAMPLE_NAME == cell_line)

        if (nrow(data_sub) == 0){
            print(paste0("No data for ", cell_line, ", skipping..."))
            next
        }

        ## select unexpressed genes
        data_sub <- data_sub %>%
            left_join(avg_lfc_unexp(data_sub), by = c("CHROMOSOME", "ARM")) %>%
            mutate(DEPENDENCY_SCORE_CORRECTED = DEPENDENCY_SCORE - MEAN_LFC_UNEXP) %>%
            dplyr::select(GENE_NAME, DEPENDENCY_SCORE_CORRECTED) %>%
            group_by(GENE_NAME) %>%
            summarise(DEPENDENCY_SCORE_CORRECTED = median(DEPENDENCY_SCORE_CORRECTED, na.rm = TRUE)) %>%
            rename(Gene = GENE_NAME, !!cell_line := DEPENDENCY_SCORE_CORRECTED) %>%
            ungroup()

        if (i == 1){
            df <- data_sub
        } else {
            df <- left_join(df, data_sub, by = "GENE_NAME")
        }
    }

    return(df)
}


# main function
main <- function(combined_data_path, cytoband_path, output_dir, lib){
    data <- readRDS(combined_data_path)
    data <- arm_info(data, cytoband_path)

    ## run geometric correction
    df <- run_geometric(data)
    write_csv(df, paste0(output_dir, "/", lib, "_gene_Geometric.csv"))
}


# Run main function
## For instance:
## Rscript src/method/run_geometric.r data/Avana_lfc_exp_cn.rds data/corrected/ Avana
if (!interactive()){
    sys.args <- commandArgs(trailingOnly = TRUE)

    if (length(sys.args) < 4){
        stop(cat("Please provide these args in the following order: 
            \n- combined_data_path: LFC + GuideMap + CN + Exp data
            \n- cytoband_path: path to the cytoband mapping file
            \n- output_dir: path to the output directory
            \n- lib: library name of the output files"))
    }

    combined_data_path <- sys.args[1]
    cytoband_path <- sys.args[2]
    output_dir <- sys.args[3]
    lib <- sys.args[4]

    main(combined_data_path, cytoband_path, output_dir, lib)
}