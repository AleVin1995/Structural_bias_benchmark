library(tidyverse)

# load absolute CN data
CN_abs <- read_csv('data/cnv_abs_copy_number_picnic_20191101.csv') %>%
    select(-model_id) %>%
    rename(Gene = colnames(.)[1]) %>%
    ## remove rows
    dplyr::slice(-2) %>%
    replace_na(list(Gene = "Gene")) %>%
    ## replace all with character
    mutate(across(where(is.numeric), as.character)) %>%
    set_names(.[1,]) %>%
    dplyr::slice(-1) %>%
    pivot_longer(-Gene, names_to = "CellLineName", values_to = "CN_abs") %>%
    ## convert cell line names to model IDs
    inner_join(read_csv('data/Model.csv') %>%
        select(ModelID, CellLineName) %>%
        distinct()) %>%
    select(-CellLineName) %>%
    mutate(CN_abs = as.numeric(CN_abs)) %>%
    drop_na()

CN_abs_tpm <- CN_abs %>%
    inner_join(read_csv("data/OmicsExpressionProteinCodingGenesTPMLogp1.csv") %>%
        dplyr::rename(ModelID = colnames(.)[1]) %>%
        rename_with(~sub(" \\(.*$", "", .x)) %>%
        pivot_longer(-ModelID, names_to = "Gene", values_to = "TPM")) %>%
    drop_na() %>%
    group_by(Gene) %>%
    mutate(mean_TPM = mean(TPM, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(mean_TPM < 1) %>%
    select(-mean_TPM)


# Define list of algorithms and libraries
algos <- c("CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK")
libs <- c("Avana", "KY")

# iterate over algorithms and libraries
for (lib in libs){
    dfs_names <- c("Uncorrected", algos)

    dfs <- paste0("data/corrected/", lib, "_gene_", algos, ".csv") %>%
        c(paste0("data/raw/", lib, "_gene_raw_LFC.csv"), .) %>%
        map(~.x %>%
            read_csv %>%
            mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
            dplyr::rename(Gene = colnames(.)[1])) %>% ## fill na with 0
        set_names(dfs_names)

    ## merge with CN abs data
    dfs <- map(dfs, ~.x %>%
            pivot_longer(-1, names_to = "ModelID", values_to = "LFC")) %>%
        bind_rows(.id = "Algorithm") %>%
        inner_join(CN_abs, by = c("ModelID", "Gene"))
    
    ## merge with TPM data
    dfs_tpm <- dfs %>%
            inner_join(CN_abs_tpm, by = c("ModelID", "Gene", "CN_abs"))

    ## save results
    saveRDS(dfs, paste0("results/analyses/cn_correction/", lib, "_cn_abs.rds"))
    saveRDS(dfs_tpm, paste0("results/analyses/cn_correction/", lib, "_cn_abs_tpm.rds"))
}
