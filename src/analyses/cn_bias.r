library(tidyverse)

# load absolute CN data
CN_abs <- read_csv('data/cnv_abs_copy_number_picnic_20191101.csv') %>%
    select(-model_id) %>%
    dplyr::slice(-2) %>% ## remove rows
    t() %>%
    as_tibble() %>%
    mutate_all(~replace_na(., "0")) %>%
    column_to_rownames(colnames(.)[1]) %>%
    t() %>%
    as_tibble() %>%
    column_to_rownames(colnames(.)[1]) %>%
    t() %>%
    as.data.frame() %>%
    mutate_all(as.numeric) %>%
    rownames_to_column("CellLineName") %>%
    pivot_longer(-CellLineName, names_to = "Gene", values_to = "CN_abs") %>%
    ## convert cell line names to model IDs
    inner_join(read_csv('data/Model.csv') %>%
        select(ModelID, CellLineName)) %>%
    select(-CellLineName)


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
    
    ## common cell lines/genes
    common_cells <- Reduce(intersect, map(dfs, ~colnames(.)[2:length(colnames(.))]))
    common_cells <- intersect(common_cells, CN_abs$ModelID)

    common_genes <- Reduce(intersect, map(dfs, ~.$Gene))
    common_genes <- intersect(common_genes, CN_abs$Gene)

    ## select only common cell lines and first columns
    dfs <- map(dfs, ~select(., colnames(.)[1], all_of(common_cells))) %>% 
        map(~filter(., Gene %in% common_genes) %>%
            pivot_longer(-1, names_to = "ModelID", values_to = "LFC")) %>%
        bind_rows(.id = "Algorithm")
    
    ## merge with CN abs data
    dfs <- inner_join(dfs, CN_abs, by = c("ModelID", "Gene"))

    ## save results
    saveRDS(dfs, paste0("results/analyses/cn_correction/", lib, "_cn_abs.rds"))
}