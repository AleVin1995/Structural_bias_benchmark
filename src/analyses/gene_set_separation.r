library(tidyverse)

# build predefined sets of genes
ess_genes <- read_csv("data/AchillesCommonEssentialControls.csv") %>%
    separate(Gene, into = c("Gene", "code"), sep = " \\(") %>%
    pull(Gene) %>%
    unique()

noness_genes <- read_csv("data/AchillesNonessentialControls.csv") %>%
    separate(Gene, into = c("Gene", "code"), sep = " \\(") %>%
    pull(Gene) %>%
    unique()


# compute gene sets separation
compute_gene_sets_separation <- function(gene_lfc, ess_genes, noness_genes){
    ess_lfc <- gene_lfc[intersect(names(gene_lfc), ess_genes)]
    noness_lfc <- gene_lfc[intersect(names(gene_lfc), noness_genes)]

    ## compute separation between essential and non-essential genes
    separation <- (median(ess_lfc, na.rm = TRUE) - median(noness_lfc, na.rm = TRUE))/mad(noness_lfc, na.rm = TRUE)
    return(separation)
}


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
            mutate(across(where(is.numeric), ~replace_na(., 0)))) %>% ## fill na with 0
        set_names(dfs_names)
    
    ## common cell lines
    common_cells <- Reduce(intersect, map(dfs, ~colnames(.)[2:length(colnames(.))]))

    ## select only common cell lines and first columns
    dfs <- map(dfs, ~select(., colnames(.)[1], all_of(common_cells)))

    ## compute separation between essential and non-essential genes
    gene_sep <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                rename(Gene = colnames(.)[1]) %>%
                split(.$ModelID) %>%
                map(~.x %>% 
                    pull(LFC, name = "Gene") %>% 
                    compute_gene_sets_separation(., ess_genes, noness_genes) %>%
                    as.numeric()) %>%
                bind_rows() %>%
                pivot_longer(everything(.), names_to = "ModelID", values_to = "Separation")) %>%
        set_names(dfs_names) %>%
        bind_rows(.id = "Algorithm") %>%
        ### replace NA or Inf with 0
        mutate(Separation = replace_na(Separation, 0)) %>%
        mutate(Separation = replace(Separation, is.infinite(Separation), 0))
    

    ## save results
    saveRDS(gene_sep, paste0("results/analyses/impact_data_quality/", lib, "_gene_sets_separation.rds"))
}
