library(CRISPRcleanR)
library(tidyverse)

# OncoKB list of oncogenes
oncogenes <- read_tsv('data/biomarkers/cancerGeneList.tsv') %>%
    filter(`Is Oncogene` == 'Yes') %>%
    pull(`Hugo Symbol`) %>%
    unique()

# get list of mutated genes from somatic and fusion data
onco_mut_stat <- read_csv('data/OmicsSomaticMutationsMatrixDamaging.csv') %>%
    dplyr::rename(ModelID = colnames(.)[1]) %>%
    rename_with(~sub(" \\(.*$", "", .x)) %>%
    ## filter oncogenes
    pivot_longer(-ModelID, names_to = "Gene", values_to = "Mutation") %>%
    filter(Gene %in% oncogenes) %>%
    ## merge with fusion data
    full_join(read_csv('data/OmicsFusionFiltered.csv') %>%
        select(ModelID, FusionName) %>%
        separate(FusionName, into = c("Gene1", "Gene2"), sep = "--") %>%
        ## collapse Gene1 and Gene2 into one column
        pivot_longer(-ModelID, names_to = "Group", values_to = "Gene") %>%
        select(-Group) %>%
        ## filter oncogenes
        filter(Gene %in% oncogenes) %>%
        mutate(Fusion = 1)) %>%
    ## fill na with 0
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    mutate(Status = ifelse(Mutation >= 1, 1, ifelse(Fusion == 1, 1, 0))) %>%
    select(-Mutation, -Fusion) %>%
    distinct()


# compute recall
get_recall <- function(x, FDRth = 0.05){
    x <- x %>%
        unite("Gene_ModelID", Gene:ModelID)

    pos <- x %>%
        filter(Status == 1) %>%
        pull(Gene_ModelID)
    
    neg <- x %>%
        filter(Status == 0) %>%
        pull(Gene_ModelID)

    vec <- x %>%
        pull(LFC, name = "Gene_ModelID")
    
    res <- ccr.ROC_Curve(vec, pos, neg, display = FALSE, FDRth = FDRth)$AUC %>% 
        as.numeric()
    
    return(res)
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
            mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
            dplyr::rename(Gene = colnames(.)[1])) %>% ## fill na with 0
        set_names(dfs_names)
    
    ## common cell lines/genes
    common_cells <- Reduce(intersect, map(dfs, ~colnames(.)[2:length(colnames(.))]))

    common_genes <- Reduce(intersect, map(dfs, ~.$Gene))
    common_genes <- intersect(common_genes, oncogenes)

    ## select only common cell lines and first columns
    dfs <- map(dfs, ~select(., colnames(.)[1], all_of(common_cells))) %>% 
        map(~filter(., Gene %in% common_genes) %>%
            pivot_longer(-1, names_to = "ModelID", values_to = "LFC")) %>%
        bind_rows(.id = "Algorithm")
    
    ## merge with somatic and fusion data
    dfs <- inner_join(dfs, onco_mut_stat, by = c("ModelID", "Gene")) %>%
        ## compute mean LFC
        group_by(Algorithm, Gene) %>%
        mutate(mean_LFC = mean(LFC)) %>%
        ungroup() %>%
        ## filter genes with mean LFC >= -0.5
        filter(mean_LFC >= -0.5) %>%
        group_by(Algorithm, Gene) %>%
        ## occurrence of mutated/non-mutated genes
        mutate(wt_occ = sum(Status == 0),
               mut_occ = sum(Status == 1)) %>%
        ungroup() %>%
        ## filter genes with at least 3 mutated and 3 non-mutated cell lines
        filter(wt_occ >= 3 & mut_occ >= 3)
    
    ## compute AUROC oncogenes
    onco_auc <- dfs %>%
        group_split(Algorithm) %>%
        map(~.x %>%
            mutate(AUROC = get_recall(.))) %>%
        bind_rows() %>%
        select(Algorithm, AUROC) %>%
        distinct()
    

    ## save results
    saveRDS(onco_auc, paste0("results/analyses/impact_data_quality/", lib, "_onco_auc.rds"))
}
