library(tidyverse)

get_sig_biomarkers <- function(biomarkers, dfs){    
    ## library-specific biomarkers
    biomarkers_lib <- biomarkers %>%
        filter(ModelID %in% common_cells) %>% 
        ## at least 3 samples per CFE with present/absent status
        group_by(CFE, DepmapModelType, Status) %>%
        mutate(sample_size = n()) %>%
        ungroup() %>%
        group_by(CFE, DepmapModelType) %>%
        mutate(Occurrence = n()) %>%
        filter(Occurrence-sample_size >= 3 & sample_size >= 3) %>%
        ungroup() %>%
        select(-sample_size, -Occurrence)

    ## join biomarkers with gene LFC
    biomarkers_lib <- map(dfs, ~.x %>%
            inner_join(biomarkers_lib, by = "ModelID", relationship = "many-to-many")) %>%
        bind_rows(.id = "Algorithm") %>%
        ## split by factor levels
        group_split(Algorithm, Gene, CFE, DepmapModelType) 
    
    n_associations <- length(biomarkers_lib)
    res_lib_sig <- list()

    for (i in seq(n_associations)){
        res_lib_sig[[i]] <- biomarkers_lib[[i]] %>%
            mutate(pvalue = t.test(LFC~Status, data=.)$p.value)

        if (i %% 10000 == 0){
            print(paste0(i, "/", n_associations))
        }
    }
    
    sig_biomarkers <- res_lib_sig %>%
        bind_rows() %>%
        group_by(Algorithm) %>%
        mutate(FDR = p.adjust(pvalue, method = "fdr")) %>%
        filter(FDR < 0.05) %>%
        mutate(n_sig_biomark = n()) %>%
        ungroup() %>%
        select(Algorithm, n_sig_biomark) %>%
        distinct()
    
    return(sig_biomarkers)
}

# load data
biomarkers <- readRDS("data/biomarkers/biomarkers.rds")
oncogenes <- read_tsv("data/biomarkers/cancerGeneList.tsv") %>%
    filter(`Is Oncogene` == "Yes") %>%
    pull(`Hugo Symbol`) %>%
    unique()
ampl_racs <- biomarkers %>%
    filter(grepl("gain", CFE)) %>%
    select(CFE) %>%
    distinct() %>%
    separate(CFE, c("CFE", "Gene"), sep = " \\(") %>%
    select(Gene) %>%
    separate(Gene, c("Gene", "Else"), sep = "\\)") %>%
    select(Gene) %>%
    na.omit() %>%
    separate_rows(Gene, sep = ",") %>%
    pull(Gene) %>%
    unique()
ssd <- read_csv("data/biomarkers/skewed_tdist.csv") %>%
    ## select strong selective dependencies
    filter(Skewed_tdist > 100) %>%
    pull(HUGO)

# Define list of algorithms and libraries
algos <- c("CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO")
libs <- c("Avana", "KY")

# iterate over algorithms and libraries
for (lib in libs){
    dfs_names <- c("Uncorrected", algos)

    dfs <- paste0("data/corrected/", lib, "_gene_", algos, ".csv") %>%
        c(paste0("data/raw/", lib, "_gene_raw_LFC.csv"), .) %>%
        map(~.x %>%
            read_csv %>%
            mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
            rename(Gene = colnames(.)[1])) %>% ## fill na with 0
        set_names(dfs_names)
    
    ## common cell lines/genes
    common_cells <- Reduce(intersect, map(dfs, ~colnames(.)[2:length(colnames(.))]))
    common_cells <- intersect(common_cells, unique(biomarkers$ModelID))

    common_genes <- Reduce(intersect, map(dfs, ~.$Gene))
    ssd <- intersect(common_genes, ssd)
    
    ## select only common cell lines and first columns
    dfs_ssd <- map(dfs, ~select(., colnames(.)[1], all_of(common_cells))) %>% 
        map(~filter(., Gene %in% ssd) %>%
            pivot_longer(-1, names_to = "ModelID", values_to = "LFC"))
    
    dfs_onco <- map(dfs, ~select(., colnames(.)[1], all_of(common_cells))) %>% 
        map(~filter(., Gene %in% oncogenes | Gene %in% ampl_racs) %>%
            pivot_longer(-1, names_to = "ModelID", values_to = "LFC"))
    
    ## library-specific biomarkers
    sig_biomarkers_ssd <- get_sig_biomarkers(biomarkers, dfs_ssd)
    sig_biomarkers_onco <- get_sig_biomarkers(biomarkers, dfs_onco)

    ## save results
    saveRDS(sig_biomarkers_ssd, paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_ssd.rds"))
    saveRDS(sig_biomarkers_onco, paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_onco.rds"))
}
