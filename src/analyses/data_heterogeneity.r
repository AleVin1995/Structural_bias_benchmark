library(tidyverse)

# load data
biomarkers <- readRDS("data/biomarkers/biomarkers.rds")
ssd <- read_csv("data/biomarkers/skewed_tdist.csv") %>%
    ## select strong selective dependencies
    filter(Skewed_tdist > 100) %>%
    pull(HUGO)

# Define list of algorithms and libraries
algos <- c("CCR", "CERES", "Chronos", "Crispy", "GAM", "Geometric", "LDO")
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
    common_genes <- intersect(common_genes, ssd)

    ssd <- common_genes

    ## select only common cell lines and first columns
    dfs <- map(dfs, ~select(., colnames(.)[1], all_of(common_cells))) %>% 
        map(~filter(., Gene %in% common_genes) %>%
            pivot_longer(-1, names_to = "ModelID", values_to = "LFC"))
    
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
        bind_rows(.id = "Algorithm")
    
    ## identification of biomarkers
    res_lib <- biomarkers_lib %>%
        group_split(Algorithm, Gene, CFE, DepmapModelType) 
    
    n_associations <- length(res_lib)
    res_lib_sig <- list()

    for (i in seq(n_associations)){
        res_lib_sig[[i]] <- res_lib[[i]] %>%
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
    

    ## save results
    p_sig_biomark <- ggplot(sig_biomarkers, aes(x = Algorithm, y = n_sig_biomark, fill = Algorithm)) +
        geom_bar(stat = "identity") +
        labs(x = "Method", y = "Nº significant associations", 
            title = "Nº significant associations") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_sig_biomark, filename = paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers.pdf"), width = 10, height = 10)
    saveRDS(sig_biomarkers, paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers.rds"))
}
