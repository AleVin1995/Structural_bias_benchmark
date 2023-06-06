library(brunnermunzel)
library(tidyverse)

# load cytoband information
cytoband <- read_csv('data/cytoband_mapping.csv') %>%
    rename(Gene = symbol) %>%
    mutate(Arm = ifelse(grepl("p", map), "p", ifelse(grepl("q", map), "q", NA))) %>%
    separate(map, into = c("Chromosome", "Band"), sep = "p|q") %>%
    na.omit() %>%
    select(-Band) %>%
    distinct()


# perform brunner-munzel test
perform_brunnermunzel <- function(df, cytoband_info){
    df_corr <- df %>% pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
            group_by(ModelID) %>%
            mutate(LFC = LFC - mean(LFC, na.rm = TRUE)) %>%
            ungroup() %>%
            pivot_wider(names_from = ModelID, values_from = LFC) %>%
            column_to_rownames("Gene") %>%
            t() %>%
            WGCNA::cor(use = "pairwise.complete.obs") %>%
            as_tibble(rownames = "Gene1") %>%
            pivot_longer(-Gene1, names_to = "Gene2", values_to = "Cor")
    
    ## merge with cytoband information
    df_corr <- df_corr %>%
        inner_join(cytoband_info, by = c("Gene1" = "Gene")) %>%
        rename(Chromosome1 = Chromosome, Arm1 = Arm) %>%
        inner_join(cytoband_info, by = c("Gene2" = "Gene")) %>%
        rename(Chromosome2 = Chromosome, Arm2 = Arm) %>%
        mutate(Class = ifelse(Chromosome1 == Chromosome2 & Arm1 == Arm2, "Intra", "Inter"))
    
    ## perform brunner-munzel test
    chrms <- unique(df_corr$Chromosome1)
    arms <- unique(df_corr$Arm1)

    res <- data.frame(Chromosome = rep(chrms, length(arms)), Arm = rep(arms, each = length(chrms)), est = NA) %>%
        arrange(Chromosome, Arm)

    for (chrm in chrms){
        for (arm in arms){
            print(paste0("Processing ", chrm, " ", arm))
            df_sub <- df_corr %>%
                filter(Chromosome1 == chrm | Chromosome2 == chrm) %>%
                filter(Arm1 == arm | Arm2 == arm)
            
            if (nrow(df_sub) == 0){
                print(paste0("No data for ", chrm, ", skipping..."))
                next
            }

            if (length(unique(df_sub$Class)) != 2){
                print(paste0("Not two classes for chromosome ", chrm, " and arm ", arm, ", skipping..."))
                next
            }
            
            res[res$Chromosome == chrm & res$Arm == arm, "est"] <- brunnermunzel.test(Cor ~ Class, data = df_sub)$estimate %>%
                as.numeric()
        }
    }

    res <- res %>%
        na.omit()

    return(res)
}


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
            dplyr::rename(Gene = colnames(.)[1])) %>% ## fill na with 0
        set_names(dfs_names)
    
    ## common cell lines/genes
    common_cells <- Reduce(intersect, map(dfs, ~colnames(.)[2:length(colnames(.))]))

    common_genes <- Reduce(intersect, map(dfs, ~.$Gene))

    ## select only common cell lines/genes and compute pairwise correlation
    dfs <- map(dfs, ~select(., colnames(.)[1], all_of(common_cells))) %>% 
        map(~filter(., Gene %in% common_genes))
    
    ## perform brunner-munzel test
    bm_pool <- dfs %>%
        map(~.x %>%
            perform_brunnermunzel(., cytoband)) %>%
        bind_rows(.id = "Algorithm")
    
    saveRDS(bm_pool, paste0('test/', lib, "_bm_pool.rds"))

    ## perform brunner-munzel test for repair genes
    bm_TP53 <- dfs %>%
        map(~.x %>%
            filter(Gene == "TP53") %>%
            perform_brunnermunzel(., cytoband)) %>%
        bind_rows(.id = "Algorithm")
    
    saveRDS(bm_TP53, paste0('test/', lib, "_bm_TP53.rds"))
    
    bm_CDKN2A <- dfs %>%
        map(~.x %>%
            filter(Gene == "CDKN2A") %>%
            perform_brunnermunzel(., cytoband)) %>%
        bind_rows(.id = "Algorithm")
    
    saveRDS(bm_CDKN2A, paste0('test/', lib, "_bm_CDKN2A.rds"))
    
    bm_CDKN2B <- dfs %>%
        map(~.x %>%
            filter(Gene == "CDKN2B") %>%
            perform_brunnermunzel(., cytoband)) %>%
        bind_rows(.id = "Algorithm")
    
    saveRDS(bm_CDKN2B, paste0('test/', lib, "_bm_CDKN2B.rds"))
    
    bm_CDKN2C <- dfs %>%
        map(~.x %>%
            filter(Gene == "CDKN2C") %>%
            perform_brunnermunzel(., cytoband)) %>%
        bind_rows(.id = "Algorithm")
    
    saveRDS(bm_CDKN2C, paste0('test/', lib, "_bm_CDKN2C.rds"))
    
    bm_BTG2 <- dfs %>%
        map(~.x %>%
            filter(Gene == "BTG2") %>%
            perform_brunnermunzel(., cytoband)) %>%
        bind_rows(.id = "Algorithm")
    
    saveRDS(bm_BTG2, paste0('test/', lib, "_bm_BTG2.rds"))
}