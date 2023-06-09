library(brunnermunzel)
library(tidyverse)

# genes involved in DNA repair
repair_genes <- c("TP53", "CDKN2A", "CDKN2B", "CDKN2C", "BTG2")

# load cytoband information
cytoband <- read_csv('data/cytoband_mapping.csv') %>%
    rename(Gene = symbol) %>%
    mutate(Arm = ifelse(grepl("p", map), "p", ifelse(grepl("q", map), "q", NA))) %>%
    separate(map, into = c("Chromosome", "Band"), sep = "p|q") %>%
    na.omit() %>%
    select(-Band) %>%
    distinct()

# load somatic mutations
mut_stat <- read_csv('data/OmicsSomaticMutationsMatrixDamaging.csv') %>%
    rename(ModelID = colnames(.)[1]) %>%
    column_to_rownames("ModelID") %>%
    t() %>%
    as_tibble(rownames = "Gene") %>%
    separate(Gene, into = c("Gene", "code"), sep = " \\(") %>%
    select(-code) %>%
    ## filter repair genes
    filter(Gene %in% repair_genes) %>%
    pivot_longer(-Gene, names_to = "ModelID", values_to = "Mutation") %>%
    ## fill na with 0
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    mutate(Status = ifelse(Mutation >= 1, 1, 0)) %>%
    select(-Mutation) %>%
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
            # print(paste0("Processing ", chrm, " ", arm))
            df_sub <- df_corr %>%
                filter(Chromosome1 == chrm | Chromosome2 == chrm) %>%
                filter(Arm1 == arm | Arm2 == arm)
            
            if (nrow(df_sub) == 0){
                # print(paste0("No data for ", chrm, ", skipping..."))
                next
            }

            if (length(unique(df_sub$Class)) != 2){
                # print(paste0("Not two classes for chromosome ", chrm, " and arm ", arm, ", skipping..."))
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


# split based on mutation status of a gene
mut_stat_split <- function(df, mut, cytoband_info){
    ## get gene of interest
    goi <- unique(mut$Gene)

    ## split based on mutation status
    df_wt <- df %>%
        select(Gene, any_of(mut %>% 
            filter(Status == 0) %>% 
            pull(ModelID)))
    
    df_mut <- df %>%
        select(Gene, any_of(mut %>% 
            filter(Status == 1) %>% 
            pull(ModelID)))
    
    ## perform brunner-munzel test (at least one cell line per group)
    if (ncol(df_wt) >= 2 & ncol(df_mut) >= 2){
        res_wt <- perform_brunnermunzel(df_wt, cytoband_info) %>%
            mutate(Status = "WT")
        res_mut <- perform_brunnermunzel(df_mut, cytoband_info) %>%
            mutate(Status = "Mut")

        res <- bind_rows(res_wt, res_mut)

        return(res)
    } else {
        print(paste0("No data for ", goi, ", skipping..."))
    }
}


# Define list of algorithms and libraries
algos <- c("CCR", "CERES", "Chronos", "Crispy", "GAM", "Geometric", "LDO")

# parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0){
    libs <- c("Avana", "KY")
} else {
    libs <- args[[1]]
    print(paste0("Running for ", libs))
}

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
    print(paste0("Performing brunner-munzel test for ", lib, " on pooled..."))
    bm_pool <- dfs %>%
        map(~.x %>%
            perform_brunnermunzel(., cytoband)) %>%
        bind_rows(.id = "Algorithm")
    
    ## save results
    p_pool <- bm_pool %>%
        unite("Coord", c("Chromosome", "Arm"), sep = "") %>%
        ggplot(aes(x = Algorithm, y = est, fill = Algorithm)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "Genome-wide") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1) +
        facet_wrap(~Coord, scales = "free_x")
    ggsave(p_pool, paste0("results/analyses/proximity_bias/", lib, "_bm_pool.pdf"), width = 10, height = 10, dpi = 300)
    saveRDS(bm_pool, paste0("results/analyses/proximity_bias/", lib, "_bm_pool.rds"))

    ## perform brunner-munzel test for repair genes
    if ("TP53" %in% common_genes){
        print(paste0("Performing brunner-munzel test for ", lib, " on TP53..."))
        bm_TP53 <- dfs %>%
            map(~.x %>%
                mut_stat_split(., 
                    mut_stat %>%
                    filter(Gene == "TP53"), 
                    cytoband)) %>%
            bind_rows(.id = "Algorithm")
        
        ### save results
        p_TP53 <- bm_TP53 %>%
            unite("Coord", c("Chromosome", "Arm"), sep = "") %>%
            ggplot(aes(x = Algorithm, y = est, fill = Status)) +
            geom_bar(stat = "identity", position=position_dodge()) +
            geom_hline(yintercept = 0.5, linetype = "dashed") +
            labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "TP53 mutational status") +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size = 10, color = 'black'),
                axis.title = element_text(size = 12),
                plot.title = element_text(size = 14, hjust = 0.5),
                plot.subtitle = element_text(size = 12, hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1),
                aspect.ratio = 1) +
            facet_wrap(~Coord, scales = "free_x")
        ggsave(p_TP53, paste0("results/analyses/proximity_bias/", lib, "_bm_TP53.pdf"), width = 10, height = 10, dpi = 300)
        saveRDS(bm_TP53, paste0("results/analyses/proximity_bias/", lib, "_bm_TP53.rds"))
    }
    
    if ("CDKN2A" %in% common_genes){
        print(paste0("Performing brunner-munzel test for ", lib, " on CDKN2A..."))
        bm_CDKN2A <- dfs %>%
            map(~.x %>%
                mut_stat_split(., 
                    mut_stat %>%
                    filter(Gene == "CDKN2A"), 
                    cytoband)) %>%
            bind_rows(.id = "Algorithm")
        
        ### save results
        p_CDKN2A <- bm_CDKN2A %>%
            unite("Coord", c("Chromosome", "Arm"), sep = "") %>%
            ggplot(aes(x = Algorithm, y = est, fill = Status)) +
            geom_bar(stat = "identity", position=position_dodge()) +
            geom_hline(yintercept = 0.5, linetype = "dashed") +
            labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "CDKN2A mutational status") +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size = 10, color = 'black'),
                axis.title = element_text(size = 12),
                plot.title = element_text(size = 14, hjust = 0.5),
                plot.subtitle = element_text(size = 12, hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1),
                aspect.ratio = 1) +
            facet_wrap(~Coord, scales = "free_x")
        ggsave(p_CDKN2A, paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2A.pdf"), width = 10, height = 10, dpi = 300)
        saveRDS(bm_CDKN2A, paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2A.rds"))
    }
    
    if ("CDKN2B" %in% common_genes){
        print(paste0("Performing brunner-munzel test for ", lib, " on CDKN2B..."))
        bm_CDKN2B <- dfs %>%
            map(~.x %>%
                mut_stat_split(., 
                    mut_stat %>%
                    filter(Gene == "CDKN2B"), 
                    cytoband)) %>%
            bind_rows(.id = "Algorithm")
        
        ### save results
        p_CDKN2B <- bm_CDKN2B %>%
            unite("Coord", c("Chromosome", "Arm"), sep = "") %>%
            ggplot(aes(x = Algorithm, y = est, fill = Status)) +
            geom_bar(stat = "identity", position=position_dodge()) +
            geom_hline(yintercept = 0.5, linetype = "dashed") +
            labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "CDKN2B mutational status") +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size = 10, color = 'black'),
                axis.title = element_text(size = 12),
                plot.title = element_text(size = 14, hjust = 0.5),
                plot.subtitle = element_text(size = 12, hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1),
                aspect.ratio = 1) +
            facet_wrap(~Coord, scales = "free_x")
        ggsave(p_CDKN2B, paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2B.pdf"), width = 10, height = 10, dpi = 300)
        saveRDS(bm_CDKN2B, paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2B.rds"))
    }
    
    if ("CDKN2C" %in% common_genes){
        print(paste0("Performing brunner-munzel test for ", lib, " on CDKN2C..."))
        bm_CDKN2C <- dfs %>%
            map(~.x %>%
                mut_stat_split(., 
                    mut_stat %>%
                    filter(Gene == "CDKN2C"), 
                    cytoband)) %>%
            bind_rows(.id = "Algorithm")
        
        p_CDKN2C <- bm_CDKN2C %>%
            unite("Coord", c("Chromosome", "Arm"), sep = "") %>%
            ggplot(aes(x = Algorithm, y = est, fill = Status)) +
            geom_bar(stat = "identity", position=position_dodge()) +
            geom_hline(yintercept = 0.5, linetype = "dashed") +
            labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "CDKN2C mutational status") +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size = 10, color = 'black'),
                axis.title = element_text(size = 12),
                plot.title = element_text(size = 14, hjust = 0.5),
                plot.subtitle = element_text(size = 12, hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1),
                aspect.ratio = 1) +
            facet_wrap(~Coord, scales = "free_x")
        ggsave(p_CDKN2C, paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2C.pdf"), width = 10, height = 10, dpi = 300)
        saveRDS(bm_CDKN2C, paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2C.rds"))
    }
    
    if ("BTG2" %in% common_genes){
        print(paste0("Performing brunner-munzel test for ", lib, " on BTG2..."))
        bm_BTG2 <- dfs %>%
            map(~.x %>%
                mut_stat_split(., 
                    mut_stat %>%
                    filter(Gene == "BTG2"), 
                    cytoband)) %>%
            bind_rows(.id = "Algorithm")
        
        ### save results
        p_BTG2 <- bm_BTG2 %>%
            unite("Coord", c("Chromosome", "Arm"), sep = "") %>%
            ggplot(aes(x = Algorithm, y = est, fill = Status)) +
            geom_bar(stat = "identity", position=position_dodge()) +
            geom_hline(yintercept = 0.5, linetype = "dashed") +
            labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "BTG2 mutational status") +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size = 10, color = 'black'),
                axis.title = element_text(size = 12),
                plot.title = element_text(size = 14, hjust = 0.5),
                plot.subtitle = element_text(size = 12, hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1),
                aspect.ratio = 1) +
            facet_wrap(~Coord, scales = "free_x")
        ggsave(p_BTG2, paste0("results/analyses/proximity_bias/", lib, "_bm_BTG2.pdf"), width = 10, height = 10, dpi = 300)
        saveRDS(bm_BTG2, paste0("results/analyses/proximity_bias/", lib, "_bm_BTG2.rds"))
    }
}