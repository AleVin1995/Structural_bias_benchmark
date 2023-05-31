library(CRISPRcleanR)
library(tidyverse)


# load essential genes from CRISPRcleanR (MSigDB)
data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.HISTONES) 
data(EssGenes.KEGG_rna_polymerase)   
data(EssGenes.PROTEASOME_cons)    
data(EssGenes.SPLICEOSOME_cons)     
data(EssGenes.ribosomalProteins)


# build predefined sets of genes
ess_genes <- read_csv("data/AchillesCommonEssentialControls.csv") %>%
    separate(Gene, into = c("Gene", "code"), sep = " \\(") %>%
    pull(Gene) %>%
    unique()
noness_genes <- read_csv("data/AchillesNonessentialControls.csv") %>%
    separate(Gene, into = c("Gene", "code"), sep = " \\(") %>%
    pull(Gene) %>%
    unique()
msigdb_genes <- c(EssGenes.DNA_REPLICATION_cons,
        EssGenes.HISTONES,
        EssGenes.KEGG_rna_polymerase,
        EssGenes.PROTEASOME_cons,
        EssGenes.SPLICEOSOME_cons,
        EssGenes.ribosomalProteins) %>%
    unique()
cn_ampl_genes <- read_csv("data/OmicsCNGene.csv") %>%
    rename(ModelID = colnames(.)[1]) %>%
    pivot_longer(-ModelID, names_to = "Gene", values_to = "CN_ratio") %>%
    group_by(Gene) %>%
    mutate(mean_CN_ratio = mean(CN_ratio, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(mean_CN_ratio)) %>%
    select(Gene) %>%
    distinct() %>%
    slice(1:round(nrow(.)*0.01)) %>% ## filter top 1% amplified genes
    separate(Gene, into = c("Gene", "code"), sep = " \\(") %>%
    pull(Gene) %>%
    unique()
cn_ampl_noexpr_genes <- read_csv("data/OmicsExpressionProteinCodingGenesTPMLogp1.csv") %>%
    rename(ModelID = colnames(.)[1]) %>%
    pivot_longer(-ModelID, names_to = "Gene", values_to = "TPM") %>%
    group_by(Gene) %>%
    mutate(mean_TPM = mean(TPM, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(mean_TPM < 1) %>% ## filter genes with no expression
    select(Gene) %>%
    distinct() %>%
    separate(Gene, into = c("Gene", "code"), sep = " \\(") %>%
    pull(Gene) %>%
    unique() %>%
    intersect(cn_ampl_genes) ## filter genes with no expression and amplified


# compute recall
get_recall <- function(dfs_corr, dfs_sig, gene_set){
    ## bind dfs
    dfs_corr <- map(dfs_corr, ~.x %>%
                    pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                    rename(Gene = colnames(.)[1]) %>%
                    filter(Gene %in% gene_set)) %>%
                bind_rows(.id = "Algorithm")
    
    res <- full_join(dfs_corr, dfs_sig) %>% 
        group_by(Algorithm, ModelID) %>% 
        mutate(tot=n()) %>% 
        mutate(recall = sum(LFC < Sig_Threshold)/tot) %>% 
        select(Algorithm, ModelID, recall) %>% 
        ungroup() %>% 
        distinct()
    
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
            mutate(across(where(is.numeric), ~replace_na(., 0)))) %>% ## fill na with 0
        set_names(dfs_names)
    
    ## common cell lines
    common_cells <- Reduce(intersect, map(dfs, ~colnames(.)[2:length(colnames(.))]))

    ## select only common cell lines and first columns
    dfs <- map(dfs, ~select(., colnames(.)[1], common_cells))

    ## compute significant threshold at 5% FDR for each algorithm across cell lines
    sigthreshold <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                rename(Gene = colnames(.)[1]) %>%
                split(.$ModelID) %>%
                map(~.x %>% 
                    pull(LFC, name = "Gene") %>% 
                    ccr.ROC_Curve(., ess_genes, noness_genes, display = FALSE, FDRth = 0.05) %>% 
                    .$sigthreshold %>% 
                    as.numeric()) %>%
                bind_rows() %>%
                pivot_longer(everything(.), names_to = "ModelID", values_to = "Sig_Threshold")) %>%
        set_names(dfs_names) %>%
        bind_rows(.id = "Algorithm") %>%
        ### replace NA or Inf with 0
        mutate(Sig_Threshold = replace_na(Sig_Threshold, 0)) %>%
        mutate(Sig_Threshold = replace(Sig_Threshold, is.infinite(Sig_Threshold), 0))
    
    ## for each algorithm and cell line, compute recall of gene sets
    recall_gene_sets <- sigthreshold %>%
        mutate(ess_genes = full_join(., get_recall(dfs, ., ess_genes)) %>%
                pull(recall),
            noness_genes = full_join(., get_recall(dfs, ., noness_genes)) %>%
                pull(recall),
            msigdb_genes = full_join(., get_recall(dfs, ., msigdb_genes)) %>%
                pull(recall),
            cn_ampl_genes = full_join(., get_recall(dfs, ., cn_ampl_genes)) %>%
                pull(recall),
            cn_ampl_noexpr_genes = full_join(., get_recall(dfs, ., cn_ampl_noexpr_genes)) %>%
                pull(recall)) %>%
        select(-Sig_Threshold) %>%
        pivot_longer(-c(Algorithm, ModelID), names_to = "Gene_Set", values_to = "Recall")

    ## compute AUROC for each algorithm across cell lines
    aurocs <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                rename(Gene = colnames(.)[1]) %>%
                split(.$ModelID) %>%
                map(~.x %>% 
                    pull(LFC, name = "Gene") %>% 
                    ccr.ROC_Curve(., ess_genes, noness_genes, display = FALSE) %>% 
                    .$AUC %>% 
                    as.numeric()) %>%
                bind_rows() %>%
                pivot_longer(everything(.), names_to = "ModelID", values_to = "AUROC")) %>%
        set_names(dfs_names) %>%
        bind_rows(.id = "Algorithm") %>%
        ### replace NA or Inf with 0
        mutate(AUROC = replace_na(AUROC, 0)) %>%
        mutate(AUROC = replace(AUROC, is.infinite(AUROC), 0))
    
    ## compute AUPRC for each algorithm across cell lines
    auprcs <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                rename(Gene = colnames(.)[1]) %>%
                split(.$ModelID) %>%
                map(~.x %>% 
                    pull(LFC, name = "Gene") %>% 
                    ccr.PrRc_Curve(., ess_genes, noness_genes, display = FALSE) %>%
                    .$AUC %>% 
                    as.numeric()) %>%
                bind_rows() %>%
                pivot_longer(everything(.), names_to = "ModelID", values_to = "AUPRC")) %>%
        set_names(dfs_names) %>%
        bind_rows(.id = "Algorithm") %>%
        ### replace NA or Inf with 0
        mutate(AUPRC = replace_na(AUPRC, 0)) %>%
        mutate(AUPRC = replace(AUPRC, is.infinite(AUPRC), 0))
    

    ## save results
    p_gene_sets <- ggplot(recall_gene_sets, aes(x = Algorithm, y = Recall, fill = Algorithm)) +
        geom_violin() +
        labs(x = "Method", y = "Recall at 5% FDR", title = "Recall of gene sets") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1) +
        facet_wrap(~Gene_Set, scales = "free_y") +
        theme(strip.background = element_blank(),
            strip.text = element_text(size = 10, face = "bold"))
    ggsave(p_gene_sets, filename = paste0("results/analyses/impact_data_quality/", lib, "_recall_gene_sets.pdf"), width = 10, height = 10, dpi = 300)

    p_aurocs <- ggplot(aurocs, aes(x = Algorithm, y = AUROC, fill = Algorithm)) +
        geom_violin() +
        labs(x = "Method", y = "AUROC", title = "AUROC") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_aurocs, filename = paste0("results/analyses/impact_data_quality/", lib, "_AUROC.pdf"), width = 10, height = 10, dpi = 300)
    
    p_auprcs <- ggplot(auprcs, aes(x = Algorithm, y = AUPRC, fill = Algorithm)) +
        geom_violin() +
        labs(x = "Method", y = "AUPRC", title = "AUPRC") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_auprcs, filename = paste0("results/analyses/impact_data_quality/", lib, "_AUPRC.pdf"), width = 10, height = 10, dpi = 300)
}
