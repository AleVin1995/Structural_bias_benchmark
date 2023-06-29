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

msigdb_genes <- c(
        EssGenes.DNA_REPLICATION_cons,
        EssGenes.HISTONES,
        EssGenes.KEGG_rna_polymerase,
        EssGenes.PROTEASOME_cons,
        EssGenes.SPLICEOSOME_cons,
        EssGenes.ribosomalProteins
        ) %>%
    unique()

cn_ratio <- read_csv("data/OmicsCNGene.csv") %>%
    dplyr::rename(ModelID = colnames(.)[1]) %>%
    pivot_longer(-ModelID, names_to = "Gene", values_to = "CN_ratio") %>%
    drop_na()

cn_ratio_tpm <- cn_ratio %>%
    inner_join(read_csv("data/OmicsExpressionProteinCodingGenesTPMLogp1.csv") %>%
        dplyr::rename(ModelID = colnames(.)[1]) %>%
        pivot_longer(-ModelID, names_to = "Gene", values_to = "TPM")) %>%
    drop_na()

cn_ampl_genes <- cn_ratio %>%
    group_split(ModelID) %>%
    map(~.x %>%
        arrange(desc(CN_ratio)) %>%
        ## filter top 1% amplified genes per cell line
        dplyr::slice(1:round(nrow(.)*0.01))) %>%
    bind_rows()
    
cn_ampl_noexpr_genes <- cn_ratio_tpm %>%
    inner_join(cn_ampl_genes) %>%
    group_by(Gene) %>%
    mutate(mean_TPM = mean(TPM, na.rm = TRUE)) %>%
    ungroup() %>%
    ## filter genes with no expression
    filter(mean_TPM < 1)


# compute recall
get_recall <- function(dfs_corr, dfs_sig, gene_info){
    if (is.vector(gene_info)){
        gene_set <- gene_info

        ## bind dfs and filter genes
        dfs_corr <- map(dfs_corr, ~.x %>%
                        pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                        dplyr::rename(Gene = colnames(.)[1]) %>%
                        filter(Gene %in% gene_set)) %>%
                    bind_rows(.id = "Algorithm")
    } else {
        ## bind dfs and filter genes
        dfs_corr <- map(dfs_corr, ~.x %>%
                        pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                        dplyr::rename(Gene = colnames(.)[1]) %>%
                        inner_join(gene_info, by = "Gene")) %>%
                    bind_rows(.id = "Algorithm")
    }

    res <- full_join(dfs_corr, dfs_sig) %>%
            group_by(Algorithm, ModelID) %>%
            mutate(tot=n()) %>%
            mutate(recall = sum(LFC < Sig_Threshold)/tot) %>%
            select(Algorithm, ModelID, recall) %>%
            ungroup() %>%
            distinct()
    
    return(res)
}


# run ROC curve
run_Curve <- function(
    FCsprofile,
    positives,
    negatives,
    display = TRUE,
    FDRth = NULL,
    expName = NULL,
    type = NULL
) {
    ## turn positive genes into a vector if it is a dataframe
    if (!is.vector(positives)){
        positives <- positives %>% 
            inner_join(FCsprofile, by = "Gene") %>% 
            pull(Gene) %>%
            unique()
    }

    FCsprofile <- pull(.data = FCsprofile, var = LFC, name = "Gene")

    if (is.null(type) | type == "ROC"){
        res <- ccr.ROC_Curve(FCsprofile, 
            positives, 
            negatives, 
            display = display, 
            FDRth = FDRth, 
            expName = expName)
    } else if (type == "PR"){
        res <- ccr.PrRc_Curve(FCsprofile, 
            positives, 
            negatives, 
            display = display, 
            FDRth = FDRth, 
            expName = expName)
    }
    
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
            dplyr::dplyr::rename(Gene = colnames(.)[1])) %>% ## fill na with 0
        set_names(dfs_names)
    
    ## common cell lines/genes
    common_cells <- Reduce(intersect, map(dfs, ~colnames(.)[2:length(colnames(.))]))

    common_genes <- Reduce(intersect, map(dfs, ~.$Gene))

    ## select only common cell lines and first columns
    dfs <- map(dfs, ~select(., colnames(.)[1], all_of(common_cells))) %>%
        map(~filter(., Gene %in% common_genes))

    ## compute significant threshold at 5% FDR for each algorithm across cell lines
    sigthreshold <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                dplyr::rename(Gene = colnames(.)[1]) %>%
                group_split(ModelID) %>%
                map(~.x %>% 
                    run_Curve(., 
                        ess_genes, 
                        noness_genes, 
                        display = FALSE, 
                        FDRth = 0.05) %>% 
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

    ## compute AUROC for each algorithm across cell lines (ess vs noness genes)
    aurocs <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                dplyr::rename(Gene = colnames(.)[1]) %>%
                group_split(ModelID) %>%
                map(~.x %>%
                    run_Curve(.,
                        ess_genes, 
                        noness_genes, 
                        display = FALSE) %>% 
                    .$AUC %>% 
                    as.numeric()) %>%
                bind_rows() %>%
                pivot_longer(everything(.), names_to = "ModelID", values_to = "AUROC")) %>%
        set_names(dfs_names) %>%
        bind_rows(.id = "Algorithm") %>%
        ### replace NA or Inf with 0
        mutate(AUROC = replace_na(AUROC, 0)) %>%
        mutate(AUROC = replace(AUROC, is.infinite(AUROC), 0))
    
    ## compute AUROC for each algorithm across cell lines (ampl vs noness genes)
    aurocs_ampl <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                dplyr::rename(Gene = colnames(.)[1]) %>%
                group_split(ModelID) %>%
                map(~.x %>%
                    run_Curve(., 
                        cn_ampl_genes, 
                        noness_genes, 
                        display = FALSE) %>% 
                    .$AUC %>% 
                    as.numeric()) %>%
                bind_rows() %>%
                pivot_longer(everything(.), names_to = "ModelID", values_to = "AUROC")) %>%
        set_names(dfs_names) %>%
        bind_rows(.id = "Algorithm") %>%
        ### replace NA or Inf with 0
        mutate(AUROC = replace_na(AUROC, 0)) %>%
        mutate(AUROC = replace(AUROC, is.infinite(AUROC), 0))
    
    ## compute AUROC for each algorithm across cell lines (ampl noexpr vs noness genes)
    aurocs_ampl_noexpr <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                dplyr::rename(Gene = colnames(.)[1]) %>%
                group_split(ModelID) %>%
                map(~.x %>%
                    run_Curve(., 
                        cn_ampl_noexpr_genes, 
                        noness_genes, 
                        display = FALSE) %>% 
                    .$AUC %>% 
                    as.numeric()) %>%
                bind_rows() %>%
                pivot_longer(everything(.), names_to = "ModelID", values_to = "AUROC")) %>%
        set_names(dfs_names) %>%
        bind_rows(.id = "Algorithm") %>%
        ### replace NA or Inf with 0
        mutate(AUROC = replace_na(AUROC, 0)) %>%
        mutate(AUROC = replace(AUROC, is.infinite(AUROC), 0))
    
    ## compute AUPRC for each algorithm across cell lines (ess vs noness genes)
    auprcs <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                dplyr::rename(Gene = colnames(.)[1]) %>%
                group_split(ModelID) %>%
                map(~.x %>% 
                    run_Curve(., 
                        ess_genes, 
                        noness_genes, 
                        display = FALSE) %>%
                    .$AUC %>% 
                    as.numeric()) %>%
                bind_rows() %>%
                pivot_longer(everything(.), names_to = "ModelID", values_to = "AUPRC")) %>%
        set_names(dfs_names) %>%
        bind_rows(.id = "Algorithm") %>%
        ### replace NA or Inf with 0
        mutate(AUPRC = replace_na(AUPRC, 0)) %>%
        mutate(AUPRC = replace(AUPRC, is.infinite(AUPRC), 0))
    
    ## compute AUPRC for each algorithm across cell lines (ampl vs noness genes)
    auprcs_ampl <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                dplyr::rename(Gene = colnames(.)[1]) %>%
                group_split(ModelID) %>%
                map(~.x %>% 
                    run_Curve(., 
                        cn_ampl_genes,
                        noness_genes, 
                        display = FALSE) %>%
                    .$AUC %>% 
                    as.numeric()) %>%
                bind_rows() %>%
                pivot_longer(everything(.), names_to = "ModelID", values_to = "AUPRC")) %>%
        set_names(dfs_names) %>%
        bind_rows(.id = "Algorithm") %>%
        ### replace NA or Inf with 0
        mutate(AUPRC = replace_na(AUPRC, 0)) %>%
        mutate(AUPRC = replace(AUPRC, is.infinite(AUPRC), 0))
    
    ## compute AUROC for each algorithm across cell lines (ampl noexpr vs noness genes)
    auprcs_ampl_noexpr <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                dplyr::rename(Gene = colnames(.)[1]) %>%
                group_split(ModelID) %>%
                map(~.x %>%
                    run_Curve(., 
                        cn_ampl_noexpr_genes, 
                        noness_genes, 
                        display = FALSE) %>% 
                    .$AUC %>% 
                    as.numeric()) %>%
                bind_rows() %>%
                pivot_longer(everything(.), names_to = "ModelID", values_to = "AUROC")) %>%
        set_names(dfs_names) %>%
        bind_rows(.id = "Algorithm") %>%
        ### replace NA or Inf with 0
        mutate(AUPRC = replace_na(AUPRC, 0)) %>%
        mutate(AUPRC = replace(AUPRC, is.infinite(AUPRC), 0))
    

    ## save results
    saveRDS(recall_gene_sets, paste0("results/analyses/impact_data_quality/", lib, "_recall_gene_sets.rds"))
    saveRDS(aurocs, paste0("results/analyses/impact_data_quality/", lib, "_AUROC.rds"))
    saveRDS(aurocs_ampl, paste0("results/analyses/impact_data_quality/", lib, "_AUROC_ampl.rds"))
    saveRDS(aurocs_ampl_noexpr, paste0("results/analyses/impact_data_quality/", lib, "_AUROC_ampl_noexpr.rds"))
    saveRDS(auprcs, paste0("results/analyses/impact_data_quality/", lib, "_AUPRC.rds"))
    saveRDS(auprcs_ampl, paste0("results/analyses/impact_data_quality/", lib, "_AUPRC_ampl.rds"))
    saveRDS(auprcs_ampl_noexpr, paste0("results/analyses/impact_data_quality/", lib, "_AUPRC_ampl_noexpr.rds"))
}
