library(CRISPRcleanR)
library(tidyverse)


# load essential genes from CRISPRcleanR (MSigDB)
load("data/gene_sets/histones.rdata")
load("data/gene_sets/Kegg.DNArep.rdata")
load("data/gene_sets/Kegg.Proteasome.rdata")
load("data/gene_sets/Kegg.Ribosome.rdata") 
load("data/gene_sets/Kegg.RNApoly.rdata")   
load("data/gene_sets/Kegg.Spliceosome.rdata")


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
        Kegg.DNArep,
        Kegg.Proteasome,
        Kegg.Ribosome,
        Kegg.RNApoly,
        Kegg.Spliceosome,
        Histones
        ) %>%
    unique()

cn_ratio <- read_csv("data/OmicsCNGene.csv") %>%
    dplyr::rename(ModelID = colnames(.)[1]) %>%
    rename_with(~sub(" \\(.*$", "", .x)) %>%
    pivot_longer(-ModelID, names_to = "Gene", values_to = "CN_ratio") %>%
    drop_na()

cn_ratio_tpm <- cn_ratio %>%
    inner_join(read_csv("data/OmicsExpressionProteinCodingGenesTPMLogp1.csv") %>%
        dplyr::rename(ModelID = colnames(.)[1]) %>%
        rename_with(~sub(" \\(.*$", "", .x)) %>%
        pivot_longer(-ModelID, names_to = "Gene", values_to = "TPM")) %>%
    drop_na()

cn_ampl_genes <- cn_ratio %>%
    split(.$ModelID) %>%
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


# compute recall at significant threshold
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
                        dplyr::rename(Gene = colnames(.)[1])) %>%
                    bind_rows(.id = "Algorithm") %>%
                    inner_join(gene_info, by = c("Gene", "ModelID"))
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


# compute recall curve for a given gene set
get_recall_curve <- function(
    FCsprofile,
    geneSet
) {
    ## turn positive genes into a vector if it is a dataframe
    if (!is.vector(geneSet)){
        geneSet <- geneSet %>% 
            inner_join(FCsprofile) %>% 
            pull(Gene) %>%
            unique()
    }
    
    FCsprofile <- pull(.data = FCsprofile, var = LFC, name = "Gene")
    FCsprofile <- sort(FCsprofile)

    geneSet <- intersect(geneSet, names(FCsprofile))

    predictions <- rep(0, length(FCsprofile))
    predictions[which(names(FCsprofile) %in% geneSet)] <- 1
    recall <- cumsum(predictions) / length(geneSet)

    AUC <- trapz(seq_along(predictions) / length(predictions), recall)

    return(AUC)
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
            inner_join(FCsprofile) %>% 
            pull(Gene) %>%
            unique()
    }

    FCsprofile <- pull(.data = FCsprofile, var = LFC, name = "Gene")

    if (is.null(type) || type == "ROC"){
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
            dplyr::rename(Gene = colnames(.)[1])) %>% ## fill na with 0
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
                split(.$ModelID) %>%
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
                split(.$ModelID) %>%
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
    
    ## compute AUPRC for each algorithm across cell lines (ess vs noness genes)
    auprcs <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                dplyr::rename(Gene = colnames(.)[1]) %>%
                split(.$ModelID) %>%
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
    
    ## compute recall curve for each algorithm across cell lines (ampl genes)
    rec_ampl <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                dplyr::rename(Gene = colnames(.)[1]) %>%
                split(.$ModelID) %>%
                map(~.x %>%
                    get_recall_curve(., 
                        cn_ampl_genes) %>% 
                    as.numeric()) %>%
                bind_rows() %>%
                pivot_longer(everything(.), names_to = "ModelID", values_to = "Recall")) %>%
        set_names(dfs_names) %>%
        bind_rows(.id = "Algorithm") %>%
        ### replace NA or Inf with 0
        mutate(Recall = replace_na(Recall, 0)) %>%
        mutate(Recall = replace(Recall, is.infinite(Recall), 0))
    
    ## compute recall curve for each algorithm across cell lines (ampl noexpr genes)
    rec_ampl_noexpr <- map(dfs, ~.x %>%
                pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
                dplyr::rename(Gene = colnames(.)[1]) %>%
                split(.$ModelID) %>%
                map(~.x %>%
                    get_recall_curve(.,
                        cn_ampl_noexpr_genes) %>%
                    as.numeric()) %>%
                bind_rows() %>%
                pivot_longer(everything(.), names_to = "ModelID", values_to = "Recall")) %>%
        set_names(dfs_names) %>%
        bind_rows(.id = "Algorithm") %>%
        ### replace NA or Inf with 0
        mutate(Recall = replace_na(Recall, 0)) %>%
        mutate(Recall = replace(Recall, is.infinite(Recall), 0))
    

    ## save results
    saveRDS(recall_gene_sets, paste0("results/analyses/impact_data_quality/", lib, "_recall_gene_sets.rds"))
    saveRDS(aurocs, paste0("results/analyses/impact_data_quality/", lib, "_AUROC.rds"))
    saveRDS(auprcs, paste0("results/analyses/impact_data_quality/", lib, "_AUPRC.rds"))
    saveRDS(rec_ampl, paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl.rds"))
    saveRDS(rec_ampl_noexpr, paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl_noexpr.rds"))
}
