library(forcats)
library(tidyverse)
library(writexl)

algos <- c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK")
libs <- c("Avana", "KY")

# iterate over algorithms and libraries
for (lib in libs){
    ## All genes
    dfs <- readRDS(paste0("results/analyses/cn_correction/", lib, "_cn_abs.rds"))
    dfs$Algorithm <- factor(dfs$Algorithm, levels = algos)

    ## LFC vs CN effect size (all genes)
    dfs_es <- dfs %>%
        group_by(Algorithm, CN_abs) %>%
        mutate(es = mean(LFC)/sd(LFC)) %>%
        ungroup() %>%
        select(Algorithm, es) %>%
        distinct() %>%
        group_by(Algorithm) %>%
        summarise(mean_cn_all = mean(es), sd_cn_all = sd(es)) %>%
        ungroup() %>% 
        mutate(mean_cn_all = signif(mean_cn_all, 3), 
          sd_cn_all = signif(sd_cn_all, 3)) %>% 
        unite(cn_all, c("mean_cn_all", "sd_cn_all"), sep = "±")

    ## Consider only unexpressed genes (TPM < 1)
    dfs_unexpr <- readRDS(paste0("results/analyses/cn_correction/", lib, "_cn_abs_tpm.rds"))
    dfs_unexpr$Algorithm <- factor(dfs_unexpr$Algorithm, levels = algos)
    
    ## LFC vs CN effect size (unexpressed genes)
    dfs_es_unexpr <- dfs_unexpr %>%
        group_by(Algorithm, CN_abs) %>%
        mutate(es = mean(LFC)/sd(LFC)) %>%
        ungroup() %>%
        select(Algorithm, es) %>%
        distinct() %>%
        group_by(Algorithm) %>%
        summarise(mean_cn_unexpr = mean(es), sd_cn_unexpr = sd(es)) %>%
        ungroup() %>%
        mutate(mean_cn_unexpr = signif(mean_cn_unexpr, 3), 
          sd_cn_unexpr = signif(sd_cn_unexpr, 3)) %>% 
        unite(cn_unexpr, c("mean_cn_unexpr", "sd_cn_unexpr"), sep = "±")
    
    ## merge results
    cn_bias <- dfs_es %>%
        left_join(dfs_es_unexpr)
    write_xlsx(cn_bias, paste0("results/tables/table_3_", lib, ".xlsx"))
    


    ## pooled results
    bm_pool <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_pool.rds"))
    bm_pool$Algorithm <- factor(bm_pool$Algorithm, levels = algos)
    bm_pool <- bm_pool %>%
        group_by(Algorithm) %>%
        summarise(mean_prox_pool = mean(est), sd_prox_pool = sd(est)) %>%
        ungroup() %>%
        mutate(mean_prox_pool = signif(mean_prox_pool, 3), 
          sd_prox_pool = signif(sd_prox_pool, 3)) %>%
        unite(prox_pool, c("mean_prox_pool", "sd_prox_pool"), sep = "±")
    
    ## TP53 results
    bm_TP53 <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_TP53.rds"))
    bm_TP53$Algorithm <- factor(bm_TP53$Algorithm, levels = algos)
    bm_TP53 <- bm_TP53 %>%
        group_by(Algorithm, Status) %>%
        summarise(mean_prox_TP53 = mean(est), sd_prox_TP53 = sd(est)) %>%
        ungroup() %>%
        mutate(mean_prox_TP53 = signif(mean_prox_TP53, 3), 
          sd_prox_TP53 = signif(sd_prox_TP53, 3)) %>%
        unite(prox_TP53, c("mean_prox_TP53", "sd_prox_TP53"), sep = "±") %>%
        pivot_wider(names_from = Status, values_from = prox_TP53) %>%
        dplyr::rename(prox_TP53_mut = Mut, prox_TP53_wt = WT)
    
    ## merge results
    prox_bias <- bm_pool %>%
        left_join(bm_TP53)
    write_xlsx(prox_bias, paste0("results/tables/table_4_", lib, ".xlsx"))



    ## AUROC ess vs noness genes
    aurocs <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_AUROC.rds"))
    aurocs$Algorithm <- factor(aurocs$Algorithm, levels = algos)
    aurocs <- aurocs %>%
        group_by(Algorithm) %>%
        summarise(mean_AUROC = mean(AUROC), sd_AUROC = sd(AUROC)) %>%
        ungroup() %>%
        mutate(mean_AUROC = signif(mean_AUROC, 3), 
          sd_AUROC = signif(sd_AUROC, 3)) %>%
        unite(AUROC_ess, c("mean_AUROC", "sd_AUROC"), sep = "±")
    
    ## AUPRC ess vs noness genes
    auprcs <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_AUPRC.rds"))
    auprcs$Algorithm <- factor(auprcs$Algorithm, levels = algos)
    auprcs <- auprcs %>%
        group_by(Algorithm) %>%
        summarise(mean_AUPRC = mean(AUPRC), sd_AUPRC = sd(AUPRC)) %>%
        ungroup() %>%
        mutate(mean_AUPRC = signif(mean_AUPRC, 3), 
          sd_AUPRC = signif(sd_AUPRC, 3)) %>%
        unite(AUPRC_ess, c("mean_AUPRC", "sd_AUPRC"), sep = "±")
    
    ## Gene sets separation
    gene_sep <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_gene_sets_separation.rds"))
    gene_sep$Algorithm <- factor(gene_sep$Algorithm, levels = algos)
    gene_sep <- gene_sep %>%
        group_by(Algorithm) %>%
        summarise(mean_NNMD = mean(Separation), sd_NNMD = sd(Separation)) %>%
        ungroup() %>%
        mutate(mean_NNMD = signif(mean_NNMD, 3), 
          sd_NNMD = signif(sd_NNMD, 3)) %>%
        unite(NNMD, c("mean_NNMD", "sd_NNMD"), sep = "±")
    
    ## Recall of oncogenes
    onco_auc <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_onco_auc.rds"))
    onco_auc$Algorithm <- factor(onco_auc$Algorithm, levels = algos)
    onco_auc <- onco_auc %>%
        dplyr::rename(AUROC_onco = AUROC)
    
    ## merge results
    ess_sum <- aurocs %>%
        left_join(auprcs) %>%
        left_join(gene_sep) %>%
        left_join(onco_auc)
    write_xlsx(ess_sum, paste0("results/tables/table_5_", lib, ".xlsx"))



    ## Recall 5% FDR gene sets
    recall_gene_sets <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_gene_sets.rds")) %>%
        filter(Gene_Set %in% c("ess_genes", "noness_genes", "msigdb_genes")) %>%
        mutate(Gene_Set = ifelse(Gene_Set == "ess_genes", "Common essential genes",
            ifelse(Gene_Set == "noness_genes", "Nonessential genes", "MsigDB genes")))
    recall_gene_sets$Algorithm <- factor(recall_gene_sets$Algorithm, levels = algos)
    recall_gene_sets <- recall_gene_sets %>%
        group_by(Algorithm, Gene_Set) %>%
        summarise(mean_recall = mean(Recall), sd_recall = sd(Recall)) %>%
        ungroup() %>%
        mutate(mean_recall = signif(mean_recall, 3), 
          sd_recall = signif(sd_recall, 3)) %>%
        unite(recall_5, c("mean_recall", "sd_recall"), sep = "±") %>%
        pivot_wider(names_from = Gene_Set, values_from = recall_5)
    
    ## Recall curve amplified genes
    rec_ampl <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl.rds"))
    rec_ampl$Algorithm <- factor(rec_ampl$Algorithm, levels = algos)
    rec_ampl <- rec_ampl %>%
        group_by(Algorithm) %>%
        summarise(mean_rec_ampl = mean(Recall), sd_rec_ampl = sd(Recall)) %>%
        ungroup() %>%
        mutate(mean_rec_ampl = signif(mean_rec_ampl, 3), 
          sd_rec_ampl = signif(sd_rec_ampl, 3)) %>%
        unite(Recall_ampl, c("mean_rec_ampl", "sd_rec_ampl"), sep = "±")
    
    ## Recall curve amplified genes (unexpressed)
    rec_ampl_noexpr <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl_noexpr.rds"))
    rec_ampl_noexpr$Algorithm <- factor(rec_ampl_noexpr$Algorithm, levels = algos)
    rec_ampl_noexpr <- rec_ampl_noexpr %>%
        group_by(Algorithm) %>%
        summarise(mean_rec_ampl_noexpr = mean(Recall), sd_rec_ampl_noexpr = sd(Recall)) %>%
        ungroup() %>%
        mutate(mean_rec_ampl_noexpr = signif(mean_rec_ampl_noexpr, 3), 
          sd_rec_ampl_noexpr = signif(sd_rec_ampl_noexpr, 3)) %>%
        unite(Recall_ampl_noexpr, c("mean_rec_ampl_noexpr", "sd_rec_ampl_noexpr"), sep = "±")
      
    ## merge results
    recall_sum <- recall_gene_sets %>%
        left_join(rec_ampl) %>%
        left_join(rec_ampl_noexpr)
    write_xlsx(recall_sum, paste0("results/tables/table_6_", lib, ".xlsx"))



    ## Nº significant biomarkers (all CFEs on strongly selective dependencies)
    sig_biomarkers_ssd <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_ssd.rds"))
    sig_biomarkers_ssd$Algorithm <- factor(sig_biomarkers_ssd$Algorithm, levels = algos)
    sig_biomarkers_ssd <- sig_biomarkers_ssd %>%
        dplyr::rename(n_sig_biomark_ssd = n_sig_biomark)

    ## Nº significant biomarkers (gain-of-function CFEs within oncogene)
    sig_biomarkers_onco <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_onco.rds"))
    sig_biomarkers_onco$Algorithm <- factor(sig_biomarkers_onco$Algorithm, levels = algos)
    sig_biomarkers_onco <- sig_biomarkers_onco %>%
        dplyr::rename(n_sig_biomark_onco = n_sig_biomark)

    ## merge results
    biomarkers_res <- sig_biomarkers_ssd %>%
        left_join(sig_biomarkers_onco)
    write_xlsx(biomarkers_res, paste0("results/tables/table_7_", lib, ".xlsx"))
}
