library(ggradar)
library(extrafont)
library(forcats)
library(patchwork)
library(RColorBrewer)
library(tidyverse)

# provide datapoints for circle
circleFun <- function(center = c(0,0), rx = 1, ry = 1, npoints = 1000){
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + rx * cos(tt)
    yy <- center[2] + ry * sin(tt)
    return(data.frame(x = xx, y = yy))
}

font_import(paths = "/group/iorio/Alessandro/CN_benchmark/arial", prompt = FALSE)

algos <- c("CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK")
libs <- c("Avana", "KY")
cols <- brewer.pal(n = 7, name = "Dark2")

# iterate over algorithms and libraries
for (lib in libs){
    ## All genes
    dfs <- readRDS(paste0("results/analyses/cn_correction/", lib, "_cn_abs.rds"))
    dfs$Algorithm <- factor(dfs$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    ## LFC vs CN effect size (all genes)
    dfs_es <- dfs %>%
        group_by(Algorithm, CN_abs) %>%
        mutate(es = mean(LFC)/sd(LFC)) %>%
        ungroup() %>%
        select(Algorithm, es) %>%
        group_by(Algorithm) %>%
        summarise(es = median(abs(es))) %>%
        ungroup() %>%
        dplyr::rename(cn_bias_all = es)

    ## Consider only unexpressed genes (TPM < 1)
    dfs_unexpr <- readRDS(paste0("results/analyses/cn_correction/", lib, "_cn_abs_tpm.rds"))
    dfs_unexpr$Algorithm <- factor(dfs_unexpr$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    
    ## LFC vs CN effect size (unexpressed genes)
    dfs_es_unexpr <- dfs_unexpr %>%
        group_by(Algorithm, CN_abs) %>%
        mutate(es = mean(LFC)/sd(LFC)) %>%
        ungroup() %>%
        select(Algorithm, es) %>%
        group_by(Algorithm) %>%
        summarise(es = median(abs(es))) %>%
        ungroup() %>%
        dplyr::rename(cn_bias_unexpr = es)
    
    ## pooled results
    bm_pool <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_pool.rds"))
    bm_pool$Algorithm <- factor(bm_pool$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    bm_pool <- bm_pool %>%
        group_by(Algorithm) %>%
        summarise(est = median(est)) %>%
        ungroup() %>%
        dplyr::rename(proximity_bias = est)
    
    ## merge results
    res <- dfs_es %>%
        left_join(dfs_es_unexpr) %>%
        left_join(bm_pool) %>%
        filter(Algorithm != "Uncorrected") %>%
        mutate(Type = c("Unsupervised", "Supervised", "Supervised", 
            "Supervised", "Unsupervised", "Unsupervised", "Supervised")) %>%
        mutate(MoA = c("Single-screen", "Multi-screen", "Single-screen", 
            "Single-screen", "Single-screen", "Single-screen", "Multi-screen"))
    
    ## plot
    radius_1 <- circleFun(c(0,0.5), 0.1, 0.05)
    radius_2 <- circleFun(c(0,0.5), 0.2, 0.1)
    radius_3 <- circleFun(c(0,0.5), 0.3, 0.15)
    radius_4 <- circleFun(c(0,0.5), 0.4, 0.2)
    radius_5 <- circleFun(c(0,0.5), 0.5, 0.25)
    radius_6 <- circleFun(c(0,0.5), 0.6, 0.3)
    radius_7 <- circleFun(c(0,0.5), 0.7, 0.35)
    radius_8 <- circleFun(c(0,0.5), 0.8, 0.4)
    radius_9 <- circleFun(c(0,0.5), 0.9, 0.45)
    radius_10 <- circleFun(c(0,0.5), 1, 0.5)

    ### all genes
    p_es <- ggplot(res, aes(x = cn_bias_all, y = proximity_bias)) +
        geom_point(aes(x = cn_bias_all, y = proximity_bias,
            shape = Type, color = MoA),
            size = 9) +
        geom_text(
            label = res$Algorithm, 
            nudge_x = 0, nudge_y = 0.012, 
            size = 8,
            check_overlap = FALSE) +
        geom_path(data = radius_1, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_2, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_3, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_4, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_5, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_6, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        labs(x = "Median correction \nfor CN bias (all genes)", 
            y = "",
            title = "All genes") +
        xlim(c(0, 0.35)) +
        ylim(c(0.5, 0.75)) +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.title = element_text(size = 35, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 18),
            text = element_text(family = "Arial"),
            aspect.ratio = 1,
            plot.title = element_text(size = 32, hjust = 0.5, face = "bold"),
            plot.margin = grid::unit(c(1,1,1,1), "cm")) +
        scale_color_manual(values = c("red", "blue"))

    ### unexpressed genes
    p_es_unexpr <- ggplot(res, aes(x = cn_bias_unexpr, y = proximity_bias)) +
        geom_point(aes(x = cn_bias_unexpr, y = proximity_bias,
            shape = Type, color = MoA),
            size = 9) +
        geom_text(
            label = res$Algorithm, 
            nudge_x = 0, nudge_y = 0.012, 
            size = 9,
            check_overlap = FALSE) +
        geom_path(data = radius_1, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_2, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_3, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_4, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_5, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_6, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_7, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_8, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_9, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        geom_path(data = radius_10, aes(x = x, y = y), linetype = "dashed", color = alpha("black", 0.5)) +
        labs(x = "Median correction \nfor CN bias (unexpressed genes)", 
            y = "Median correction \nfor proximity bias",
            title = "Unexpressed genes") +
        xlim(c(0, 0.9)) +
        ylim(c(0.5, 0.75)) +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.title = element_text(size = 35, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            text = element_text(family = "Arial"),
            aspect.ratio = 1,
            plot.title = element_text(size = 32, hjust = 0.5, face = "bold"),
            plot.margin = grid::unit(c(1,1,1,1), "cm")) +
        scale_color_manual(values = c("red", "blue"))
    

    ## AUPRC ess vs noness genes
    auprcs <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_AUPRC.rds"))
    auprcs$Algorithm <- factor(auprcs$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    auprcs <- auprcs %>%
        group_by(Algorithm) %>%
        summarise(AUPRC = median(AUPRC)) %>%
        ungroup() %>%
        dplyr::rename(AUPRC_ess = AUPRC)
    
    ## Gene sets separation
    gene_sep <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_gene_sets_separation.rds"))
    gene_sep$Algorithm <- factor(gene_sep$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    gene_sep <- gene_sep %>%
        group_by(Algorithm) %>%
        summarise(Separation = median(Separation)) %>%
        ungroup() %>%
        dplyr::rename(NNMD = Separation)
    
    ## Recall of oncogenes
    onco_auc <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_onco_auc.rds"))
    onco_auc$Algorithm <- factor(onco_auc$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    onco_auc <- onco_auc %>%
        dplyr::rename(AUROC_onco = AUROC)

    ## Recall curve amplified genes
    rec_ampl <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl.rds"))
    rec_ampl$Algorithm <- factor(rec_ampl$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    rec_ampl <- rec_ampl %>%
        group_by(Algorithm) %>%
        summarise(Recall = median(Recall)) %>%
        ungroup() %>%
        dplyr::rename(Recall_ampl = Recall)
    
    ## Recall curve amplified genes (unexpressed)
    rec_ampl_noexpr <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl_noexpr.rds"))
    rec_ampl_noexpr$Algorithm <- factor(rec_ampl_noexpr$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    rec_ampl_noexpr <- rec_ampl_noexpr %>%
        group_by(Algorithm) %>%
        summarise(Recall = median(Recall)) %>%
        ungroup() %>%
        dplyr::rename(Recall_ampl_noexpr = Recall)

    ## Nº significant biomarkers (all CFEs on strongly selective dependencies)
    sig_biomarkers_ssd <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_ssd.rds"))
    sig_biomarkers_ssd$Algorithm <- factor(sig_biomarkers_ssd$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    sig_biomarkers_ssd <- sig_biomarkers_ssd %>%
        dplyr::rename(n_sig_biomark_ssd = n_sig_biomark)

    ## Nº significant biomarkers (gain-of-function CFEs within oncogene)
    sig_biomarkers_onco <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_onco.rds"))
    sig_biomarkers_onco$Algorithm <- factor(sig_biomarkers_onco$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    sig_biomarkers_onco <- sig_biomarkers_onco %>%
        dplyr::rename(n_sig_biomark_onco = n_sig_biomark)

    ## merge results
    params <- auprcs %>%
        left_join(onco_auc) %>%
        left_join(rec_ampl) %>%
        left_join(rec_ampl_noexpr) %>%
        left_join(sig_biomarkers_ssd) %>%
        left_join(sig_biomarkers_onco) %>%
        left_join(gene_sep) %>%
        mutate(Recall_ampl = 1 - abs(0.5 - Recall_ampl)*2,
            Recall_ampl_noexpr = 1 - abs(0.5 - Recall_ampl_noexpr)*2,
            n_sig_biomark_ssd = n_sig_biomark_ssd/max(n_sig_biomark_ssd),
            n_sig_biomark_onco = n_sig_biomark_onco/max(n_sig_biomark_onco),
            NNMD = NNMD/min(NNMD)) %>%
        filter(Algorithm != "Uncorrected")
    
    ## plot
    radar <- ggradar(params,
        group.point.size = 0,
        group.line.width = 2,
        values.radar = c("", "", ""),
        fill = FALSE,
        fill.alpha = 0.7,
        gridline.min.colour = "black",
        gridline.mid.colour = "black",
        gridline.max.colour = "black",
        background.circle.colour = "white",
        legend.position = "right",
        font.radar = "Arial",
        axis.labels = c("AUPRC common\n essential genes", 
            "AUROC\noncogenes", "Recall\n(amplified\ngenes)", 
            "Recall (amplified\nunexpressed genes)", 
            "Nº significant\nbiomarkers\n(all CFEs)", 
            "Nº significant\nbiomarkers\n(GOF CFEs)",
            "NNMD"),
        axis.label.size = 12) +
        scale_color_manual(values = cols) +
        theme(legend.text = element_text(size = 18))

    ## Assemble panel
    panel <- (p_es_unexpr | p_es ) / radar +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag.position = c(0, 1),
            plot.tag = element_text(size = 40, face = "bold", family = "Arial"))
    ggsave(panel, filename = paste0("results/panels/summary_panel_", lib, ".pdf"), width = 40, height = 40)
}
