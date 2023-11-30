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

font_import(paths = "arial", prompt = FALSE)

algos <- c("CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK")
libs <- c("Avana", "KY")
cols <- brewer.pal(n = 7, name = "Dark2")

# iterate over algorithms and libraries
for (lib in libs){    
    ## Nº significant biomarkers (all CFEs on strongly selective dependencies)
    sig_biomarkers_ssd <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_ssd.rds"))
    sig_biomarkers_ssd$Algorithm <- factor(sig_biomarkers_ssd$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    p_sig_biomark_ssd <- ggplot(sig_biomarkers_ssd, aes(x = Algorithm, y = n_sig_biomark_ssd, fill = Algorithm)) +
        geom_bar(stat = "identity") +
        labs(x = "", y = "Nº significant associations") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 25, color = 'black'),
            axis.title = element_text(size = 30),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1,
            plot.margin = grid::unit(c(2,2,2,2), "cm"),
            legend.position = "none") +
        scale_fill_manual(values = c("#B3B3B3", cols))

    ## pooled results
    bm_pool <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_pool.rds"))
    bm_pool$Algorithm <- factor(bm_pool$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    bm_pool <- bm_pool %>%
        group_by(Algorithm) %>%
        summarise(est = median(est)) %>%
        ungroup() %>%
        dplyr::rename(proximity_bias = est)
    
    ## Recall curve amplified genes (unexpressed)
    rec_ampl_noexpr <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl_noexpr.rds"))
    rec_ampl_noexpr$Algorithm <- factor(rec_ampl_noexpr$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    rec_ampl_noexpr <- rec_ampl_noexpr %>%
        group_by(Algorithm) %>%
        summarise(Recall = median(Recall)) %>%
        ungroup() %>%
        dplyr::rename(Recall_ampl_noexpr = Recall)
    
    ## merge results
    res <- rec_ampl_noexpr %>%
        left_join(bm_pool) %>%
        filter(Algorithm != "Uncorrected") %>%
        mutate(Type = c("Unsupervised", "Supervised", "Supervised", 
            "Supervised", "Unsupervised", "Unsupervised", "Supervised")) %>%
        mutate(MoA = c("Single-screen", "Multi-screen", "Single-screen", 
            "Single-screen", "Single-screen", "Single-screen", "Multi-screen"))
    colnames(res)[2] <- "cn_bias_unexpr"
    
    ## plot
    radius_1 <- circleFun(c(0.5,0.5), 0.1, 0.05)
    radius_2 <- circleFun(c(0.5,0.5), 0.2, 0.1)
    radius_3 <- circleFun(c(0.5,0.5), 0.3, 0.15)
    radius_4 <- circleFun(c(0.5,0.5), 0.4, 0.2)
    radius_5 <- circleFun(c(0.5,0.5), 0.5, 0.25)
    radius_6 <- circleFun(c(0.5,0.5), 0.6, 0.3)
    radius_7 <- circleFun(c(0.5,0.5), 0.7, 0.35)
    radius_8 <- circleFun(c(0.5,0.5), 0.8, 0.4)
    radius_9 <- circleFun(c(0.5,0.5), 0.9, 0.45)
    radius_10 <- circleFun(c(0.5,0.5), 1, 0.5)

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
        geom_vline(xintercept = 0.5, linetype = 'longdash') +
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
            y = "Median correction \nfor proximity bias") +
        xlim(c(0.35, 0.65)) +
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

    ## Nº significant biomarkers (all CFEs on strongly selective dependencies)
    sig_biomarkers_ssd <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_ssd.rds"))
    sig_biomarkers_ssd$Algorithm <- factor(sig_biomarkers_ssd$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    sig_biomarkers_ssd <- sig_biomarkers_ssd %>%
        dplyr::rename(n_sig_biomark_ssd = n_sig_biomark)

    ## merge results
    params <- auprcs %>%
        left_join(onco_auc) %>%
        left_join(sig_biomarkers_ssd) %>%
        left_join(gene_sep) %>%
        mutate(n_sig_biomark_ssd = n_sig_biomark_ssd/max(n_sig_biomark_ssd),
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
            "AUROC\noncogenes", 
            "Nº significant\nbiomarkers\n(all CFEs)",
            "NNMD"),
        axis.label.size = 12) +
        scale_color_manual(values = cols) +
        theme(legend.text = element_text(size = 18))
    
    ## Area
    params <- params %>% 
        mutate(Area = 0.5**2 * (AUPRC_ess * AUROC_onco + 
            AUROC_onco * n_sig_biomark_ssd + 
            n_sig_biomark_ssd * NNMD +
            NNMD * AUPRC_ess))

    area <- ggplot(params, aes(x = Algorithm, y = Area, fill = Algorithm)) +
        geom_bar(stat = "identity") +
        labs(x = "", y = "Area (radar plot)") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 25, color = 'black'),
            axis.title = element_text(size = 30),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1,
            plot.margin = grid::unit(c(2,2,2,2), "cm"),
            legend.position = "none") +
        scale_fill_manual(values = cols) +
        coord_cartesian(ylim = c(0.5, 0.8))

    ## Assemble panel
    panel <- (p_sig_biomark_ssd | p_es_unexpr) / (radar | area) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag.position = c(0, 1),
            plot.tag = element_text(size = 40, face = "bold", family = "Arial"))
    ggsave(panel, filename = paste0("results/panels/summary_panel_", lib, ".pdf"), width = 40, height = 40)
}
