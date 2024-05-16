library(extrafont)
library(forcats)
library(patchwork)
library(RColorBrewer)
library(tidyverse)

font_import(paths = "arial", prompt = FALSE)

source("src/figures/utils.r")

libs <- c("Avana", "KY")
cols <- c("#B3B3B3", brewer.pal(n = 7, name = "Dark2"))
cols <- c(cols[1:3], "#EE8208", cols[4:8])

# iterate over algorithms and libraries
for (lib in libs){
    ## All genes
    # dfs <- readRDS(paste0("results/analyses/cn_correction/", lib, "_cn_abs.rds"))
    # dfs$Algorithm <- factor(dfs$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    # p_cn_abs <- ggplot(dfs, aes(x = as.factor(CN_abs), y = LFC, fill = Algorithm)) +
    #     stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + 
    #     labs(x = "Absolute copy number", y = "LFC") +
    #     theme_bw() +
    #     theme(
    #         axis.text = element_text(size = 25, color = 'black'),
    #         axis.text.x = element_text(angle = 45, hjust = 1),
    #         axis.ticks.length = unit(0.5, "cm"),
    #         axis.title = element_text(size = 35, color = 'black'),
    #         strip.background = element_blank(), 
    #         strip.placement = "outside",
    #         strip.text = element_text(size = 35, color = 'black'),
    #         panel.border = element_blank(),
    #         panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(),
    #         panel.background = element_blank(),
    #         legend.position = "none",
    #         text = element_text(family = "Arial")) +
    #     geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    #     facet_wrap(~Algorithm, scales = "free_x", ncol = 3) +
    #     scale_fill_manual(values = cols)

    ## Recall curve amplified genes
    # rec_ampl <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl.rds"))
    # rec_ampl$Algorithm <- factor(rec_ampl$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    # p_rec_ampl <- ggplot(rec_ampl, aes(x = Algorithm, y = Recall, fill = Algorithm)) +
    #     geom_boxplot() +
    #     geom_hline(yintercept = 0.5, linetype = "dashed") +
    #     labs(x = "", y = "AURC") +
    #     theme_bw() +
    #     theme(
    #         panel.border = element_blank(),
    #         panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(),
    #         axis.text = element_text(size = 25, color = 'black'),
    #         axis.ticks.length = unit(0.5, "cm"),
    #         axis.title = element_text(size = 30),
    #         axis.text.x = element_text(angle = 45, hjust = 1),
    #         plot.title = element_text(size = 32, hjust = 0.5, face = "bold"),
    #         legend.position = "none") +
    #     scale_fill_manual(values = cols)


    ## Consider only unexpressed genes (TPM < 1)
    dfs_unexpr <- readRDS(paste0("results/analyses/cn_correction/", lib, "_cn_abs_tpm.rds"))
    dfs_unexpr$Algorithm <- factor(dfs_unexpr$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    
    p_cn_abs_unexpr <- ggplot(dfs_unexpr, aes(x = as.factor(CN_abs), y = LFC, fill = Algorithm)) +
        stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + 
        labs(x = "Unexpressed genes (TPM <1) absolute copy number", y = "Gene depletion log fold-change") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
            axis.ticks.length = unit(0.5, "cm"),
            axis.title = element_text(size = 35, color = 'black', hjust = 0.5),
            strip.background = element_blank(), 
            strip.placement = "outside",
            strip.text = element_text(size = 35, color = 'black'),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            text = element_text(family = "Arial")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        facet_wrap(~Algorithm, scales = "free_x", ncol = 3) +
        scale_fill_manual(values = cols) +
        scale_x_discrete(breaks = seq(0, 15, by = 2))

    ## Recall curve amplified genes (unexpressed)
    # rec_ampl_noexpr <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl_noexpr.rds"))
    # rec_ampl_noexpr$Algorithm <- factor(rec_ampl_noexpr$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    
    # p_rec_ampl_noexpr <- ggplot(rec_ampl_noexpr, aes(x = Algorithm, y = Recall, fill = Algorithm)) +
    #     geom_boxplot() +
    #     geom_hline(yintercept = 0.5, linetype = "dashed") +
    #     labs(x = "", y = "AURC") +
    #     theme_bw() +
    #     theme(
    #         panel.border = element_blank(),
    #         panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(),
    #         axis.text = element_text(size = 25, color = 'black'),
    #         axis.ticks.length = unit(0.5, "cm"),
    #         axis.title = element_text(size = 30),
    #         axis.text.x = element_text(angle = 45, hjust = 1),
    #         legend.position = "none") +
    #     scale_fill_manual(values = cols)
    
    ## Recall curve amplified genes (unexpressed + background unexpressed)
    rec_ampl_noexpr_bg_noexpr <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl_noexpr_bg_noexpr.rds"))
    rec_ampl_noexpr_bg_noexpr$Algorithm <- factor(rec_ampl_noexpr_bg_noexpr$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    
    p_rec_ampl_noexpr_bg_noexpr <- ggplot(rec_ampl_noexpr_bg_noexpr, aes(x = Algorithm, y = Recall, fill = Algorithm)) +
        geom_boxplot() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "AURC") +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 25, color = 'black'),
            axis.title = element_text(size = 30),
            axis.ticks.length = unit(0.5, "cm"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
        scale_fill_manual(values = cols)
    
    ## pooled results
    bm_pool <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_pool.rds"))
    bm_pool$Algorithm <- factor(bm_pool$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    bm_pool$Coord <- factor(bm_pool$Coord, levels = c(paste0(rep(c(1:23, "X", "Y"), each = 2), c("p", "q"))))

    ## perform a systematic t-test between Uncorrected vs other algorithms
    # bm_pool %>%
    #     filter(Algorithm != "Uncorrected") %>%
    #     inner_join(bm_pool %>% 
    #         filter(Algorithm == "Uncorrected") %>%
    #         select(-Algorithm), by = c("Coord")) %>%
    #     group_by(Algorithm) %>%
    #     summarise(p = t.test(est.x, est.y)$p.value)

    p_pool_sum <- ggplot(bm_pool, aes(x = Algorithm, y = est, fill = Algorithm)) +
        geom_boxplot() +
        geom_jitter(width = 0.1, size = 2) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "p(intra-arm > inter-arm cosine similarity)") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1),
            axis.ticks.length = unit(0.5, "cm"),
            axis.title = element_text(size = 30, color = 'black'),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio = 1,
            text = element_text(family = "Arial"),
            legend.position = "none") +
        scale_fill_manual(values = cols)
    
    ## TP53 results
    bm_TP53 <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_TP53.rds"))
    bm_TP53$Algorithm <- factor(bm_TP53$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    bm_TP53$Coord <- factor(bm_TP53$Coord, levels = c(paste0(rep(c(1:23, "X", "Y"), each = 2), c("p", "q"))))

    ## perform a systematic t-test between Uncorrected vs other algorithms
    # bm_TP53 %>%
    #     filter(Algorithm != "Uncorrected" & Status == "WT") %>%
    #     inner_join(bm_TP53 %>% 
    #         filter(Algorithm == "Uncorrected" & Status == "WT") %>%
    #         select(-Algorithm), by = c("Coord")) %>%
    #     group_by(Algorithm) %>%
    #     summarise(p = t.test(est.x, est.y)$p.value)
    
    # bm_TP53 %>%
    #     filter(Algorithm != "Uncorrected" & Status == "Mut") %>%
    #     inner_join(bm_TP53 %>% 
    #         filter(Algorithm == "Uncorrected" & Status == "Mut") %>%
    #         select(-Algorithm), by = c("Coord")) %>%
    #     group_by(Algorithm) %>%
    #     summarise(p = t.test(est.x, est.y)$p.value)

    p_TP53_sum <- ggplot(bm_TP53, aes(x = Algorithm, y = est, fill = Status)) +
        geom_boxplot(outlier.shape=NA) +
        geom_point(position=position_jitterdodge(jitter.width = 0.1)) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1),
            axis.ticks.length = unit(0.5, "cm"),
            axis.title = element_text(size = 30, color = 'black'),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio = 1,
            text = element_text(family = "Arial"),
            legend.key.size = unit(1.5, 'cm'),
            legend.text = element_text(size = 20),
            legend.position = c(0.15, 0.95),
            legend.background = element_rect(fill = NA, color = NA)) +
        scale_fill_manual(labels = c("TP53 mut", "TP53 wt"), 
            values = c("#1F78B4", "#A6CEE3"), name = "")

    # Create panel
    panel_bias <- p_cn_abs_unexpr + p_rec_ampl_noexpr_bg_noexpr +
        p_pool_sum + p_TP53_sum +
        plot_layout(widths = c(1, 1), ncol = 2) +
        plot_annotation(tag_levels = 'A') &
        theme(plot.tag.position = c(0, 1),
            plot.tag = element_text(size = 50, face = "bold", family = "Arial"))
    ggsave(panel_bias, filename = paste0("results/panels/bias/bias_unexpr_", lib, "_bg_unexpr.pdf"), width = 25, height = 25, units = "in", dpi = 300)
}
