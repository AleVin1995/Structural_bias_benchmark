library(extrafont)
library(forcats)
library(patchwork)
library(RColorBrewer)
library(tidyverse)

font_import(paths = "/group/iorio/Alessandro/CN_benchmark/arial", prompt = FALSE)

source("src/figures/utils.r")

libs <- c("Avana", "KY")
cols <- c("#B3B3B3", brewer.pal(n = 7, name = "Dark2"))

# iterate over algorithms and libraries
for (lib in libs){
    ## All genes
    dfs <- readRDS(paste0("results/analyses/cn_correction/", lib, "_cn_abs.rds"))
    dfs$Algorithm <- factor(dfs$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    p_cn_abs <- ggplot(dfs, aes(x = as.factor(CN_abs), y = LFC, fill = Algorithm)) +
        stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + 
        labs(x = "Absolute CN", y = "LFC") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_text(hjust = 0.5),
            axis.title = element_text(size = 35, color = 'black'),
            strip.background = element_blank(), 
            strip.placement = "outside",
            strip.text = element_text(size = 35, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            text = element_text(family = "Arial"),
            plot.margin = grid::unit(c(2,2,2,2), "cm")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        facet_wrap(~Algorithm, scales = "free_x", ncol = 2) +
        scale_fill_manual(values = cols)

    ## LFC vs CN effect size (all genes)
    dfs_es <- dfs %>%
        group_by(Algorithm, CN_abs) %>%
        mutate(es = mean(LFC)/sd(LFC)) %>%
        ungroup() %>%
        select(Algorithm, es) %>%
        distinct()
    
    p_es <- ggplot(dfs_es, aes(x = Algorithm, y = es, color = Algorithm)) +
        geom_jitter(width = 0.15, size = 5) +
        geom_point(aes(x = Algorithm, y = es), 
            data = dfs_es %>% group_by(Algorithm) %>% summarize(es = mean(es)), 
            size = 7,
            shape = 23,
            fill = "black",
            color = "black") +
        labs(x = "", y = "Effect size") +
        theme_bw() +
        theme(
            axis.ticks.x = element_blank(),
            axis.text = element_text(size = 32, color = 'black'),
            axis.title = element_text(size = 35, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1),
            text = element_text(family = "Arial"),
            plot.margin = grid::unit(c(2,2,2,2), "cm"),
            legend.position = "none") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        scale_color_manual(values = cols)


    ## Consider only unexpressed genes (TPM < 1)
    dfs_unexpr <- readRDS(paste0("results/analyses/cn_correction/", lib, "_cn_abs_tpm.rds"))
    dfs_unexpr$Algorithm <- factor(dfs_unexpr$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    
    p_cn_abs_unexpr <- ggplot(dfs_unexpr, aes(x = as.factor(CN_abs), y = LFC, fill = Algorithm)) +
        stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + 
        labs(x = "Absolute CN (mean TPM < 1)", y = "LFC") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 30, color = 'black'),
            axis.text.x = element_text(hjust = 1),
            axis.title = element_text(size = 35, color = 'black', hjust = 0.5),
            strip.background = element_blank(), 
            strip.placement = "outside",
            strip.text = element_text(size = 35, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            text = element_text(family = "Arial"),
            plot.margin = grid::unit(c(2,2,2,2), "cm")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        facet_wrap(~Algorithm, scales = "free_x", ncol = 2) +
        scale_fill_manual(values = cols)

    ## LFC vs CN effect size (unexpressed genes)
    dfs_es_unexpr <- dfs_unexpr %>%
        group_by(Algorithm, CN_abs) %>%
        mutate(es = mean(LFC)/sd(LFC)) %>%
        ungroup() %>%
        select(Algorithm, es) %>%
        distinct()
    
    p_es_unexpr <- ggplot(dfs_es_unexpr, aes(x = Algorithm, y = es, color = Algorithm)) +
        geom_jitter(width = 0.15, size = 5) +
        geom_point(aes(x = Algorithm, y = es), 
            data = dfs_es_unexpr %>% group_by(Algorithm) %>% summarize(es = mean(es)), 
            size = 7,
            shape = 23,
            fill = "black",
            color = "black") +
        labs(x = "", y = "Effect size") +
        theme_bw() +
        theme(
            axis.ticks.x = element_blank(),
            axis.text = element_text(size = 32, color = 'black'),
            axis.title = element_text(size = 35, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1),
            text = element_text(family = "Arial"),
            plot.margin = grid::unit(c(2,2,2,2), "cm"),
            legend.position = "none") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        scale_color_manual(values = cols)
    

    # Create panel
    panel_all <- p_cn_abs + p_es +
        plot_layout(widths = c(1.3, 1)) +
        plot_annotation(tag_levels = 'A') &
        theme(plot.tag.position = c(0, 1),
            plot.tag = element_text(size = 40, face = "bold", family = "Arial"))
    ggsave(panel_all, filename = paste0("results/panels/cn_bias/cn_bias_all_", lib, ".pdf"), width = 35, height = 20, units = "in", dpi = 300)
    
    panel_unexpr <- p_cn_abs_unexpr + p_es_unexpr +
        plot_layout(widths = c(1.3, 1)) +
        plot_annotation(tag_levels = 'A') &
        theme(plot.tag.position = c(0, 1),
            plot.tag = element_text(size = 50, face = "bold", family = "Arial"))
    ggsave(panel_unexpr, filename = paste0("results/panels/cn_bias/cn_bias_unexpr_", lib, ".pdf"), width = 35, height = 20, units = "in", dpi = 300)
}
