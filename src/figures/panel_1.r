library(cowplot)
library(extrafont)
library(forcats)
library(RColorBrewer)
library(tidyverse)

font_import(paths = "/group/iorio/Alessandro/CN_benchmark/arial", prompt = FALSE)

source("src/figures/utils.r")

libs <- c("Avana", "KY")
cols <- c("#B3B3B3", brewer.pal(n = 7, name = "Paired"))

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
            axis.title = element_text(size = 30, color = 'black'),
            strip.background = element_blank(), 
            strip.placement = "outside",
            strip.text = element_text(size = 28, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            text = element_text(family = "Arial"),
            plot.margin = grid::unit(c(5,5,5,5), "mm")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        facet_wrap(~Algorithm, scales = "free_x", ncol = 4) +
        scale_fill_manual(values = cols)

    ## LFC vs CN effect size (all genes)
    dfs_es <- dfs %>%
        group_by(Algorithm, CN_abs) %>%
        mutate(es = mean(LFC)/sd(LFC)) %>%
        ungroup() %>%
        select(Algorithm, es) %>%
        distinct()
    
    p_es <- ggplot(dfs_es, aes(x = Algorithm, y = es, color = Algorithm)) +
        geom_jitter(width = 0.15, size = 3) +
        geom_point(aes(x = Algorithm, y = es), 
            data = dfs_es %>% group_by(Algorithm) %>% summarize(es = mean(es)), 
            size = 5,
            shape = 23,
            fill = "black",
            color = "black") +
        labs(x = "", y = "Effect size") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.title = element_text(size = 30, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1),
            text = element_text(family = "Arial"),
            plot.margin = grid::unit(c(5,5,5,5), "mm")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        scale_color_manual(values = cols)


    ## Consider only unexpressed genes (TPM < 1)
    dfs_unexpr <- readRDS(paste0("results/analyses/cn_correction/", lib, "_cn_abs_tpm.rds"))
    dfs_unexpr$Algorithm <- factor(dfs_unexpr$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    
    p_cn_abs_unexpr <- ggplot(dfs_unexpr, aes(x = as.factor(CN_abs), y = LFC, fill = Algorithm)) +
        stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + 
        labs(x = "Absolute CN (TPM < 1)", y = "LFC") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_text(hjust = 1),
            axis.title = element_text(size = 30, color = 'black', hjust = 0.5),
            strip.background = element_blank(), 
            strip.placement = "outside",
            strip.text = element_text(size = 28, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            text = element_text(family = "Arial"),
            plot.margin = grid::unit(c(5,5,5,5), "mm")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        facet_wrap(~Algorithm, scales = "free_x", ncol = 4) +
        scale_fill_manual(values = cols)

    ## LFC vs CN effect size (unexpressed genes)
    dfs_es_unexpr <- dfs_unexpr %>%
        group_by(Algorithm, CN_abs) %>%
        mutate(es = mean(LFC)/sd(LFC)) %>%
        ungroup() %>%
        select(Algorithm, es) %>%
        distinct()
    
    p_es_unexpr <- ggplot(dfs_es_unexpr, aes(x = Algorithm, y = es, color = Algorithm)) +
        geom_jitter(width = 0.15, size = 3) +
        geom_point(aes(x = Algorithm, y = es), 
            data = dfs_es %>% group_by(Algorithm) %>% summarize(es = mean(es)), 
            size = 5,
            shape = 23,
            fill = "black",
            color = "black") +
        labs(x = "", y = "Effect size") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.title = element_text(size = 30, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1),
            text = element_text(family = "Arial"),
            plot.margin = grid::unit(c(5,5,5,5), "mm")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        scale_fill_manual(values = cols)
    

    # Create panel
    legend <- get_legend(p_es +
        theme(legend.box.margin = margin(0,0,0,0),
            legend.key.size = unit(2, 'cm'),
            legend.key.height = unit(2, 'cm'),
            legend.key.width = unit(2, 'cm'),
            legend.title = element_text(size=24),
            legend.text = element_text(size=22)))

    panel_all <- plot_grid(p_cn_abs, 
        p_es + theme(legend.position="none"), 
        nrow = 1, ncol = 2, align = 'h', axis = 'bt', 
        rel_widths = c(2, 1), rel_heights = c(1, 1),
        labels = "AUTO", label_size = 40, family = "Arial")
    panel_all <- plot_grid(panel_all, legend, rel_widths = c(4, .5))
    ggsave(panel_all, filename = paste0("results/panels/cn_bias_all_", lib, ".pdf"), width = 45, height = 15, units = "in", dpi = 300)
    
    panel_unexpr <- plot_grid(p_cn_abs_unexpr, 
        p_es_unexpr + theme(legend.position="none"), 
        nrow = 1, ncol = 2, align = 'h', axis = 'bt', 
        rel_widths = c(2, 1), rel_heights = c(1, 1),
        labels = "AUTO", label_size = 40, family = "Arial")
    panel_unexpr <- plot_grid(panel_unexpr, legend, rel_widths = c(4, .5))
    ggsave(panel_unexpr, filename = paste0("results/panels/cn_bias_unexpr_", lib, ".pdf"), width = 45, height = 15, units = "in", dpi = 300)
}
