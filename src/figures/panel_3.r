library(extrafont)
library(forcats)
library(patchwork)
library(RColorBrewer)
library(tidyverse)

font_import(paths = "/group/iorio/Alessandro/CN_benchmark/arial", prompt = FALSE)

libs <- c("Avana", "KY")
cols <- c("#B3B3B3", brewer.pal(n = 7, name = "Dark2"))

# iterate over algorithms and libraries
for (lib in libs){
    ## pooled results
    bm_pool <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_pool.rds"))
    bm_pool$Algorithm <- factor(bm_pool$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    bm_pool$Coord <- factor(bm_pool$Coord, levels = c(paste0(rep(c(1:23, "X", "Y"), each = 2), c("p", "q"))))

    p_pool <- ggplot(bm_pool, aes(x = Algorithm, y = est, fill = Algorithm)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "P(intra-arm cosine > inter)") +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            strip.text = element_text(size = 12, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1,
            panel.spacing.y = unit(4, "lines")) +
        facet_wrap(~Coord, ncol = 8) +
        scale_fill_manual(values = cols)
    ggsave(p_pool, filename = paste0("results/panels/proximity_bias/proximity_bias_", lib, "_bm_pool.pdf"), width = 15, height = 15, dpi = 300)

    p_pool_sum <- ggplot(bm_pool, aes(x = Algorithm, y = est, fill = Algorithm)) +
        geom_boxplot() +
        geom_jitter(width = 0.1, size = 2) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "P(intra-arm cosine > inter)") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1),
            axis.title = element_text(size = 30, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio = 1,
            text = element_text(family = "Arial"),
            plot.margin = unit(c(2, 2, 2, 2), "cm"),
            legend.position = "none") +
        scale_fill_manual(values = cols)
        
    ## TP53 results
    bm_TP53 <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_TP53.rds"))
    bm_TP53$Algorithm <- factor(bm_TP53$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    bm_TP53$Coord <- factor(bm_TP53$Coord, levels = c(paste0(rep(c(1:23, "X", "Y"), each = 2), c("p", "q"))))

    p_TP53 <- ggplot(bm_TP53, aes(x = Algorithm, y = est, fill = Status)) +
        geom_bar(stat = "identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "P(intra-arm cosine > inter)") +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            strip.text = element_text(size = 12, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1,
            panel.spacing.y = unit(4, "lines")) +
        facet_wrap(~Coord) +
        scale_fill_manual(labels = c("TP53 mut", "TP53 wt"), 
            values = c("#1F78B4", "#A6CEE3"), name = "")
    ggsave(p_TP53, filename = paste0("results/panels/proximity_bias/proximity_bias_", lib, "_bm_TP53.pdf"), width = 15, height = 15, dpi = 300)

    p_TP53_sum <- ggplot(bm_TP53, aes(x = Algorithm, y = est, fill = Status)) +
        geom_boxplot() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1),
            axis.title = element_text(size = 30, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio = 1,
            text = element_text(family = "Arial"),
            plot.margin = unit(c(2, 2, 2, 2), "cm"),
            legend.key.size = unit(1, 'cm'),
            legend.text = element_text(size = 12),
            legend.position = c(0.15, 0.9),
            legend.background = element_rect(fill = "white", color = NA)) +
        scale_fill_manual(labels = c("TP53 mut", "TP53 wt"), 
            values = c("#1F78B4", "#A6CEE3"), name = "")
    
    ## CDKN2A results
    bm_CDKN2A <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2A.rds"))
    bm_CDKN2A$Algorithm <- factor(bm_CDKN2A$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    bm_CDKN2A$Coord <- factor(bm_CDKN2A$Coord, levels = c(paste0(rep(c(1:23, "X", "Y"), each = 2), c("p", "q"))))

    p_CDKN2A <- ggplot(bm_CDKN2A, aes(x = Algorithm, y = est, fill = Status)) +
        geom_bar(stat = "identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "P(intra-arm cosine > inter)") +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            strip.text = element_text(size = 12, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1,
            panel.spacing.y = unit(4, "lines")) +
        facet_wrap(~Coord) +
        scale_fill_manual(labels = c("CDKN2A mut", "CDKN2A wt"), 
            values = c("#33A02C", "#B2DF8A"), name = "")
    ggsave(p_CDKN2A, filename = paste0("results/panels/proximity_bias/proximity_bias_", lib, "_bm_CDKN2A.pdf"), width = 15, height = 15, dpi = 300)

    p_CDKN2A_sum <- ggplot(bm_CDKN2A, aes(x = Algorithm, y = est, fill = Status)) +
        geom_boxplot() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1),
            axis.title = element_text(size = 30, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio = 1,
            text = element_text(family = "Arial"),
            plot.margin = unit(c(2, 2, 2, 2), "cm"),
            legend.key.size = unit(1, 'cm'),
            legend.text = element_text(size = 12),
            legend.position = c(0.15, 0.9),
            legend.background = element_rect(fill = "white", color = NA)) +
        scale_fill_manual(labels = c("CDKN2A mut", "CDKN2A wt"), 
            values = c("#33A02C", "#B2DF8A"), name = "")
    
    ## Assemble panel
    panel <- p_pool_sum + p_TP53_sum + p_CDKN2A_sum +
        plot_annotation(tag_levels = 'A') &
        theme(plot.tag.position = c(0, 1.05),
            plot.tag = element_text(size = 40, face = "bold", family = "Arial"))
    ggsave(panel, filename = paste0("results/panels/proximity_bias/proximity_bias_", lib, "_bm_all.pdf"), width = 30, height = 15, dpi = 300)
}
