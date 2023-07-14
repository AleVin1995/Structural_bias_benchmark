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
    ## Recall at 5% FDR
    recall_gene_sets <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_gene_sets.rds")) %>%
        filter(Gene_Set %in% c("ess_genes", "noness_genes", "msigdb_genes")) %>%
        mutate(Gene_Set = ifelse(Gene_Set == "ess_genes", "Essential genes",
            ifelse(Gene_Set == "noness_genes", "Nonessential genes", "MSigDB genes")))
    recall_gene_sets$Algorithm <- factor(recall_gene_sets$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    p_gene_sets <- ggplot(recall_gene_sets, aes(x = Algorithm, y = Recall, fill = Algorithm)) +
            geom_boxplot() +
            labs(x = "", y = "Recall at 5% FDR") +
            theme_bw() +
            theme(
                strip.text.x = element_text(size = 25, color = 'black'),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size = 25, color = 'black'),
                axis.title = element_text(size = 30),
                plot.title = element_text(size = 32, hjust = 0.5),
                axis.text.x = element_blank(),
                plot.margin = grid::unit(c(2,5,2,5), "cm"),
                legend.position = "none") +
            facet_wrap(~Gene_Set, scales = "free_y") +
            theme(strip.background = element_blank(),
                strip.text = element_text(size = 10, face = "bold")) +
            scale_fill_manual(values = cols)
    
    ## Recall curve amplified genes
    rec_ampl <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl.rds"))
    rec_ampl$Algorithm <- factor(rec_ampl$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    p_rec_ampl <- ggplot(rec_ampl, aes(x = Algorithm, y = Recall, fill = Algorithm)) +
        geom_boxplot() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "Recall curve") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 25, color = 'black'),
            axis.title = element_text(size = 30),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.margin = grid::unit(c(1,1,1,1), "cm"),
            legend.position = "none") +
        scale_fill_manual(values = cols)

    ## Recall curve amplified genes (unexpressed)
    rec_ampl_noexpr <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl_noexpr.rds"))
    rec_ampl_noexpr$Algorithm <- factor(rec_ampl_noexpr$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    
    p_rec_ampl_noexpr <- ggplot(rec_ampl_noexpr, aes(x = Algorithm, y = Recall, fill = Algorithm)) +
        geom_boxplot() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 25, color = 'black'),
            axis.title = element_text(size = 30),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.margin = grid::unit(c(1,1,1,1), "cm"),
            legend.title = element_text(size = 24),
            legend.key.size = unit(1, 'cm'),
            legend.text = element_text(size = 14)) +
        scale_fill_manual(values = cols)
    
    ## Assemble panel
    panel <- p_gene_sets / (p_rec_ampl | p_rec_ampl_noexpr) + 
        plot_annotation(tag_levels = "A") &
        theme(plot.tag.position = c(0, 1),
            plot.tag = element_text(size = 40, face = "bold", family = "Arial"))
    ggsave(panel, filename = paste0("results/panels/data_quality/gene_sets_", lib, ".pdf"), width = 25, height = 15, dpi = 300)
}
