library(extrafont)
library(forcats)
library(patchwork)
library(RColorBrewer)
library(tidyverse)

font_import(paths = "arial", prompt = FALSE)

libs <- c("Avana", "KY")
cols <- c("#B3B3B3", brewer.pal(n = 7, name = "Dark2"))
cols <- c(cols[1:3], "#EE8208", cols[4:8])

# iterate over algorithms and libraries
for (lib in libs){
    ## AUROC ess vs noness genes
    aurocs <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_AUROC.rds"))
    aurocs$Algorithm <- factor(aurocs$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    
    p_aurocs <- ggplot(aurocs, aes(x = Algorithm, y = AUROC, fill = Algorithm)) +
        geom_violin() +
        geom_boxplot(width=0.1, color="black", outlier.shape = NA) +
        labs(x = "", y = "AUROC\n(Common essential/nonessential genes)") +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_blank(),
            axis.title = element_text(size = 30, color = 'black'),
            axis.ticks.length = unit(0.5, "cm"),
            aspect.ratio = 1,
            legend.position = "none") +
        scale_fill_manual(values = cols)

    ## AUPRC ess vs noness genes
    auprcs <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_AUPRC.rds"))
    auprcs$Algorithm <- factor(auprcs$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    p_auprcs <- ggplot(auprcs, aes(x = Algorithm, y = AUPRC, fill = Algorithm)) +
        geom_violin() +
        geom_boxplot(width=0.1, color="black", outlier.shape = NA) +
        labs(x = "", y = "AUPRC\n(Common essential/nonessential genes)") +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_blank(),
            axis.title = element_text(size = 30, color = 'black'),
            axis.ticks.length = unit(0.5, "cm"),
            aspect.ratio = 1,
            legend.position = "none") +
        scale_fill_manual(values = cols)

    ## Ess - noness gene sets separation
    gene_sep <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_gene_sets_separation.rds"))
    gene_sep$Algorithm <- factor(gene_sep$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    p_gene_sep <- ggplot(gene_sep, aes(x = Algorithm, y = Separation, fill = Algorithm)) +
        geom_boxplot() +
        labs(x = "", y = "NNMD\n(Common essential/nonessential genes)") +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_blank(),
            axis.title = element_text(size = 30, color = 'black'),
            axis.ticks.length = unit(0.5, "cm"),
            aspect.ratio = 1,
            legend.position = "none") +
        scale_fill_manual(values = cols)

    # Nº significant biomarkers (all CFEs on strongly selective dependencies)
    sig_biomarkers_ssd <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_ssd.rds"))
    sig_biomarkers_ssd$Algorithm <- factor(sig_biomarkers_ssd$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    p_sig_biomark_ssd <- ggplot(sig_biomarkers_ssd, aes(x = Algorithm, y = n_sig_biomark, fill = Algorithm)) +
        geom_bar(stat = "identity") +
        labs(x = "", y = "Nº significant associations") +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.length = unit(0.5, "cm"),
            axis.text = element_text(size = 25, color = 'black'),
            axis.title = element_text(size = 30),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1,
            legend.position = "none") +
        scale_fill_manual(values = cols)

    ## Recall at 5% FDR
    recall_gene_sets <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_gene_sets.rds")) %>%
        filter(Gene_Set %in% c("ess_genes", "noness_genes", "msigdb_genes")) %>%
        mutate(Gene_Set = ifelse(Gene_Set == "ess_genes", "Common essential genes",
            ifelse(Gene_Set == "noness_genes", "Nonessential genes", "MsigDB genes")))
    recall_gene_sets$Algorithm <- factor(recall_gene_sets$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    p_ess_gene_sets <- ggplot(recall_gene_sets %>%
        filter(Gene_Set != "Nonessential genes"), 
        aes(x = Algorithm, y = Recall, fill = Algorithm)) +
            geom_boxplot() +
            labs(x = "", y = "Recall at 5% FDR") +
            theme_bw() +
            theme(
                strip.text.x = element_text(size = 25, color = 'black'),
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size = 25, color = 'black'),
                axis.title = element_text(size = 30),
                axis.ticks.length = unit(0.5, "cm"),
                plot.title = element_text(size = 32, hjust = 0.5),
                axis.text.x = element_blank(),
                legend.position = "none") +
            facet_wrap(~Gene_Set, scales = "free_y") +
            theme(strip.background = element_blank(),
                strip.text = element_text(size = 10, face = "bold")) +
            scale_fill_manual(values = cols)
    
    p_noness_gene_sets <- ggplot(recall_gene_sets %>%
        filter(Gene_Set == "Nonessential genes"), 
        aes(x = Algorithm, y = Recall, fill = Algorithm)) +
            geom_boxplot() +
            labs(x = "", y = "Recall at 5% FDR") +
            theme_bw() +
            theme(
                strip.text.x = element_text(size = 25, color = 'black'),
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size = 25, color = 'black'),
                axis.title = element_text(size = 30),
                axis.ticks.length = unit(0.5, "cm"),
                plot.title = element_text(size = 32, hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none") +
            facet_wrap(~Gene_Set, scales = "free_y") +
            theme(strip.background = element_blank(),
                strip.text = element_text(size = 10, face = "bold")) +
            scale_fill_manual(values = cols)
    
    ## Assemble panel
    panel <- p_aurocs + p_auprcs + p_gene_sep + p_ess_gene_sets +
        p_noness_gene_sets + p_sig_biomark_ssd + 
        plot_layout(ncol = 2) + plot_annotation(tag_levels = "A") &
        theme(plot.tag.position = c(0, 1),
            plot.tag = element_text(size = 40, face = "bold", family = "Arial"))
    ggsave(panel, filename = paste0("results/panels/data_quality/general_", lib, ".pdf"), width = 20, height = 30, dpi = 300)
}
