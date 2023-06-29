library(tidyverse)


# Define list of algorithms and libraries
algos <- c("CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK")
libs <- c("Avana", "KY")

# iterate over algorithms and libraries
for (lib in libs){
    recall_gene_sets <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_gene_sets.rds"))

    p_gene_sets <- ggplot(recall_gene_sets, aes(x = Algorithm, y = Recall, fill = Algorithm)) +
            geom_boxplot() +
            labs(x = "Method", y = "Recall at 5% FDR", title = "Recall of gene sets") +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size = 10, color = 'black'),
                axis.title = element_text(size = 12),
                plot.title = element_text(size = 14, hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1),
                aspect.ratio = 1) +
            facet_wrap(~Gene_Set, scales = "free_y") +
            theme(strip.background = element_blank(),
                strip.text = element_text(size = 10, face = "bold"))
        ggsave(p_gene_sets, filename = paste0("results/analyses/impact_data_quality/", lib, "_recall_gene_sets.pdf"), width = 10, height = 10, dpi = 300)
    

    aurocs <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_AUROC.rds"))

    p_aurocs <- ggplot(aurocs, aes(x = Algorithm, y = AUROC, fill = Algorithm)) +
        geom_boxplot() +
        labs(x = "Method", y = "AUROC", title = "AUROC ess vs noness genes") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_aurocs, filename = paste0("results/analyses/impact_data_quality/", lib, "_AUROC.pdf"), width = 10, height = 10, dpi = 300)


    auprcs <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_AUPRC.rds"))

    p_auprcs <- ggplot(auprcs, aes(x = Algorithm, y = AUPRC, fill = Algorithm)) +
        geom_boxplot() +
        labs(x = "Method", y = "AUPRC", title = "AUPRC ess vs noness genes") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_auprcs, filename = paste0("results/analyses/impact_data_quality/", lib, "_AUPRC.pdf"), width = 10, height = 10, dpi = 300)


    rec_ampl <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl.rds"))

    p_rec_ampl <- ggplot(rec_ampl, aes(x = Algorithm, y = Recall, fill = Algorithm)) +
        geom_boxplot() +
        labs(x = "Method", y = "Recall curve", title = "Recall curve amplified genes") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_rec_ampl, filename = paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl.pdf"), width = 10, height = 10, dpi = 300)


    rec_ampl_noexpr <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl_noexpr.rds"))
    
    p_rec_ampl_noexpr <- ggplot(rec_ampl, aes(x = Algorithm, y = Recall, fill = Algorithm)) +
        geom_boxplot() +
        labs(x = "Method", y = "Recall", title = "Amplified (unexpressed) genes") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_rec_ampl, filename = paste0("results/analyses/impact_data_quality/", lib, "_recall_ampl_noexpr.pdf"), width = 10, height = 10, dpi = 300)
}