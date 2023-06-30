library(tidyverse)

# Define list of algorithms and libraries
libs <- c("Avana", "KY")

# iterate over algorithms and libraries
for (lib in libs){
    gene_sep <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_gene_sets_separation.rds"))

    p_gene_sep <- ggplot(gene_sep, aes(x = Algorithm, y = Separation, fill = Algorithm)) +
        geom_boxplot() +
        labs(x = "Method", y = "Gene sets difference", title = "Gene sets difference",
            subtitle = "median(ess) - median(noness)") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_gene_sep, filename = paste0("results/analyses/impact_data_quality/", lib, "_gene_sets_separation.pdf"), width = 10, height = 10, dpi = 300)
}