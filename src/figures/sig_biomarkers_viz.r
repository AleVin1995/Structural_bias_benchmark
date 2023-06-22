library(tidyverse)

libs <- c("Avana", "KY")

# iterate over algorithms and libraries
for (lib in libs){
    biomarkers_lib_ssd <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_ssd.rds"))

    p_sig_biomark_ssd <- ggplot(biomarkers_lib_ssd, aes(x = Algorithm, y = n_sig_biomark, fill = Algorithm)) +
        geom_bar(stat = "identity") +
        labs(x = "Method", y = "Nº significant associations (Strong Selective Dependencies)", 
            title = "Nº significant associations") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_sig_biomark_ssd, filename = paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_ssd.pdf"), width = 10, height = 10)

    biomarkers_lib_onco <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_onco.rds"))

    p_sig_biomark_onco <- ggplot(biomarkers_lib_onco, aes(x = Algorithm, y = n_sig_biomark, fill = Algorithm)) +
        geom_bar(stat = "identity") +
        labs(x = "Method", y = "Nº significant associations (gain-of-function CFEs)", 
            title = "Nº significant associations") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_sig_biomark_onco, filename = paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_onco.pdf"), width = 10, height = 10)
}