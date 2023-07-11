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
    ## Nº significant biomarkers (all CFEs on strongly selective dependencies)
    sig_biomarkers_ssd <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_ssd.rds"))
    sig_biomarkers_ssd$Algorithm <- factor(sig_biomarkers_ssd$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    p_sig_biomark_ssd <- ggplot(sig_biomarkers_ssd, aes(x = Algorithm, y = n_sig_biomark, fill = Algorithm)) +
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
        scale_fill_manual(values = cols)

    ## Nº significant biomarkers (gain-of-function CFEs within oncogene)
    sig_biomarkers_onco <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_sig_biomarkers_onco.rds"))
    sig_biomarkers_onco$Algorithm <- factor(sig_biomarkers_onco$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

    p_sig_biomark_onco <- ggplot(sig_biomarkers_onco, aes(x = Algorithm, y = n_sig_biomark, fill = Algorithm)) +
        geom_bar(stat = "identity") +
        labs(x = "", y = "Nº significant associations\n(gain-of-function CFEs)") +
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
        scale_fill_manual(values = cols)
    
    ## Assemble panel
    panel <- p_sig_biomark_ssd + p_sig_biomark_onco + 
        plot_layout(ncol = 2) + plot_annotation(tag_levels = "A") &
        theme(plot.tag.position = c(0, 1),
            plot.tag = element_text(size = 40, face = "bold", family = "Arial"))
    ggsave(panel, filename = paste0("results/panels/data_quality/biomarkers_", lib, ".pdf"), width = 25, height = 12)
}
