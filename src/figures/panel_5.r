library(extrafont)
library(forcats)
library(patchwork)
library(RColorBrewer)
library(tidyverse)

font_import(paths = "arial", prompt = FALSE)

cols <- c("#B3B3B3", brewer.pal(n = 7, name = "Dark2"))

# Nº significant biomarkers (all CFEs on strongly selective dependencies)
## Avana
sig_biomarkers_ssd <- readRDS("results/analyses/impact_data_quality/Avana_sig_biomarkers_ssd.rds")
sig_biomarkers_ssd$Algorithm <- factor(sig_biomarkers_ssd$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

p_sig_biomark_ssd_avana <- ggplot(sig_biomarkers_ssd, aes(x = Algorithm, y = n_sig_biomark, fill = Algorithm)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "Nº significant associations") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(0.5, "cm"),
        axis.text = element_text(size = 25, color = 'black'),
        axis.title = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1,
        plot.margin = grid::unit(c(2,2,2,2), "cm"),
        legend.position = "none") +
    scale_fill_manual(values = c("#B3B3B3", cols))

## KY
sig_biomarkers_ssd <- readRDS("results/analyses/impact_data_quality/KY_sig_biomarkers_ssd.rds")
sig_biomarkers_ssd$Algorithm <- factor(sig_biomarkers_ssd$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

p_sig_biomark_ssd_ky <- ggplot(sig_biomarkers_ssd, aes(x = Algorithm, y = n_sig_biomark, fill = Algorithm)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "Nº significant associations") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(0.5, "cm"),
        axis.text = element_text(size = 25, color = 'black'),
        axis.title = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1,
        plot.margin = grid::unit(c(2,2,2,2), "cm"),
        legend.position = "none") +
    scale_fill_manual(values = c("#B3B3B3", cols))

## Assemble panel
panel <- p_sig_biomark_ssd_avana + p_sig_biomark_ssd_ky +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 40, face = "bold", family = "Arial"))
ggsave(panel, filename = "results/panels/data_quality/sig_biomarkers.pdf", width = 20, height = 10)
