library(extrafont)
library(forcats)
library(patchwork)
library(RColorBrewer)
library(tidyverse)

font_import(paths = "arial", prompt = FALSE)

cols <- c("#B3B3B3", brewer.pal(n = 7, name = "Dark2"))
cols <- c(cols[1:3], "#EE8208", cols[4:8])

# Recall of oncogenetic addictions
## Avana
onco_auc <- readRDS(paste0("results/analyses/impact_data_quality/Avana_onco_auc.rds"))
onco_auc$Algorithm <- factor(onco_auc$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

p_oncogenes_avana <- ggplot(onco_auc, aes(x = Algorithm, y = AUROC, fill = Algorithm)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "AUROC") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 25, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 30, color = 'black'),
        axis.ticks.length = unit(0.5, "cm"),
        aspect.ratio = 1,
        plot.margin = grid::unit(c(1,1,1,1), "cm"),
        legend.position = "none") +
    scale_fill_manual(values = cols) +
    coord_cartesian(ylim = c(0.6, max(onco_auc$AUROC)+0.01))

## KY
onco_auc <- readRDS(paste0("results/analyses/impact_data_quality/KY_onco_auc.rds"))
onco_auc$Algorithm <- factor(onco_auc$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))

p_oncogenes_ky <- ggplot(onco_auc, aes(x = Algorithm, y = AUROC, fill = Algorithm)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "AUROC") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 25, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 30, color = 'black'),
        axis.ticks.length = unit(0.5, "cm"),
        aspect.ratio = 1,
        plot.margin = grid::unit(c(1,1,1,1), "cm"),
        legend.position = "none") +
    scale_fill_manual(values = cols) +
    coord_cartesian(ylim = c(0.6, max(onco_auc$AUROC)+0.01))


## Assemble panel
panel <- p_oncogenes_avana + p_oncogenes_ky +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 40, face = "bold", family = "Arial"))
ggsave(panel, filename = "results/panels/data_quality/onco_add.pdf", width = 20, height = 10)
