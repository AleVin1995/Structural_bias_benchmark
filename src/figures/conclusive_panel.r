library(extrafont)
library(patchwork)
library(readxl)
library(RColorBrewer)
library(tidyverse)

font_import(paths = "arial", prompt = FALSE)
cols <- c("#B3B3B3", brewer.pal(n = 7, name = "Dark2"))

# Read the data
dt <- readRDS("results/analyses/impact_data_quality/Avana_recall_ampl_noexpr.rds") %>%
  group_by(Algorithm) %>%
  summarise(
    `CN Bias (A)` = mean(Recall)
  ) %>% 
  ungroup() %>%
  inner_join(
    readRDS("results/analyses/impact_data_quality/KY_recall_ampl.rds") %>%
      group_by(Algorithm) %>%
      summarise(
        `CN Bias (S)` = mean(Recall)
      ) %>% 
      ungroup()
  ) %>%
  inner_join(
    readRDS("results/analyses/proximity_bias/Avana_bm_pool.rds") %>%
      group_by(Algorithm) %>%
      summarise(
        `Proximity Bias (A)` = mean(est)
      ) %>% 
      ungroup()
  ) %>%
  inner_join(
    readRDS("results/analyses/proximity_bias/KY_bm_pool.rds") %>%
      group_by(Algorithm) %>%
      summarise(
        `Proximity Bias (S)` = mean(est)
      ) %>% 
      ungroup()
  ) %>%
  inner_join(
    readRDS("results/analyses/impact_data_quality/Avana_AUROC.rds") %>%
      group_by(Algorithm) %>%
      summarise(
        `AUROC (A)` = mean(AUROC)
      ) %>% 
      ungroup()
  ) %>%
  inner_join(
    readRDS("results/analyses/impact_data_quality/KY_AUROC.rds") %>%
      group_by(Algorithm) %>%
      summarise(
        `AUROC (S)` = mean(AUROC)
      ) %>% 
      ungroup()
  ) %>%
  inner_join(
    readRDS("results/analyses/impact_data_quality/Avana_gene_sets_separation.rds") %>%
      group_by(Algorithm) %>%
      summarise(
        `NNMD (A)` = mean(Separation)
      ) %>% 
      ungroup()
  ) %>%
  inner_join(
    readRDS("results/analyses/impact_data_quality/KY_gene_sets_separation.rds") %>%
      group_by(Algorithm) %>%
      summarise(
        `NNMD (S)` = mean(Separation)
      ) %>% 
      ungroup()
  ) %>%
  inner_join(
    readRDS("results/analyses/impact_data_quality/Avana_onco_auc.rds") %>%
    rename("AUROC oncogenes (A)" = AUROC)
  ) %>%
  inner_join(
    readRDS("results/analyses/impact_data_quality/KY_onco_auc.rds") %>%
    rename("AUROC oncogenes (S)" = AUROC)
  ) %>%
  inner_join(
    readRDS("results/analyses/impact_data_quality/Avana_sig_biomarkers_ssd.rds") %>%
    rename("Biomarkers (A)" = n_sig_biomark)
  ) %>%
  inner_join(
    readRDS("results/analyses/impact_data_quality/KY_sig_biomarkers_ssd.rds") %>%
    rename("Biomarkers (S)" = n_sig_biomark)
  ) %>%
  .[c(8,1:7),]

dt <- dt %>%
  mutate(
    Mean_CNBias = -100+100*rowMeans(cbind(.[[2]],.[[3]]))/0.5,
    Mean_ProximityBias = -100+100*rowMeans(cbind(.[[4]]/.[[4]][1],.[[5]]/.[[5]][1])),
    AUROC_EssNonEss = rowMeans(cbind(100*.[[6]]/.[[6]][1],100*.[[7]]/.[[7]][1])),
    NNMD_EssNonEss = rowMeans(cbind(100*.[[8]]/.[[8]][1],100*.[[9]]/.[[9]][1])),
    AUROC_OncoAddictions = rowMeans(cbind(100*.[[10]]/.[[10]][1],100*.[[11]]/.[[11]][1])),
    Biomarkers = rowMeans(cbind(100*.[[12]]/.[[12]][1],100*.[[13]]/.[[13]][1]))
  ) %>%
  select(
    Algorithm,
    Mean_CNBias, 
    Mean_ProximityBias, 
    AUROC_EssNonEss, 
    NNMD_EssNonEss, 
    AUROC_OncoAddictions, 
    Biomarkers
  ) %>%
  ## add Type and MoA columns
  mutate(Type = c("Uncorrected", "Unsupervised", "Supervised", "Supervised", "Supervised", 
    "Unsupervised", "Unsupervised", "Supervised")) %>%
  mutate(MoA = c("Uncorrected", "Single-screen", "Multi-screen", "Single-screen", "Single-screen", 
    "Single-screen", "Single-screen", "Multi-screen"))

# Plotting for Average bias correction
p1 <- ggplot(dt, aes(x = Mean_CNBias, y = Mean_ProximityBias)) +
  geom_point(aes(shape = Type, color = MoA), size = 9) +
  geom_text(
    label = dt$Algorithm, 
    nudge_x = 0, nudge_y = 0.9, 
    size = 7,
    check_overlap = FALSE) +
  geom_hline(yintercept = 0, color = alpha("black", 0.5), linetype = "longdash") +
  geom_vline(xintercept = 0, color = alpha("black", 0.5), linetype = "longdash") +
  theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(0.5, "cm"),
      axis.text = element_text(size = 25, color = 'black'),
      axis.title = element_text(size = 30),
      aspect.ratio = 1,
      text = element_text(family = "Arial"),
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      plot.margin = grid::unit(c(2,2,2,2), "cm"),
      legend.position = "none") +
  labs(x = 'Mean CN bias (%)', y = 'Mean proximity bias (%)', title = 'Bias correction') +
  scale_color_manual(values = c("red", "blue", "#B3B3B3"))

# Plotting for Average Impact on data quality
p2 <- ggplot(dt, aes(x = AUROC_EssNonEss, y = NNMD_EssNonEss)) +
  geom_point(aes(shape = Type, color = MoA), size = 9) +
  geom_text(
    label = dt$Algorithm, 
    nudge_x = 0, nudge_y = 0.9, 
    size = 7,
    check_overlap = FALSE) +
  geom_hline(yintercept = 100, color = alpha("black", 0.5), linetype = "longdash") +
  geom_vline(xintercept = 100, color = alpha("black", 0.5), linetype = "longdash") +
  theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(0.5, "cm"),
      axis.text = element_text(size = 25, color = 'black'),
      axis.title = element_text(size = 30),
      aspect.ratio = 1,
      text = element_text(family = "Arial"),
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      plot.margin = grid::unit(c(2,2,2,2), "cm"),
      legend.position = "none") +
  labs(x = 'Mean % of AUROC Ess/nonEss', y = 'Mean % of NNMD Ess/nonEss', title = 'Impact on data quality') +
  scale_color_manual(values = c("red", "blue", "#B3B3B3"))

# Plotting for Average Impact on data heterogeneity
p3 <- ggplot(dt, aes(x = AUROC_OncoAddictions, y = Biomarkers)) +
  geom_point(aes(shape = Type, color = MoA), size = 9) +
  geom_text(
    label = dt$Algorithm, 
    nudge_x = 0, nudge_y = 0.9, 
    size = 7,
    check_overlap = FALSE) +
  geom_hline(yintercept = 100, color = alpha("black", 0.5), linetype = "longdash") +
  geom_vline(xintercept = 100, color = alpha("black", 0.5), linetype = "longdash") +
  theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(0.5, "cm"),
      axis.text = element_text(size = 25, color = 'black'),
      axis.title = element_text(size = 30),
      aspect.ratio = 1,
      text = element_text(family = "Arial"),
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      plot.margin = grid::unit(c(2,2,2,2), "cm"),
      legend.title = element_text(size = 20, color = 'black'),
      legend.text = element_text(size = 15, color = 'black')) +
  labs(x = 'Mean % of AUROC OncoAddictions', y = 'Mean % of Biomarkers', title = 'Impact on data heterogeneity') +
  scale_color_manual(values = c("red", "blue", "#B3B3B3"))

## Assemble panel
panel <- p1 + p2 + p3 +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 40, face = "bold", family = "Arial"))
ggsave(panel, filename = "results/panels/summary_panel.pdf", width = 45, height = 15)
