library(tidyverse)

# Load data
Model <- read_csv('data/Model.csv') %>%
    group_by(OncotreeLineage) %>%
    mutate(occurrence_lineage = n()) %>%
    ungroup() %>%
    filter(occurrence_lineage >= 15) %>% ## select lineages with >= 15 cell lines
    select(ModelID, CellLineName, OncotreeLineage) %>%
    filter(!is.na(ModelID) & !is.na(CellLineName) & !is.na(OncotreeLineage)) %>%
    distinct()

## Preprocess absolute copy number data
CN_abs <- read_csv('data/cnv_abs_copy_number_picnic_20191101.csv') %>%
    select(-model_id) %>%
    dplyr::slice(-2) %>% ## remove rows
    t() %>%
    as_tibble() %>%
    mutate_all(~replace_na(., "0")) %>%
    column_to_rownames(colnames(.)[1]) %>%
    t() %>%
    as_tibble() %>%
    column_to_rownames(colnames(.)[1]) %>%
    t() %>%
    as.data.frame() %>%
    mutate_all(as.numeric) %>%
    rownames_to_column('CellLineName') %>%
    pivot_longer(-CellLineName, names_to = 'Gene', values_to = 'CN_abs')

# Avana gene-level LFC dataset
avana_gene <- read_csv('data/raw/Avana_gene_raw_LFC.csv') %>%
  pivot_longer(-Gene, names_to = 'ModelID', values_to = 'LFC')

# Merge data
res <- inner_join(avana_gene, Model, by = 'ModelID') %>%
    inner_join(CN_abs, c('Gene', 'CellLineName')) %>%
    group_by(CellLineName) %>%
    mutate(rank = rank(-CN_abs, ties.method = 'first')) %>%
    filter(rank <= 100) %>% ## select the top (100) amplified genes
    mutate(average_CN_abs = mean(CN_abs)) %>%
    mutate(average_LFC = mean(LFC)) %>%
    ungroup() %>%
    mutate(Class = case_when(
        average_CN_abs < 2 ~ 'CN: 0-2',
        average_CN_abs < 4 ~ 'CN: 2-4',
        average_CN_abs < 6 ~ 'CN: 4-6',
        average_CN_abs < 8 ~ 'CN: 6-8',
        TRUE ~ 'CN: 8+'
    )) %>%
    left_join(res %>%  ## count number of cell lines in each lineage
        distinct(CellLineName, OncotreeLineage) %>% 
        group_by(OncotreeLineage) %>% 
        mutate(occurrence_lineage = n()) %>%
        ungroup()) 

# Plot
p <- res %>%
    filter(occurrence_lineage >= 10) %>%
    mutate(CellLineName = factor(CellLineName, levels = unique(CellLineName[order(average_LFC)]))) %>%
    ggplot(aes(x = CellLineName, y = LFC)) +
    geom_point(aes(x = CellLineName, y = average_LFC, color = Class), size = 0.5, alpha = 0.5) +
    geom_errorbar(aes(ymin = average_LFC - sd(LFC), ymax = average_LFC + sd(LFC), color = Class), width = 0.1) +
    ylim(-3, 1) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = 'bottom'
    ) +
    geom_hline(yintercept = 0, linetype = 'longdash') +
    geom_hline(yintercept = -1, linetype = 'longdash') +
    labs(
        title = paste0('Depletion of the top 100 amplified genes in each lineage'),
        x = 'Genes Ranked by Amplification',
        y = 'Gene-level LFC',
        color = 'Mean Copy Number',
    ) +
      scale_color_manual(values = c('yellow', alpha('orange', 0.7), 
        alpha('red', 0.7), alpha('brown', 0.7), alpha('black', 0.7))) +
    facet_grid(.~OncotreeLineage, scales = 'free_x', space = 'free_x') +
    theme(strip.background = element_blank(),
        strip.text.x = element_text(angle = 45, hjust = 0.5, size = 10, color = 'black'),
        strip.clip = 'off')

ggsave('results/EDA/Summary_CN_bias.pdf', p, width = 16, height = 8, dpi = 300)
