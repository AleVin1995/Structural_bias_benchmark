library(CRISPRcleanR)
library(data.table)
library(dtplyr)
library(ggpubr)
library(tidyverse)

# Pick a model
get_model_id <- function(df, ModelName){
  CellID <- df %>%
    filter(CellLineName == ModelName) %>%
    select(ModelID) %>%
    pull()

  return(CellID)
}

# Cell line specific CN bias
plot_CN_bias_CellLine <- function(model_info, CN_df, LFC_df, ModelName = 'HT-29',
  axis.text.size = 10, axis.title.size = 12, plot.title.size = 14){
  ## referece essential genes
  data(EssGenes.DNA_REPLICATION_cons)
  data(EssGenes.HISTONES)
  data(EssGenes.KEGG_rna_polymerase)
  data(EssGenes.PROTEASOME_cons)
  data(EssGenes.SPLICEOSOME_cons)
  data(EssGenes.ribosomalProteins)

  ## get LFC and CN data for specific cell line
  model_id <- get_model_id(model_info, ModelName)

  model_CN <- CN_df[ModelName, , drop = FALSE] %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(model = colnames(.)[1])

  model_LFC <- avana_gene[, model_id, drop = FALSE] %>%
    dplyr::rename(model = colnames(.)[1]) %>%
    mutate(rank = rank(model))

  common_genes <- intersect(rownames(model_LFC), rownames(model_CN))
  model_CN <- model_CN[common_genes, , drop = FALSE]
  model_LFC <- model_LFC[common_genes, , drop = FALSE]

  ## add categorical information
  model <- model_LFC %>%
    rownames_to_column('Gene') %>%
    left_join(model_CN %>% 
      rownames_to_column('Gene') %>%
      dplyr::rename(CN = model), by = 'Gene') %>%
    dplyr::rename(LFC = model) %>%
    mutate(Class = ifelse(CN <= 3, 'CN: 0-3',
            ifelse(CN <= 7, 'CN: 4-7', 'CN: 8+'))) %>%
    mutate(
      Class = ifelse(
        rownames(model_LFC) %in% EssGenes.DNA_REPLICATION_cons, 'DNA replication',
        ifelse(
          rownames(model_LFC) %in% EssGenes.HISTONES, 'Ribosomal proteins',
          ifelse(
            rownames(model_LFC) %in% EssGenes.KEGG_rna_polymerase, 'RNA polymerase',
            ifelse(
              rownames(model_LFC) %in% EssGenes.PROTEASOME_cons, 'Proteasome',
              ifelse(
                rownames(model_LFC) %in% EssGenes.SPLICEOSOME_cons, 'Spliceosome',
                .$Class
              )
            )
          )
        )
      )
    ) %>%
    mutate(fill = ifelse(
      Class %in% c('DNA replication', 'Ribosomal proteins', 'RNA polymerase', 
        'Proteasome', 'Spliceosome'), 'c1', 
        ifelse(
          Class == 'CN: 0-3', 'c2',
          ifelse(
            Class == 'CN: 4-7', 'c3', 'c4'
          )
        )
      )
    )

  ## plot ranked LFC
  p1 <- ggplot(model, aes(x = rank, y = LFC)) +
    geom_point() +
    theme_bw() +
    theme(
      axis.text = element_text(size = axis.text.size, color = 'black'),
      axis.title = element_text(size = axis.title.size),
      legend.position = 'none',
      aspect.ratio = 1
    ) +
    geom_hline(yintercept = 0, linetype = 'longdash') +
    geom_hline(yintercept = -1, linetype = 'longdash') +
    labs(
      x = 'Genes Ranked by Depletion',
      y = 'Gene-level LFC'
    ) 

  ## plot CN vs LFC
  p2 <- ggplot(model, aes(x = Class, y = LFC, fill = fill)) +
      geom_boxplot(notch = TRUE, outlier.shape = NA) +
      theme_bw() +
      theme(
        axis.text = element_text(size = axis.text.size, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = axis.title.size),
        plot.title = element_text(size = plot.title.size, hjust = 0.5),
        aspect.ratio = 1,
        legend.position = 'none'
      ) +
      geom_hline(yintercept = 0, linetype = 'longdash') +
      geom_hline(yintercept = -1, linetype = 'longdash') +
      labs(
        title = paste0(ModelName, ' cell line'),
        x = 'Gene class',
        y = 'Gene-level LFC'
      ) +
      scale_fill_manual(values = c('light blue', 'yellow', 
        alpha('orange', 0.7), alpha('red', 0.7)))

  ## merge plots
  p <- ggarrange(p1, NULL, p2, nrow = 3, heights = c(0.7, 0.1, 1))
  
  return(p)
}

# Load data
Model <- read_csv('data/Model.csv')

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
    mutate_all(as.numeric)

# Avana gene-level LFC dataset
avana_gene <- read_csv('data/raw/Avana_gene_raw_LFC.csv') %>%
  column_to_rownames('Gene')

# Plot cell line specific CN bias
p <- plot_CN_bias_CellLine(Model, CN_abs, avana_gene, ModelName = 'HT-29')
ggsave('results/EDA/HT29_CN_bias.pdf', p, width = 6, height = 10)

p <- plot_CN_bias_CellLine(Model, CN_abs, avana_gene, ModelName = 'HCC1419')
ggsave('results/EDA/HCC1419_CN_bias.pdf', p, width = 6, height = 10)
