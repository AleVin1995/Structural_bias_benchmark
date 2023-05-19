library(tidyverse)


# Function to format omics data (Copy number and expression)
format_omics <- function(df, values_to){
  df <- df %>% 
    pivot_longer(-colnames(.)[1], names_to = 'GENE_NAME', values_to = values_to) %>%
    rename(SAMPLE_NAME = colnames(.)[1]) %>%
    separate(col = "GENE_NAME", sep = " \\(", into = c("GENE_NAME", "code")) %>%
    select(-code)

    return(df)
}


# Function to combine GuideMap with omics data (Copy number and expression)
combine_map_with_omics <- function(lib, essential_genes){
    df <- read_csv(paste0('data/', lib, 'GuideMap.csv')) %>%
        filter(UsedByChronos == TRUE) %>%
        separate(col = "Gene", sep = " \\(", into = c("GENE_NAME", "code")) %>%
        select(-c(code, nAlignments, DropReason, UsedByChronos)) %>%
        separate(col = "GenomeAlignment", sep = '_', into = c("CHROMOSOME", "POSITION", "STRAND")) %>%
        select(-STRAND) %>%
        mutate(POSITION = as.numeric(POSITION)) %>%
        inner_join(read_csv(paste0('data/raw/', lib, '_sgrna_raw_LFC.csv')) %>%
            pivot_longer(-sgRNA, names_to = 'SAMPLE_NAME', values_to = 'DEPENDENCY_SCORE'), 
            by = c("sgRNA")) %>%
        left_join(read_csv('data/OmicsCNGene.csv') %>%
            format_omics(., values_to = 'CNA'),
            by = c("GENE_NAME", "SAMPLE_NAME")) %>%
        left_join(read_csv('data/OmicsExpressionProteinCodingGenesTPMLogp1.csv') %>%
            format_omics(., values_to = 'EXP'),
            by = c("GENE_NAME", "SAMPLE_NAME")) %>%
        mutate(ESSENTIAL = GENE_NAME %in% essential_genes) ## Genome-wide common essential genes
    
    return(df)
}


# Common essential genes
EssGenes <- read_csv("data/AchillesCommonEssentialControls.csv") %>% 
    separate(col = "Gene", sep = " \\(", into = c("Gene", "code")) %>%
    pull(Gene)

# Avana data
data <- combine_map_with_omics('Avana', EssGenes)
saveRDS(data, 'data/Avana_lfc_exp_cn.rds')

# KY data
data <- combine_map_with_omics('KY', EssGenes)
saveRDS(data, 'data/KY_lfc_exp_cn.rds')
