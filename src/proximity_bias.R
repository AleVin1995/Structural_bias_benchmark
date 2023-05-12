library(circlize)
library(ComplexHeatmap)
library(preprocessCore)
library(tidyverse)

# reduce matrix size, using a summarising function (default, mean)
redim_matrix <- function(
    mat,
    target_height = 1000,
    target_width = 1000,
    summary_func = function(x) mean(x, na.rm = TRUE),
    output_type = 0.0, # vapply style
    n_core = 1 # parallel processing
) {
  if(target_height > nrow(mat) | target_width > ncol(mat)) {
    stop("Input matrix must be bigger than target width and height.")
  }
  
  seq_height <- round(seq(1, nrow(mat), length.out = target_height + 1))
  seq_width  <- round(seq(1, ncol(mat), length.out = target_width  + 1))
  
  # complicated way to write a double for loop
  do.call(rbind, parallel::mclapply(seq_len(target_height), function(i) { # i is row
    vapply(seq_len(target_width), function(j) { # j is column
      summary_func(
        mat[
          seq(seq_height[i], seq_height[i + 1]),
          seq(seq_width[j] , seq_width[j + 1] )
        ]
      )
    }, output_type)
  }, mc.cores = n_core))
}

# gene to chromosome mapping
gene2chr_mapping <- function(df, GuideMap, transpose = FALSE){
  gene2chr <- separate(data = GuideMap, col = "GenomeAlignment", into = c("chr", "pos", "strand"), sep = "_") %>%
    separate(col = 'chr', into = c('null', 'chr'), sep = 'chr') %>%
    select(-null) %>%
    separate(col = 'Gene', into = c('Gene', 'code'), sep = ' \\(') %>%
    select(-code) %>%
    group_by(Gene) %>%
    mutate(avg_pos = mean(as.numeric(pos))) %>%
    select(Gene, chr, avg_pos) %>%
    distinct() %>%
    mutate(chr = str_replace(chr, "X", "23")) %>%
    mutate(chr = str_replace(chr, "Y", "24")) %>%
    mutate(chr = as.numeric(chr))

  if (transpose){
    df <- t(df)
  }

  gene2chr <- gene2chr %>% 
    filter(Gene %in% colnames(df)) %>%
    arrange(chr, avg_pos)

  return(gene2chr)
}

# Compute pairwise correlation matrix
pairwise_corr <- function(df, gene2chr, transpose = F){
  if (transpose){
    df <- t(df)
  }

  ## normalize matrix
  df <- df[, gene2chr$Gene]

  for (i in 1:length(ncol(df))){
    df[,i] <- df[,i] - mean(df[,i])
  }

  ## compute pearson correlation
  df_corr <- WGCNA::cor(df, use = "pairwise.complete.obs")
  df_corr <- (df_corr - mean(df_corr))/sd(df_corr)*0.2
  df_corr <- round(df_corr, digits = 2)
  df_corr <- redim_matrix(df_corr, target_height = 2000, target_width = 2000)

  return(df_corr)
}

# plot pairwise correlation matrix
plot <- function(df, title, savepath='results/EDA/'){
  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  ht <- Heatmap(df, cluster_rows = F, show_row_names = F, 
                cluster_columns = F, show_column_names = F,
                name = "Corr", col = col_fun, column_title = title)
  
  ## replace " " with "_" in title
  title <- str_replace_all(title, " ", "_")

  pdf(paste0(savepath, title, '_proximity_bias.pdf'), width = 10, height = 10)
  draw(ht)
  dev.off()
}

# execute entire pipeline
pipeline <- function(df, GuideMap, title, transpose = TRUE, savepath='results/EDA/'){
  ## gene to chromosome mapping
  gene2chr <- gene2chr_mapping(df, GuideMap, transpose = transpose)

  ## Compute pairwise correlation matrix
  df_corr <- pairwise_corr(df, gene2chr, transpose = transpose)

  ## plot pairwise correlation matrix
  plot(df_corr, title, savepath)
}

# Avana gene-level LFC dataset
avana_gene <- read_csv('data/raw/Avana_gene_raw_LFC.csv') %>%
  column_to_rownames('Gene')
AvanaGuideMap <- read_csv('data/AvanaGuideMap.csv') %>%
  filter(UsedByChronos == T)

## execute pipeline
pipeline(avana_gene, AvanaGuideMap, 'Avana gene raw LFC')

# KY gene-level LFC dataset
ky_gene <- read_csv('data/raw/KY_gene_raw_LFC.csv') %>%
  column_to_rownames('Gene')
KYGuideMap <- read_csv('data/KYGuideMap.csv') %>%
  filter(UsedByChronos == T)

## execute pipeline
pipeline(ky_gene, KYGuideMap, 'KY gene raw LFC')
