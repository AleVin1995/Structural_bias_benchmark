library(circlize)
library(ComplexHeatmap)
library(CRISPRcleanR)
library(data.table)
library(dtplyr)
library(extrafont)
library(forcats)
library(patchwork)
library(preprocessCore)
library(RColorBrewer)
library(tidyverse)

font_import(paths = "/group/iorio/Alessandro/Structural_bias_benchmark/arial", prompt = FALSE)

# reduce matrix size, using a summarising function (default, mean)
redim_matrix <- function(
    mat,
    target_height = 1000,
    target_width = 1000,
    summary_func = function(x) median(x, na.rm = TRUE),
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
plot <- function(df, savepath = NULL){
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    ht <- Heatmap(df, cluster_rows = F, show_row_names = F, 
                                cluster_columns = F, show_column_names = F,
                                name = "Corr", col = col_fun)
    
    pdf(paste0(savepath), width = 40, height = 40)
    draw(ht)
    dev.off()
}

# execute entire pipeline
pipeline <- function(df, GuideMap, transpose = TRUE, savepath = NULL){
    ## gene to chromosome mapping
    gene2chr <- gene2chr_mapping(df, GuideMap, transpose = transpose)

    ## Compute pairwise correlation matrix
    df_corr <- pairwise_corr(df, gene2chr, transpose = transpose)

    plot(df_corr, savepath)
}

# Avana gene-level LFC dataset
avana_gene <- read_csv('data/raw/Avana_gene_raw_LFC.csv') %>%
    column_to_rownames('Gene')
AvanaGuideMap <- read_csv('data/AvanaGuideMap.csv') %>%
    filter(UsedByChronos == T)

## execute pipeline
pipeline(avana_gene, AvanaGuideMap, savepath = 'results/panels/EDA/proximity_bias.pdf')



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
    mutate(Class = ifelse(CN == 0, 'CN: 0',
            ifelse(CN == 1, 'CN: 1',
            ifelse(CN == 2, 'CN: 2', 
            ifelse(CN == 3, 'CN: 3',
            ifelse(CN == 4, 'CN: 4',
            ifelse(CN == 5, 'CN: 5',
            ifelse(CN == 6, 'CN: 6',
            ifelse(CN == 7, 'CN: 7', 'CN: 8+'))))))))) %>%
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
        'Proteasome', 'Spliceosome'), 'c1', 'c2')
    )

  ## plot CN vs LFC
  p1 <- ggplot(model, aes(x = Class, y = LFC, fill = fill)) +
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
        x = '',
        y = 'LFC'
      ) +
      scale_fill_manual(values = c(alpha('orange', 0.7), 'light blue'))
  
  return(p1)
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
p <- plot_CN_bias_CellLine(Model, CN_abs, avana_gene, ModelName = 'HCC1419')
ggsave(p, filename = 'results/panels/EDA/HCC1419_CN_bias.pdf', width = 6, height = 10)


# Unknown source of bias
## CCR correction pipeline
ccr_correction <- function(sgrna_lfc, 
                          guide_map,
                          cell_label,
                          outdir = NULL){
  cell_data <- sgrna_lfc %>%
    select(sgRNA, cell_label) %>%
    left_join(guide_map, by = join_by(sgRNA == CODE)) %>%
    select(sgRNA, GENES, cell_label) %>%
    rename(LFC = cell_label, gene =  GENES) %>%
    as.data.frame()

  rownames(cell_data) <- cell_data$sgRNA
  rownames(guide_map) <- guide_map$CODE

  gwSortedFCs <- ccr.logFCs2chromPos(cell_data, guide_map)
  correctedFCs <- ccr.GWclean(gwSortedFCs, display = F, verbose = -1,
    label = cell_label, saveTO = outdir)
  
  return(list(cell_data = cell_data, correctedFCs = correctedFCs))
}

## plot of unknown bias
unknown_bias <- function(res,
                        cn_abs, 
                        guide_map,
                        model_map,
                        cell_label,
                        chr = '1',
                        start = NULL,
                        end = NULL,
                        ar = 0.3,
                        min_lfc = -6,
                        max_lfc = 3,
                        min_cn = 0,
                        max_cn = 4){
  ## get model id
  cellname <- model_map %>% 
    filter(ModelID == cell_label) %>% 
    pull(CellLineName)
  
  if (length(cellname) != 1){
    stop('Cell line not found')
  }

  if (!(cellname %in% colnames(cn_abs))){
    stop('Cell line not found')
  }

  point_lfc <- res$cell_data %>%
    left_join(guide_map, by = join_by(sgRNA == CODE, gene == GENES)) %>%
    select(sgRNA, gene, CHRM, STARTpos, LFC) %>%
    rename(BP = STARTpos)
  
  if (is.null(start)){
    start <- min(point_lfc$BP)
  }

  if (is.null(end)){
    end <- max(point_lfc$BP)
  }

  point_lfc <- point_lfc %>%
    filter(CHRM == chr & BP > start & BP < end)

  seg_lfc <- res$correctedFCs$segments %>%
    filter(CHR == chr & startp > start & endp < end)
  seg_lfc[1, "startp"] <- start
  seg_lfc[nrow(seg_lfc), "endp"] <- end

  ## LFC plot
  p_lfc <- ggplot(point_lfc, aes(y = LFC, x = BP)) +
    geom_point(col = "green") +
    ylim(min_lfc, max_lfc) +
    theme_bw() +
    theme(
        aspect.ratio = ar,
        plot.margin = grid::unit(c(1,1,1,1), "cm"),
        axis.text = element_text(size = 25, color = 'black'),
        axis.title = element_text(size = 30, color = 'black')) +
    xlim(as.numeric(start), as.numeric(end)) +
    geom_segment(data = seg_lfc, aes(y = `avg.logFC`, yend = `avg.logFC`, 
      x = startp, xend = endp), col = "black", linewidth = 1)
  
  ## CN plot
  point_cn <- CN_abs %>%
    select(CODE, Gene, CHRM, STARTpos, cellname) %>%
    rename(CN = cellname, BP = STARTpos) %>%
    filter(CHRM == chr & BP > start & BP < end) %>%
    arrange(CHRM, BP)

  p_cn <- ggplot(point_cn, aes(y = CN, x = BP)) +
    geom_point(col = "red") +
    ylim(min_cn, max_cn) +
    theme_bw() +
    theme(
        aspect.ratio = ar,
        plot.margin = grid::unit(c(1,1,1,1), "cm"),
        axis.text = element_text(size = 25, color = 'black'),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 30, color = 'black')) +
    xlim(start, end) +
    geom_segment(aes(y = CN, yend = CN, x = start, xend = end), 
      col = "black", linewidth = 1)

  ## assemble plot
  p_tot <- p_cn + p_lfc + plot_layout(nrow = 2)

  return(p_tot)
}

## Load data
avana_sgrna <- read_csv('data/raw/Avana_sgrna_raw_LFC.csv')

Model <- Model %>%
  select(CellLineName, ModelID)

guide_map <- read_csv('data/AvanaGuideMap.csv') %>%
  filter(UsedByChronos == T) %>%
  select(sgRNA, GenomeAlignment, Gene) %>%
  separate(col = 'GenomeAlignment', into = c('chr', 'pos', 'strand'), sep = '_') %>%
  separate(col = 'chr', into = c('null', 'chr'), sep = 'chr') %>%
  separate(col = 'Gene', into = c('Gene', 'code'), sep = ' \\(') %>%
  select(-c(null, code)) %>%
  rename(CODE = sgRNA, GENES = Gene, CHRM = chr, STRAND = strand, STARTpos = pos) %>%
  mutate(STARTpos = as.numeric(STARTpos)) %>%
  mutate(ENDpos = STARTpos + 23, seq = CODE) %>%
  as.data.frame()
rownames(guide_map) <- guide_map$CODE

CN_abs <- CN_abs %>% 
  t() %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "Gene") %>%
  left_join(guide_map, by = join_by(Gene == GENES))

## wild-type CN, negative bias
res <- ccr_correction(avana_sgrna, guide_map, 'ACH-000585', outdir = 'test/')

p_neg <- unknown_bias(res, CN_abs, guide_map, Model, 'ACH-000585', chr = '3', start = 0, end = 5e7, ar = 0.3)
ggsave(p_neg, filename = 'results/panels/EDA/wt_CN_neg_bias.pdf', width = 20, height = 10)

## wild-type CN, positive bias
res <- ccr_correction(avana_sgrna, guide_map, 'ACH-000098')

p_pos <- unknown_bias(res, CN_abs, guide_map, Model, 'ACH-000098', chr = '18', start = 2.3e7, end = 3.2e7, min_lfc = -2, max_lfc = 2)
ggsave(p_pos, filename = 'results/panels/EDA/wt_CN_pos_bias.pdf', width = 6, height = 10)
