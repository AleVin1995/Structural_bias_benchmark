library(tidyverse)

norm_chr_arms <- function(gene_effect, chromosome_arms) {
  
  # loop through arms
  for (arm in unique(chromosome_arms$ChrArm)) {
    # get genes on that arm
    genes <- intersect(rownames(chromosome_arms)[which(chromosome_arms == arm)], rownames(gene_effect))
    
    # skip arms with small number of genes
    if (length(genes) < 5) {
      next
    }
    
    # get the arm median gene effect
    gene_meds <- apply(gene_effect[genes,], 2, median)
    med_med <- median(gene_meds)
    
    # subtract the difference from the median of median gene effects
    gene_effect[genes,] <- gene_effect[genes,] - (gene_meds - med_med)
  }
  
  return(gene_effect)
}


# load chromosome arm coordinates
chromosome_arms <- read_csv("data/chromosome_arms.csv") %>%
  separate(Gene, into = c("Gene", "code"), sep = " \\(") %>%
  select(-code) %>%
  column_to_rownames("Gene")

# Define list of libraries
libs <- c("Avana", "KY")

# iterate over libraries
for (lib in libs){
  chronos_df <- paste0("data/corrected/", lib, "_gene_Chronos.csv") %>%
    read_csv %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    dplyr::rename(Gene = colnames(.)[1]) %>% ## fill na with 0
    column_to_rownames("Gene")
  
  chronos_df_norm <- norm_chr_arms(chronos_df, chromosome_arms) %>%
    rownames_to_column("Gene")
  write_csv(chronos_df_norm, paste0("data/corrected/", lib, "_gene_AC Chronos.csv"))
}
