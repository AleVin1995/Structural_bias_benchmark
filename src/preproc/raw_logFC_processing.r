library(data.table)
library(dtplyr)
library(tidyverse)

libraries <- c("Avana", "KY")

# Filter screens/sequences that passed QC
ScreenSequenceMap <- read_csv("data/ScreenSequenceMap.csv") %>% 
    filter(PassesQC == TRUE)
SequenceQCReport <- read_csv("data/AchillesSequenceQCReport.csv") %>% 
    filter(PassesQC == TRUE & QCStatus == "PASS" & ScreenPassesQC == TRUE)

for (lib in libraries){
    # Read in the data
    df <- read_csv(paste0("data/", lib, "LogfoldChange.csv"))
    colnames(df)[1] <- "sgRNA"

    # Filter guide map that passed QC
    GuideMap <- read_csv(paste0("data/", lib, "GuideMap.csv")) %>% 
        filter(UsedByChronos == TRUE)

    print(paste0("Processing ", lib, "..."))

    # Format the data
    df <- df %>% 
        as.data.table() %>% ## speed up the process pt 1
        lazy_dt() %>% ## speed up the process pt 2
        pivot_longer(cols = -sgRNA, names_to = "SequenceID", values_to = "LogFoldChange") %>% 
        filter(sgRNA %in% GuideMap$sgRNA) %>% ## Filter out sgRNAs that did not pass QC
        filter(SequenceID %in% SequenceQCReport$SequenceID) %>% ## Filter out sequences that did not pass QC
        left_join(ScreenSequenceMap, by = "SequenceID") %>%
        select(sgRNA, SequenceID, ModelID, LogFoldChange) %>%
        group_by(sgRNA, ModelID) %>%
        mutate(meanLogFoldChange = mean(LogFoldChange, na.rm = TRUE)) %>% ## Take the mean of replicates
        ungroup() %>%
        select(sgRNA, ModelID, meanLogFoldChange) %>%
        distinct() %>%
        as_tibble()
    
    ## sgRNA-level LFC data
    df_sgRNA <- df %>%
        pivot_wider(names_from = ModelID, values_from = meanLogFoldChange)
    
    # Gene-level LFC data
    df_gene <- df %>%
        as.data.table() %>% ## speed up the process pt 1
        lazy_dt() %>% ## speed up the process pt 2
        left_join(GuideMap, by = "sgRNA") %>%
        separate(col = "Gene", sep = " \\(", into = c("Gene", "code")) %>%
        select(Gene, ModelID, meanLogFoldChange) %>%
        rename(LogFoldChange = meanLogFoldChange) %>%
        group_by(Gene, ModelID) %>%
        mutate(medianLogFoldChange = median(LogFoldChange, na.rm = TRUE)) %>%
        ungroup() %>%
        select(Gene, ModelID, medianLogFoldChange) %>%
        distinct() %>%
        as_tibble() %>%
        pivot_wider(names_from = ModelID, values_from = medianLogFoldChange)
    
    print(paste0("Finished to process. Saving ", lib, "..."))

    # Save the data
    write_csv(df_sgRNA, paste0("data/raw/", lib, "_sgrna_raw_LFC.csv"), col_names = TRUE)
    write_csv(df_gene, paste0("data/raw/", lib, "_gene_raw_LFC.csv"), col_names = TRUE)
}
