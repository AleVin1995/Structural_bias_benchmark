library(data.table)
library(dtplyr)
library(tidyverse)

libraries <- c("Avana", "KY")

# Filter screens/sequences that passed QC
ScreenSequenceMap <- read_csv("data/ScreenSequenceMap.csv") %>% 
    filter(PassesQC == TRUE)
ScreenQCReport <- read_csv("data/AchillesScreenQCReport.csv") %>% 
    filter(PassesQC == TRUE & QCStatus == "PASS")
SequenceQCReport <- read_csv("data/AchillesSequenceQCReport.csv") %>% 
    filter(PassesQC == TRUE & QCStatus == "PASS" & ScreenPassesQC == TRUE)

for (lib in libraries){
    # Read in the data
    df <- read_csv(paste0("data/", lib, "RawReadcounts.csv"))
    colnames(df)[1] <- "sgRNA"

    # Filter guide map that passed QC
    GuideMap <- read_csv(paste0("data/", lib, "GuideMap.csv")) %>% 
        filter(UsedByChronos == TRUE)

    print(paste0("Processing ", lib, "..."))

    # Format the data
    df <- df %>% 
        as.data.table() %>% ## speed up the process pt 1
        lazy_dt() %>% ## speed up the process pt 2
        pivot_longer(cols = -sgRNA, names_to = "SequenceID", values_to = "Readcounts") %>% 
        filter(sgRNA %in% GuideMap$sgRNA) %>% ## Filter out sgRNAs that did not pass QC
        filter(SequenceID %in% SequenceQCReport$SequenceID) %>% ## Filter out sequences that did not pass QC
        left_join(ScreenSequenceMap, by = "SequenceID") %>%
        select(sgRNA, SequenceID, ModelID, pDNABatch, Replicate, ScreenType, Readcounts) %>%
        distinct() %>%
        as_tibble()
    
    print(paste0("Finished to process. Saving ", lib, "..."))

    # Save the data
    write_csv(df_sgRNA, paste0("data/raw/", lib, "_sgrna_raw_readcounts.csv"), col_names = TRUE)
}