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
        inner_join(GuideMap, by = "sgRNA") %>% ## Filter out sgRNAs that did not pass QC
        separate(col = "Gene", sep = " \\(", into = c("Gene", "code")) %>%
        select(-code) %>%
        filter(SequenceID %in% SequenceQCReport$SequenceID) %>% ## Filter out sequences that did not pass QC
        select(sgRNA, Gene, SequenceID, Readcounts) %>%
        distinct() %>%
        as_tibble() %>%
        mutate(Readcounts = replace_na(Readcounts, 0)) ## replace NA with 0
        pivot_wider(names_from = SequenceID, values_from = Readcounts)
    
    print(paste0("Finished to process. Saving ", lib, "..."))

    # Save the data
    write_csv(df, paste0("data/raw/", lib, "_sgrna_raw_readcounts.csv"), col_names = TRUE)
}
