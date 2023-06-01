library(tidyverse)

load('data/biomarkers/MoBEM.RData')


# create biomarker reference
biomarkers <- MoBEM %>%
    as.data.frame() %>%
    rownames_to_column("CFE") %>%
    pivot_longer(-CFE, names_to = "COSMICID", values_to = "Status") %>%
    mutate(COSMICID = as.numeric(COSMICID)) %>%
    ## assign model type
    inner_join(., read_csv('data/Model.csv') %>%
        select(ModelID, COSMICID, DepmapModelType), by = "COSMICID") %>%
    ## assign CFE type (somatic mutation, copy number alteration or hypermethylation)
    mutate(CFE_type = ifelse(grepl("mut", CFE), "mut", 
        ifelse(grepl("HypMET", CFE), "HypMET", "RACS"))) %>%
    ## strip CFE name from bracketed gene names and other infos 
    ## (gain/loss prefix for copy number alteration, mut suffix for mutation and HypMET suffix for hypermethylation)
    mutate(CFE_stripped = ifelse(grepl("mut", CFE), gsub("_mut", "", CFE),
        ifelse(grepl("HypMET", CFE), gsub("_HypMET", "", CFE),
            ifelse(grepl("gain", CFE), gsub("gain:", "", CFE), gsub("loss:", "", CFE))))) %>%
    separate(CFE_stripped, into = c("CFE_stripped", "code"), sep = " \\(") %>%
    select(-code) %>%
    separate(CFE_stripped, into = c("CFE_stripped", "code"), sep = "\\(") %>%
    select(-code)


# save reference
saveRDS(biomarkers, "data/biomarkers/biomarkers.rds")
