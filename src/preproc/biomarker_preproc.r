library(tidyverse)

load('data/biomarkers/CELLector.CFEs.CNAid_decode.RData')
load('data/biomarkers/CELLector.CFEs.CNAid_mapping.RData')
load('data/biomarkers/CELLector.CFEs.HMSid_decode.RData')
load('data/biomarkers/MoBEM.RData')


# create RACS reference
racs_reference <- CELLector.CFEs.CNAid_mapping %>%
    pivot_longer(-Identifier, names_to = "DepmapModelType", values_to = "TissueID") %>%
    filter(TissueID != "") %>%
    inner_join(., CELLector.CFEs.CNAid_decode %>% 
        select(Identifier, ContainedGenes) %>%
        rename(TissueID = Identifier), by = "TissueID") %>%
    filter(ContainedGenes != "") %>%
    ## split column ContainedGenes into multiple rows
    separate_rows(ContainedGenes, sep = ",") %>%
    ## nest ContainedGenes into a list
    group_by(Identifier, DepmapModelType) %>%
    nest(ContainedGenes = ContainedGenes) %>%
    ungroup() %>%
    select(-TissueID) %>%
    rename(CFE_stripped = Identifier)


# create hypermethylation reference
hypmet_reference <- CELLector.CFEs.HMSid_decode %>%
    select(Genomic.Coordinates, Cancer.Types, GN) %>%
    rename(Identifier = Genomic.Coordinates, ContainedGenes = GN, DepmapModelType = Cancer.Types) %>%
    filter(ContainedGenes != "") %>%
    ## split column ContainedGenes into multiple rows
    separate_rows(ContainedGenes, sep = "; ") %>%
    ## nest ContainedGenes into a list
    group_by(Identifier, DepmapModelType) %>%
    nest(ContainedGenes = ContainedGenes) %>%
    ungroup() %>%
    rename(CFE_stripped = Identifier)


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
    select(-code) %>%
    ## add RACS and hypermethylation references
    left_join(., bind_rows(racs_reference, hypmet_reference), by = c("CFE_stripped", "DepmapModelType")) %>%
    ## add CFEs of sumatic mutations
    unnest(cols = c(ContainedGenes), keep_empty = TRUE) %>%
    mutate(ContainedGenes = ifelse(is.na(ContainedGenes), CFE_stripped, ContainedGenes)) %>%
    group_by(CFE_stripped, DepmapModelType) %>%
    nest(ContainedGenes = ContainedGenes) %>%
    ungroup() %>%
    ## at least 5 samples per CFE with present/absent status
    group_by(CFE, DepmapModelType, Status) %>%
    mutate(sample_size = n()) %>%
    ungroup() %>%
    group_by(CFE, DepmapModelType) %>%
    mutate(Occurrence = n()) %>%
    filter(Occurrence-sample_size >= 5 & sample_size >= 5) %>%
    ungroup() %>%
    select(-sample_size, -Occurrence)


# save reference
saveRDS(biomarkers, "data/biomarkers/biomarkers.rds")
