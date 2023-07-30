library(tidyverse)

# Map screen to model IDs
ScreenGeneEffect <- read_csv("data/ScreenGeneEffectUncorrected.csv") %>%
    inner_join(read_csv("data/ScreenSequenceMap.csv") %>% 
        filter(PassesQC == TRUE) %>%
        select(ScreenID, ModelID) %>%
        distinct(), by = "ScreenID") %>%
    separate(col = "ScreenID", sep = "\\.", into = c("ScreenID", "code")) %>%
    separate(col = "code", sep = "0", into = c("code", "num")) %>%
    select(-c(ScreenID, num)) %>%
    rename_with(~sub(" \\(.*$", "", .)) %>%
    group_split(code, .keep = FALSE)

# Sort screens by library
ScreenGeneEffect_Avana <- ScreenGeneEffect[[1]] %>%
    group_by(ModelID) %>%
    summarise(across(everything(), mean)) %>%
    ungroup() %>%
    column_to_rownames("ModelID") %>%
    t() %>%
    as_tibble(rownames = "Gene") %>%
    ## remove rows of NAs
    drop_na()

ScreenGeneEffect_KY <- ScreenGeneEffect[[2]] %>%
    group_by(ModelID) %>%
    summarise(across(everything(), mean)) %>%
    ungroup() %>%
    column_to_rownames("ModelID") %>%
    t() %>% 
    as_tibble(rownames = "Gene") %>%
    ## remove rows of NAs
    drop_na()

# Save the data
write_csv(ScreenGeneEffect_Avana, "data/raw/Avana_screen_gene_effect.csv", col_names = TRUE)
write_csv(ScreenGeneEffect_KY, "data/raw/KY_screen_gene_effect.csv", col_names = TRUE)
