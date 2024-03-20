library(extrafont)
library(brunnermunzel)
library(RColorBrewer)
library(tidyverse)

font_import(paths = "arial", prompt = FALSE)

# load cytoband information
cytoband <- read_csv('data/cytoband_mapping.csv') %>%
    rename(Gene = symbol) %>%
    mutate(Arm = ifelse(grepl("p", map), "p", ifelse(grepl("q", map), "q", NA))) %>%
    separate(map, into = c("Chromosome", "Band"), sep = "p|q") %>%
    na.omit() %>%
    select(-Band) %>%
    distinct() %>%
    filter(Chromosome %in% c(1:22, "X", "Y"))


# perform brunner-munzel test
perform_brunnermunzel <- function(df, cytoband_info){
    df_corr <- df %>% 
            pivot_longer(-1, names_to = "ModelID", values_to = "LFC") %>%
            group_by(ModelID) %>%
            mutate(LFC = LFC - mean(LFC, na.rm = TRUE)) %>%
            ungroup() %>%
            pivot_wider(names_from = ModelID, values_from = LFC) %>%
            column_to_rownames("Gene") %>%
            t() %>%
            WGCNA::cor(use = "pairwise.complete.obs") %>%
            as_tibble(rownames = "Gene1") %>%
            pivot_longer(-Gene1, names_to = "Gene2", values_to = "Cor")
    
    ## merge with cytoband information
    df_corr <- df_corr %>%
        inner_join(cytoband_info, by = c("Gene1" = "Gene")) %>%
        rename(Chromosome1 = Chromosome, Arm1 = Arm) %>%
        inner_join(cytoband_info, by = c("Gene2" = "Gene")) %>%
        rename(Chromosome2 = Chromosome, Arm2 = Arm) %>%
        mutate(Class = ifelse(Chromosome1 == Chromosome2 & Arm1 == Arm2, "Intra", "Inter"))
    
    ## perform brunner-munzel test
    chrms <- unique(df_corr$Chromosome1)
    arms <- unique(df_corr$Arm1)

    res <- data.frame(Chromosome = rep(chrms, length(arms)), Arm = rep(arms, each = length(chrms)), est = NA) %>%
        arrange(Chromosome, Arm)

    for (chrm in chrms){
        for (arm in arms){
            # print(paste0("Processing ", chrm, " ", arm))
            df_sub <- df_corr %>%
                filter(Chromosome1 == chrm | Chromosome2 == chrm) %>%
                filter(Arm1 == arm | Arm2 == arm)
            
            if (nrow(df_sub) == 0){
                # print(paste0("No data for ", chrm, ", skipping..."))
                next
            }

            if (length(unique(df_sub$Class)) != 2){
                # print(paste0("Not two classes for chromosome ", chrm, " and arm ", arm, ", skipping..."))
                next
            }
            
            res[res$Chromosome == chrm & res$Arm == arm, "est"] <- brunnermunzel.test(Cor ~ Class, data = df_sub)$estimate %>%
                as.numeric()
        }
    }

    res <- res %>%
        na.omit()

    return(res)
}

# Load data
df <- read_csv('data/CRISPRGeneEffect23Q4.csv') %>%
    pivot_longer(cols = -1, names_to = 'Gene', values_to = 'LFC') %>%
    separate(col = 'Gene', into = c('Gene', 'code'), sep = ' \\(') %>%
    select(-code) %>%
    dplyr::rename(ModelID = '...1') %>%
    pivot_wider(names_from = 'ModelID', values_from = 'LFC') %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) ## fill na with 0

## perform brunner-munzel test
bm_pool <- df %>%
    perform_brunnermunzel(., cytoband) %>%
    unite("Coord", c("Chromosome", "Arm"), sep = "")

## save results
saveRDS(bm_pool, paste0("results/miscellanea/Chronos_23Q4_bm_pool.rds"))


# Plot results
bm_pool <- readRDS("results/analyses/proximity_bias/Avana_bm_pool.rds") %>%
    mutate(Algorithm = ifelse(Algorithm == "Chronos", "Chronos_23Q2", Algorithm)) %>%
    bind_rows(readRDS("results/miscellanea/Chronos_23Q4_bm_pool.rds") %>%
                mutate(Algorithm = "Chronos_23Q4")) %>%
        mutate(Algorithm = factor(Algorithm, levels = c("Uncorrected", "CCR", "Chronos_23Q2", "Chronos_23Q4", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))) %>%
    mutate(Coord = factor(Coord, levels = c(paste0(rep(c(1:23, "X", "Y"), each = 2), c("p", "q")))))

cols <- c("#B3B3B3", brewer.pal(n = 8, name = "Dark2"))

p_pool_sum <- ggplot(bm_pool, aes(x = Algorithm, y = est, fill = Algorithm)) +
        geom_boxplot() +
        geom_jitter(width = 0.1, size = 2) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "P(intra-arm cosine > inter)") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1),
            axis.title = element_text(size = 30, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio = 1,
            text = element_text(family = "Arial"),
            plot.margin = unit(c(2,1,1,1), "cm"),
            legend.position = "none") +
        scale_fill_manual(values = cols)
ggsave(p_pool_sum, filename = "results/miscellanea/Chronos_23Q4_bm_pool.pdf", width = 10, height = 10, dpi = 300)