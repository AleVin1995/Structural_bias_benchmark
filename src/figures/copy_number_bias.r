library(tidyverse)

source("src/figures/utils.r")

libs <- c("Avana", "KY")

# iterate over algorithms and libraries
for (lib in libs){
    ## All genes
    dfs <- readRDS(paste0("results/analyses/cn_correction/", lib, "_cn_abs.rds"))

    p_cn_abs <- ggplot(dfs, aes(x = as.factor(CN_abs), y = LFC, fill = Algorithm)) +
        stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + 
        labs(x = "Absolute CN", y = "LFC", title = paste0("CN correction on ", lib, " dataset")) +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(hjust = 1),
            aspect.ratio = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        facet_wrap(~Algorithm, scales = "free", ncol = 4)
    ggsave(p_cn_abs, filename = paste0("results/analyses/cn_correction/", lib, "_cn_abs.pdf"), width = 20, height = 10, units = "in", dpi = 300)

    ## LFC vs CN effect size (all genes)
    dfs_es <- dfs %>%
        group_by(Algorithm, CN_abs) %>%
        mutate(es = mean(LFC)/sd(LFC)) %>%
        ungroup() %>%
        select(Algorithm, es) %>%
        distinct()
    
    p_es <- ggplot(dfs_es, aes(x = Algorithm, y = es, fill = Algorithm)) +
        geom_boxplot() +
        labs(x = "Method", y = "Effect size", title = paste0("LFC vs CN abs effect size ", lib, " dataset")) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(hjust = 0.5),
            aspect.ratio = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black")
    ggsave(p_es, filename = paste0("results/analyses/cn_correction/", lib, "_cn_abs_es.pdf"), width = 20, height = 10, units = "in", dpi = 300)


    ## Consider only unexpressed genes (TPM < 1)
    dfs_unexpr <- dfs %>%
        inner_join(read_csv("data/OmicsExpressionProteinCodingGenesTPMLogp1.csv") %>%
            rename(ModelID = colnames(.)[1]) %>%
            pivot_longer(-ModelID, names_to = "Gene", values_to = "TPM") %>%
            separate(Gene, into = c("Gene", "Code"), sep = " \\(") %>%
            group_by(Gene) %>%
            mutate(mean_TPM = mean(TPM, na.rm = TRUE)) %>%
            ungroup() %>%
            filter(mean_TPM < 1) %>%
            select(-c(Code, mean_TPM)))
    
    p_cn_abs_unexpr <- ggplot(dfs_unexpr, aes(x = as.factor(CN_abs), y = LFC, fill = Algorithm)) +
        stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + 
        labs(x = "Absolute CN (TPM < 1)", y = "LFC", title = paste0("CN correction on ", lib, " dataset")) +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(hjust = 1),
            aspect.ratio = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        facet_wrap(~Algorithm, scales = "free", ncol = 4)
    ggsave(p_cn_abs_unexpr, filename = paste0("results/analyses/cn_correction/", lib, "_cn_abs_unexpr.pdf"), width = 20, height = 10, units = "in", dpi = 300)

    ## LFC vs CN effect size (unexpressed genes)
    dfs_es_unexpr <- dfs_unexpr %>%
        group_by(Algorithm, CN_abs) %>%
        mutate(es = mean(LFC)/sd(LFC)) %>%
        ungroup() %>%
        select(Algorithm, es) %>%
        distinct()
    
    p_es_unexpr <- ggplot(dfs_es_unexpr, aes(x = Algorithm, y = es, fill = Algorithm)) +
        geom_boxplot() +
        labs(x = "Method", y = "Effect size", title = paste0("LFC vs CN abs effect size ", lib, " dataset")) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(hjust = 0.5),
            aspect.ratio = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black")
    ggsave(p_es_unexpr, filename = paste0("results/analyses/cn_correction/", lib, "_cn_abs_es_unexpr.pdf"), width = 20, height = 10, units = "in", dpi = 300)
}