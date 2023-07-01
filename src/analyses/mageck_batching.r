library(tidyverse)

libs <- c("Avana", "KY")

# iterate over algorithms and libraries
for (lib in libs){

    for (i in seq(1:100)){
        MAGeCK_tmp <- read_csv(paste0("data/corrected/MAGeCK/", lib, "_gene_MAGeCK_", i,".csv")) %>%
            column_to_rownames("Gene")

        if (i == 1){
            MAGeCK_mean <- MAGeCK_tmp
        } else {
            MAGeCK_tmp <- MAGeCK_tmp[rownames(MAGeCK_mean), colnames(MAGeCK_mean)]
            MAGeCK_mean <- MAGeCK_mean + MAGeCK_tmp
        }
    }

    MAGeCK_mean <- MAGeCK_mean/100

    for (i in seq(1:100)){
        MAGeCK_tmp <- read_csv(paste0("data/corrected/MAGeCK/", lib, "_gene_MAGeCK_", i,".csv")) %>%
            column_to_rownames("Gene")
        MAGeCK_tmp <- MAGeCK_tmp[rownames(MAGeCK_mean), colnames(MAGeCK_mean)]

        if (i == 1){
            MAGeCK_sd <- (MAGeCK_tmp - MAGeCK_mean)^2/(100-1)
        } else {
            MAGeCK_sd <- MAGeCK_sd + (MAGeCK_tmp - MAGeCK_mean)^2/(100-1)
        }
    }

    MAGeCK_sd <- sqrt(MAGeCK_sd) %>%
        as_tibble(rownames="Gene") %>%
        pivot_longer(-Gene, names_to = "ModelID", values_to = "sd")

    ## save results
    p_sd <- ggplot(MAGeCK_sd, aes(y = log10(sd), x = factor(0))) +
        geom_boxplot() +
        labs(x = "", y = "log10(sd)", title = paste0("MAGeCK point-wise sd ", lib, " dataset")) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(size = 14, hjust = 0.5),
            aspect.ratio = 1)
    ggsave(p_sd, filename = paste0("results/analyses/mageck_batching/", lib, "_MAGeCK_sd.pdf"), width = 20, height = 10, units = "in", dpi = 300)
    saveRDS(MAGeCK_sd, paste0("results/analyses/mageck_batching/", lib, "_MAGeCK_sd.rds"))
}
