library(tidyverse)

# Define list of algorithms and libraries
algos <- c("CCR", "Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK")
libs <- c("Avana", "KY")

# iterate over algorithms and libraries
for (lib in libs){
    onco_auc <- readRDS(paste0("results/analyses/impact_data_quality/", lib, "_onco_auc.rds"))

    p_oncogenes <- ggplot(onco_auc, aes(x = Algorithm, y = AUROC, fill = Algorithm)) +
        geom_bar(stat = "identity") +
        labs(x = "Method", y = "AUROC", title = "AUROC oncogenes") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_oncogenes, filename = paste0("results/analyses/impact_data_quality/", lib, "_onco_auc.pdf"), width = 10, height = 10, dpi = 300)
}
