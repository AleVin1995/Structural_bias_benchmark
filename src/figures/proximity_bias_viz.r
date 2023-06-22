library(tidyverse)

libs <- c("Avana", "KY")

# iterate over algorithms and libraries
for (lib in libs){
    bm_pool <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_pool.rds"))

    p_pool <- ggplot(bm_pool, aes(x = Algorithm, y = est, fill = Algorithm)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "Genome-wide") +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(~Coord, scales = "free_x")
    ggsave(p_pool, filename = paste0("results/analyses/proximity_bias/", lib, "_bm_pool.pdf"), width = 15, height = 15, dpi = 300)

    p_pool_sum <- ggplot(bm_pool, aes(x = Algorithm, y = est, fill = Algorithm)) +
        geom_boxplot() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "Genome-wide") +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_pool_sum, filename = paste0("results/analyses/proximity_bias/", lib, "_bm_pool_sum.pdf"), width = 15, height = 15, dpi = 300)
        
    ## save results
    bm_TP53 <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_TP53.rds"))

    p_TP53 <- ggplot(bm_TP53, aes(x = Algorithm, y = est, fill = Status)) +
        geom_bar(stat = "identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "TP53 mutational status") +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(~Coord, scales = "free_x")
    ggsave(p_TP53, filename = paste0("results/analyses/proximity_bias/", lib, "_bm_TP53.pdf"), width = 15, height = 15, dpi = 300)

    p_TP53_sum <- ggplot(bm_TP53, aes(x = Algorithm, y = est, fill = Status)) +
        geom_boxplot() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "TP53 mutational status") +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_TP53_sum, filename = paste0("results/analyses/proximity_bias/", lib, "_bm_TP53_sum.pdf"), width = 15, height = 15, dpi = 300)
    
    ## save results
    bm_CDKN2A <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2A.rds"))

    p_CDKN2A <- ggplot(bm_CDKN2A, aes(x = Algorithm, y = est, fill = Status)) +
        geom_bar(stat = "identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "CDKN2A mutational status") +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(~Coord, scales = "free_x")
    ggsave(p_CDKN2A, filename = paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2A.pdf"), width = 15, height = 15, dpi = 300)

    p_CDKN2A_sum <- ggplot(bm_CDKN2A, aes(x = Algorithm, y = est, fill = Status)) +
        geom_boxplot() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "CDKN2A mutational status") +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_CDKN2A_sum, filename = paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2A_sum.pdf"), width = 15, height = 15, dpi = 300)
        
    ## save results
    bm_CDKN2B <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2B.rds"))

    p_CDKN2B <- ggplot(bm_CDKN2B, aes(x = Algorithm, y = est, fill = Status)) +
        geom_bar(stat = "identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "CDKN2B mutational status") +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(~Coord, scales = "free_x")
    ggsave(p_CDKN2B, filename = paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2B.pdf"), width = 15, height = 15, dpi = 300)

    p_CDKN2B_sum <- ggplot(bm_CDKN2B, aes(x = Algorithm, y = est, fill = Status)) +
        geom_boxplot() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "Method", y = "P(intra-arm cosine > inter)", title = "CDKN2B mutational status") +
        theme_bw() +
        theme(
            strip.background = element_blank(), 
            strip.placement = "outside",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10, color = 'black'),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 1)
    ggsave(p_CDKN2B_sum, filename = paste0("results/analyses/proximity_bias/", lib, "_bm_CDKN2B_sum.pdf"), width = 15, height = 15, dpi = 300)
}
