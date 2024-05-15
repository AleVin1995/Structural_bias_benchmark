library(extrafont)
library(forcats)
library(patchwork)
library(RColorBrewer)
library(tidyverse)

font_import(paths = "arial", prompt = FALSE)

libs <- c("Avana", "KY")
cols <- c("#B3B3B3", brewer.pal(n = 7, name = "Dark2"))
cols <- c(cols[1:3], "#EE8208", cols[4:8])

# iterate over algorithms and libraries
for (lib in libs){
    ## pooled results
    bm_pool <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_pool.rds"))
    bm_pool$Algorithm <- factor(bm_pool$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    bm_pool$Coord <- factor(bm_pool$Coord, levels = c(paste0(rep(c(1:23, "X", "Y"), each = 2), c("p", "q"))))

    ## perform a systematic t-test between Uncorrected vs other algorithms
    # bm_pool %>%
    #     filter(Algorithm != "Uncorrected") %>%
    #     inner_join(bm_pool %>% 
    #         filter(Algorithm == "Uncorrected") %>%
    #         select(-Algorithm), by = c("Coord")) %>%
    #     group_by(Algorithm) %>%
    #     summarise(p = t.test(est.x, est.y)$p.value)


    p_pool_sum <- ggplot(bm_pool, aes(x = Algorithm, y = est, fill = Algorithm)) +
        geom_boxplot() +
        geom_jitter(width = 0.1, size = 2) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "P(intra-arm cosine > inter)") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1),
            axis.ticks.length = unit(0.5, "cm"),
            axis.title = element_text(size = 30, color = 'black'),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio = 1,
            text = element_text(family = "Arial"),
            plot.margin = unit(c(2,1,1,1), "cm"),
            legend.position = "none") +
        scale_fill_manual(values = cols)
        
    ## TP53 results
    bm_TP53 <- readRDS(paste0("results/analyses/proximity_bias/", lib, "_bm_TP53.rds"))
    bm_TP53$Algorithm <- factor(bm_TP53$Algorithm, levels = c("Uncorrected", "CCR", "Chronos", "AC Chronos", "Crispy", "GAM", "Geometric", "LDO", "MAGeCK"))
    bm_TP53$Coord <- factor(bm_TP53$Coord, levels = c(paste0(rep(c(1:23, "X", "Y"), each = 2), c("p", "q"))))

    ## perform a systematic t-test between Uncorrected vs other algorithms
    # bm_TP53 %>%
    #     filter(Algorithm != "Uncorrected" & Status == "WT") %>%
    #     inner_join(bm_TP53 %>% 
    #         filter(Algorithm == "Uncorrected" & Status == "WT") %>%
    #         select(-Algorithm), by = c("Coord")) %>%
    #     group_by(Algorithm) %>%
    #     summarise(p = t.test(est.x, est.y)$p.value)
    
    # bm_TP53 %>%
    #     filter(Algorithm != "Uncorrected" & Status == "Mut") %>%
    #     inner_join(bm_TP53 %>% 
    #         filter(Algorithm == "Uncorrected" & Status == "Mut") %>%
    #         select(-Algorithm), by = c("Coord")) %>%
    #     group_by(Algorithm) %>%
    #     summarise(p = t.test(est.x, est.y)$p.value)

    p_TP53_sum <- ggplot(bm_TP53, aes(x = Algorithm, y = est, fill = Status)) +
        geom_boxplot(outlier.shape=NA) +
        geom_point(position=position_jitterdodge(jitter.width = 0.1)) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        labs(x = "", y = "") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1),
            axis.ticks.length = unit(0.5, "cm"),
            axis.title = element_text(size = 30, color = 'black'),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio = 1,
            text = element_text(family = "Arial"),
            plot.margin = unit(c(2,1,1,1), "cm"),
            legend.key.size = unit(1.5, 'cm'),
            legend.text = element_text(size = 20),
            legend.position = c(0.15, 0.95),
            legend.background = element_rect(fill = NA, color = NA)) +
        scale_fill_manual(labels = c("TP53 mut", "TP53 wt"), 
            values = c("#1F78B4", "#A6CEE3"), name = "")
    
    ## Assemble panel
    panel <- p_pool_sum + p_TP53_sum +
        plot_annotation(tag_levels = 'A') &
        theme(plot.tag.position = c(0, 1),
            plot.tag = element_text(size = 40, face = "bold", family = "Arial"))
    ggsave(panel, filename = paste0("results/panels/proximity_bias/proximity_bias_", lib, "_bm_all.pdf"), width = 30, height = 15, dpi = 300)
}
