library(extrafont)
library(forcats)
library(patchwork)
library(RColorBrewer)
library(tidyverse)

font_import(paths = "/group/iorio/Alessandro/CN_benchmark/arial", prompt = FALSE)

libs <- c("Avana", "KY")
projects <- c("Project Achilles", "Project Score")
cols <- c("#B3B3B3", brewer.pal(n = 7, name = "Dark2"))

panel <- list()

# iterate over algorithms and libraries
for (i in 1:2){
    lib <- libs[i]
    project <- projects[i]

    ## MAGeCK pointwise difference
    MAGeCK_diff <- readRDS(paste0("results/analyses/mageck_batching/", lib, "_MAGeCK_diff.rds")) %>%
        select(`10-5`, `20-10`, `50-20`) %>%
        pivot_longer(everything(.), names_to = "Type", values_to = "Diff")

    p_diff <- ggplot(MAGeCK_diff, aes(y = Diff, x = Type)) +
        geom_violin() +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(x = "Comparison", title = project) +
        {if (i == 1) labs(y = expression(Delta*"Beta score")) else labs(y = "")} +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 25, color = 'black'),
            axis.title = element_text(size = 28),
            plot.title = element_text(size = 32, hjust = 0.5),
            aspect.ratio = 1,
            plot.margin = grid::unit(c(1,1,1,1), "cm"))
    
    ## MAGeCK pointwise binning of the difference
    p_bin <- MAGeCK_diff %>%
        filter(Type == "50-20") %>%
        mutate(Cut = cut(Diff, breaks = c(-5, seq(-0.5, 0.5, 0.1), 5))) %>%
        group_by(Cut) %>%
        mutate(Perc = 100 * (n() / nrow(.))) %>%
        ungroup() %>%
        select(Cut, Perc) %>%
        distinct() %>%
        ggplot(aes(x = Cut, y = Perc)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        labs(x = "Bin") +
        {if (i == 1) labs(y = "Percentage (%)") else labs(y = "")} +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 25, color = 'black'),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_text(size = 28),
            aspect.ratio = 1,
            plot.margin = grid::unit(c(1,1,1,1), "cm"))

    ## MAGeCK pointwise sd
    MAGeCK_sd <- readRDS(paste0("results/analyses/mageck_batching/", lib, "_MAGeCK_sd.rds"))

    p_sd <- ggplot(MAGeCK_sd, aes(y = log10(sd), x = factor(0))) +
        geom_boxplot() +
        labs(x = "") +
        {if (i == 1) labs(y = "log10(sd)") else labs(y = "")} +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 25, color = 'black'),
            axis.title = element_text(size = 28),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            aspect.ratio = 1,
            plot.margin = grid::unit(c(1,1,1,1), "cm"))
    
    ## Assemble plots
    panel[[i]] <- p_diff + p_bin + p_sd +
        plot_layout(nrow = 3)
}

# Save results
panel <- wrap_plots(panel[[1]], panel[[2]]) +
    plot_annotation(tag_levels = list(c("A", "C", "E", "B", "D", "F"))) &
    theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 40, face = "bold", family = "Arial"))
ggsave(panel, filename = "results/panels/suppl_panel_mageck.pdf", width = 16, height = 24, dpi = 300)
