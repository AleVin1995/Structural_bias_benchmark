library(tidyverse)

libs <- c("Avana", "KY")
nums <- c(5, 10, 20, 50)

# iterate over algorithms and libraries
for (lib in libs){
    dfs_names <- paste0("MAGeCK_", nums)

    dfs <- paste0("data/corrected/", lib, "_gene_MAGeCK_", nums, ".csv") %>%
        map(~.x %>%
            read_csv %>%
            ## fill na with 0
            mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
            pivot_longer(-1, names_to = "ModelID", values_to = "LFC")) %>%
        reduce(., full_join, by = c("Gene", "ModelID")) %>%
        ## set column names
        set_names(c("Gene", "ModelID", dfs_names))
    
    ## compute pointwise difference between each dataset and the next one
    dfs <- dfs %>%
        mutate(`10-5` = `MAGeCK_10` - `MAGeCK_5`,
            `20-10` = `MAGeCK_20` - `MAGeCK_10`,
            `50-20` = `MAGeCK_50` - `MAGeCK_20`)

    ## compute pointwise mean across iterations
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

    ## compute pointwise standard deviation across iterations
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
    saveRDS(dfs, paste0("results/analyses/mageck_batching/", lib, "_MAGeCK_diff.rds"))
    saveRDS(MAGeCK_sd, paste0("results/analyses/mageck_batching/", lib, "_MAGeCK_sd.rds"))
}
