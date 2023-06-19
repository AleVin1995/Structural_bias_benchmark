library(ceres)
library(magrittr)
library(tidyverse)

# For each copy number segment in a cell line, find the sgRNAs that intersect it.
intersect_guide_cn <- function(guide, cn_seg, verbose=FALSE){
    guide_cn <- guide %>%
        mutate(segments = NA)

    for (i in 1:nrow(guide_cn)){
        guide_cn$segments[i] <- list(cn_seg %>% 
            filter(Chromosome == guide_cn$Chromosome[i] & Start <= guide_cn$Start[i] & End >= guide_cn$Start[i]) %>% 
            select(ModelID, SegmentMean) %>%
            group_by(ModelID) %>%
            summarize(SegmentMean = mean(SegmentMean)) %>% ## some segments overlap multiple sgRNAs
            ungroup())
        
        if (verbose & i %% 1000 == 0){
            print(paste0("guide ", i, "/", nrow(guide_cn)))
        }
    }

    guide_cn <- guide_cn %>% 
        hoist(segments, ModelID = 1, SegmentMean = 2) %>%
        unnest_longer(col = c(ModelID, SegmentMean)) %>%
        unite("Locus", c(Chromosome, Start, Strand), sep = "_")
    
    return(guide_cn)
}

# Wrapper to run CERES
wrap_ceres <- function(sg, cn, guide_locus,
                       locus_gene, replicate_map, 
                       run_id = "CERES",
                       fit_efficacy=TRUE){
    
    ## normalize sgRNA LFCs
    sg_data <- sg %>%
        as.matrix() %>%
        apply(., 2, function(cl){
            (cl - median(cl, na.rm = TRUE)) / mad(cl, na.rm = TRUE)
        })
        
    cn_data <- cn %>%
        as.matrix()

    ## get common cell lines
    common_cell_lines <- intersect(colnames(sg_data), colnames(cn_data))

    sg_data <- sg_data[, common_cell_lines]
    cn_data <- cn_data[, common_cell_lines]
    replicate_map <- replicate_map %>%
        filter(CellLine %in% common_cell_lines)

    ## set default parameters
    params <- list(lambda_g=0.68129207,
        lambda_o=0.01,
        lambda_s=0.001,
        n_segments=25,
        validation_set=0,
        run_name=run_id)
    
    res <- run_ceres(sg_data, cn_data, guide_locus, locus_gene,
                      replicate_map, params, fit_efficacy)
  
    return(res)
}

# main function
main <- function(raw_LFC_path, GuideMap_path, CN_seg_path, 
    profile_path, output_dir, lib){
    
    ## Load the data
    sg <- read_csv(raw_LFC_path) %>%
        column_to_rownames('sgRNA')

    guide_locus <- read_csv(GuideMap_path) %>%
        filter(UsedByChronos == TRUE) %>%
        separate(col = "Gene", sep = " \\(", into = c("Gene", "code")) %>%
        rename(Guide = sgRNA, Locus = GenomeAlignment, Value = nAlignments) %>%
        select(Guide, Locus, Gene, Value)
    locus_gene <- guide_locus %>%
        select(Locus, Gene, Value)
    guide_locus <- guide_locus %>%
        select(Guide, Locus, Value)

    replicate_map <- colnames(sg) %>%
        as_tibble() %>%
        rename(Replicate = value) %>%
        mutate(CellLine = Replicate)

    ## Intersect copy number segments with sgRNAs
    guide_cn <- read_csv(CN_seg_path) %>%
        full_join(read_csv(profile_path), by = "ProfileID") %>%
        select(ModelID, Chromosome, Start, End, SegmentMean) %>%
        mutate(SegmentMean = log2(SegmentMean)) %>%
        mutate(Chromosome = paste0("chr", Chromosome)) %>%
        intersect_guide_cn(guide_locus %>% 
            separate(col = "Locus", sep = "_", into = c("Chromosome", "Start", "Strand")) %>%
            mutate(Start = as.numeric(Start)) %>%
            select(-Value), ., verbose=TRUE) %>%
        select(-Guide) %>%
        pivot_wider(names_from = ModelID, values_from = SegmentMean) %>%
        column_to_rownames("Locus")

    ## Run CERES
    res <- wrap_ceres(sg=sg, cn=guide_cn, guide_locus=guide_locus,
        locus_gene=locus_gene, replicate_map=replicate_map)
    
    gene_correctedLFCs <- res$gene_essentiality_results$ge_fit %>%
        as.data.frame() %>%
        rownames_to_column("Gene")

    ## Save results
    write_csv(gene_correctedLFCs, paste0(output_dir, '/', paste0(lib, "_gene_CERES.csv")))
}

# Run the main function
## For instance:
## Rscript src/method/run_CERES.r data/raw/Avana_sgrna_raw_LFC.csv data/AvanaGuideMap.csv \
##      data/OmicsCNSegmentsProfile.csv data/OmicsProfiles.csv data/corrected/ Avana
if (!interactive()){
    sys.args <- commandArgs(trailingOnly = TRUE)

    if (length(sys.args) < 6){
        stop(cat("Please provide these args in the following order: 
            \n- raw_LFC_path: path to the raw LFC file
            \n- GuideMap_path: path to the GuideMap file
            \n- CN_seg_path: path to the copy number segments file
            \n- profile_path: path to the copy number profiles file
            \n- output_dir: path to the output directory
            \n- lib: library name of the output files"))
    }

    raw_LFC_path <- sys.args[1]
    GuideMap_path <- sys.args[2]
    CN_seg_path <- sys.args[3]
    profile_path <- sys.args[4]
    output_dir <- sys.args[5]
    lib <- sys.args[6]

    main(raw_LFC_path, GuideMap_path, CN_seg_path, profile_path, output_dir, lib)
}