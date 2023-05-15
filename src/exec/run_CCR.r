library(CRISPRcleanR)
library(tidyverse)

# Genome-wide sorting of LFC
GW_sorting_LFC <- function(
  foldchanges,
  libraryAnnotation
) {
    libraryAnnotation <- libraryAnnotation %>%
        mutate(CHRM = replace(CHRM, CHRM == "X", 23)) %>%
        mutate(CHRM = replace(CHRM, CHRM == "Y", 24)) %>%
        mutate(CHRM = as.numeric(CHRM))

    if (ncol(foldchanges) > 3) {
        stop("The input data frame should have 3 columns: sgRNA, gene, and foldchanges")
    }else{
        converted <- foldchanges %>%
            left_join(libraryAnnotation %>%
                select(CHRM, GENES, STARTpos, ENDpos, seq), by = c('gene' = 'GENES', 'sgRNA' = 'seq')) %>%
            column_to_rownames('sgRNA') %>%
            select(CHR = CHRM, startp = STARTpos, endp = ENDpos, genes = gene, avgFC = colnames(.)[2]) %>%
            arrange(CHR, startp) %>% ## sort by chromosome and start position
            mutate(BP = startp + (endp - startp) / 2) %>%
            na.omit()
    }

    return(converted)
}

# main function
main <- function(raw_LFC_path, GuideMap_path, output_dir, lib){
    ## Load the data
    GuideMap <- read_csv(GuideMap_path) %>%
        filter(UsedByChronos == TRUE) %>%
        select(-c(nAlignments, DropReason, UsedByChronos)) %>%
        separate(col = "Gene", sep = " \\(", into = c("GENES", "code")) %>%
        select(-code) %>%
        na.omit() %>%
        separate(col = "GenomeAlignment", sep = '_', into = c("CHRM", "STARTpos", "STRAND")) %>%
        mutate(STARTpos = as.numeric(STARTpos)) %>%
        mutate(ENDpos = STARTpos + 23) %>% ## 23 is the length of the sgRNA
        separate(col = "CHRM", sep = "chr", into = c("null", "CHRM")) %>%
        select(-null) %>%
        mutate(CODE = sgRNA, seq = sgRNA) %>%
        select(-sgRNA) %>%
        na.omit() %>%
        distinct()

    raw_LFC <- read_csv(raw_LFC_path) %>%
        inner_join(GuideMap %>%
            select(GENES, CODE), by = c("sgRNA" = "CODE")) %>%
            select(sgRNA, GENES, everything()) %>%
            select(sgRNA, gene = GENES, everything())

    ## Genome-wide sorting and cleaning of LFC
    for (i in seq(3, ncol(raw_LFC))){
        CellLine <- colnames(raw_LFC)[i]

        gwSortedLFC <- GW_sorting_LFC(raw_LFC[, c(1:2, i)], GuideMap)
        correctedLFC <- ccr.GWclean(
            gwSortedLFC, 
            display = FALSE,
            verbose = 0)

        correctedLFC <- correctedLFC$corrected_logFCs %>%
            rownames_to_column('sgRNA') %>%
            select(sgRNA, genes, correctedFC) %>%
            rename(!!CellLine := correctedFC) ## rename corrected LFC to cell line name
        
        if (i == 3){
            correctedLFCs <- correctedLFC
        }else{
            correctedLFCs <- correctedLFCs %>%
                full_join(correctedLFC, by = c("sgRNA", "genes"))
        }

        print(paste0("Finished processing: ", CellLine, ". Progress: ", i-2, "/", ncol(raw_LFC)-2))
    }

    ## Gene-level LFC
    gene_correctedLFCs <- correctedLFCs %>%
        group_by(genes) %>%
        select(-sgRNA) %>%
        summarise_all(median, na.rm = TRUE) %>%
        ungroup()

    ## Save the corrected LFCs
    write_csv(correctedLFCs, paste0(output_dir, '/', lib, '_sgrna_CCR.csv'))
    write_csv(gene_correctedLFCs, paste0(output_dir, '/', lib, '_sgrna_CCR.csv'))
}

# Run the main function
## For instance:
## Rscript src/exec/run_CCR.r data/raw/Avana_sgrna_raw_LFC.csv data/AvanaGuideMap.csv data/corrected/ Avana
if (!interactive()){
    sys.args <- commandArgs(trailingOnly = TRUE)

    if (length(sys.args) < 4){
        stop(cat("Please provide these args in the following order: 
            \n- raw_LFC_path: path to the raw LFC file
            \n- GuideMap_path: path to the GuideMap file
            \n- output_dir: path to the output directory
            \n- lib: library name of the output files"))
    }

    raw_LFC_path <- sys.args[1]
    GuideMap_path <- sys.args[2]
    output_dir <- sys.args[3]
    lib <- sys.args[4]

    main(raw_LFC_path, GuideMap_path, output_dir, lib)
}