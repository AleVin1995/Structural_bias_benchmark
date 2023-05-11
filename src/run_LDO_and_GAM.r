library(tidyverse)
library(magrittr)
library(ggsci)

setwd("")

source("support_themes.R")

source("LDO_and_GAM_func.R")

outdir <- "Output"
figure.dir <- "Figures"
avana.data.dir <- "Avana"



cat("  Load data from Muñoz paper (Muñoz & al.(2016)) \n")
data <- readRDS("InputData.RDS")

cat("  Load data from second MET screen \n")
data_screen2 <- readRDS("InputData_screen2.RDS")

cat("  Load essential genes information from TableS3 in Wang & al.(2015) \n")
# Identification and characterization of essential genes in the human genome; DOI: 10.1126/science.aac7041
essential.genes <- read.csv("EssentialGenes.csv") %>%
  group_by(Gene) %>% do({
    xx <- .
    adj.p.val.vec <- c(xx$KBM7.adjusted.p.value,xx$K562.adjusted.p.value,xx$Jiyoye.adjusted.p.value,xx$Raji.adjusted.p.value)
    cs.vec <- c(xx$KBM7.CS, xx$K562.CS, xx$Jiyoye.CS, xx$Raji.CS)
    idx <- which.max(adj.p.val.vec)
    data.frame(., merge.CS = cs.vec[idx], merge.adjusted.p.value = adj.p.val.vec[idx])
  }) %>% as.data.frame %>%
  filter(merge.CS < -1, merge.adjusted.p.value < 0.05) %>% .$Gene

# save supp Table
write.csv(essential.genes, file = file.path( outdir, "suppTable.csv"), row.names = F)

BAGEL_essential <- readRDS("BAGEL_essential.RDS")
BAGEL_nonessential <- readRDS("BAGEL_nonEssential.RDS")


#####################################
####### GET SENSITIVITY SCORE #######
#####################################

# Calculate LogFC_org

cat(" -Calculate LogFC_org \n")
data <- data %>% group_by(CLEANNAME) %>% mutate(LogFC_org = logfc_org_func(SAMPLE_COUNT, REFERENCE_COUNT)) %>% ungroup



######################################
########### LDO CORRECTION ###########
######################################

# rename columns
data <- data %>%
  mutate( ESSENTIAL = GENESYMBOL %in% essential.genes) %>% 
  dplyr::rename(DEPENDENCY_SCORE = LogFC_org, GENE_NAME = GENESYMBOL, SAMPLE_NAME = CLEANNAME, CHROMOSOME = CHR)

# run Local Drop Out correction
data <- LDO( data,
             params_weigthed_mean = list(
               omega = 1e6,
               decay_func = (function(x)exp(x)),
               in_subunit_weight = 0,
               side = "both"
             ),
             params_regression = list(
               minbucket = 2, 
               minsubunits = 3, 
               cp_init = -3, 
               cp_iteration = 0.1),
             verbose = FALSE)



######################################
########### GAM CORRECTION ###########
######################################

# Smooth out multi-alignment and CNA effect: using generalized additive models
cat("  Smooth out CNA effect using generalized additive models \n")

data2 <- GAM(data, formula = "DEPENDENCY_SCORE ~ CNA", subunit = "SAMPLE_NAME", verbose = FALSE)



#######################################
################ PLOTS ################
#######################################


arranged.gs <- data %>% dplyr::filter(!is.na(CNA)) %>% arrange(CHROMOSOME,POSITION) %>% .$GENE_NAME %>% .[!duplicated(.)]
arranged.gs.tmp <- data %>% arrange(CHROMOSOME,POSITION) %>% .$GENE_NAME %>% .[!duplicated(.)]
data <- data %>% inset("GENESYMBOLS", value = factor(.$GENE_NAME, levels = arranged.gs.tmp))


###########
# figure 1a
tmp.plot1 <- ggplot( data %>% mutate(LogFC = DEPENDENCY_SCORE), aes( x = CNA, y = LogFC)) + geom_point(size=0.5) + facet_wrap(~SAMPLE_NAME) + geom_smooth(size=2) + 
  theme_Publication(base_family="arial", base_size = 8) +
  ggtitle("Copy Number effect before correction") + xlab("log2(CN)")
tmp.plot1


file.name <- "Figure1a.tiff"
ggsave(file.path(figure.dir, file.name), width = 3.5, height = 3.2, dpi=600)


###########
# figure 1b
pos.of.interest <- which(arranged.gs == "YAP1")
gs.of.interest <- arranged.gs[(pos.of.interest-6):(pos.of.interest+6)]
data.tmp <- data %>% filter( SAMPLE_NAME == "sf268", GENESYMBOLS %in% gs.of.interest) %>%
  transmute(SAMPLE_NAME, SEQ, Original = DEPENDENCY_SCORE, LDO_corrected = DEPENDENCY_SCORE_LDO) %>%
  gather(SMethod, LogFC, -SAMPLE_NAME, -SEQ) %>% merge(data, by = c("SAMPLE_NAME", "SEQ"), all.x = T) %>%
  mutate(SMethod = factor(SMethod, levels = c("Original", "LDO_corrected")))
tmp.plot2 <- ggplot( data.tmp, aes( x = GENESYMBOLS)) + 
  geom_boxplot(aes(y = LogFC), outlier.alpha = 0) + 
  geom_jitter(aes(y = LogFC), size = 0.5) + 
  geom_line(aes(y = -0.18*CNA, group = 1), colour = "red", alpha = 0.3, size = 3) +
  facet_wrap(~SMethod) + theme_Publication(base_family="arial", base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_hline(aes(yintercept=0), colour = "dodgerblue") +
  ggtitle("YAP1 amplicon in SF268") +
  scale_y_continuous(sec.axis = sec_axis(~-./0.18, name = "log2(CN)")) +
  theme(axis.title.y.right = element_text(colour = "indianred1"))
tmp.plot2

file.name <- "Figure1b.tiff"
ggsave(file.path(figure.dir, file.name), width = 6, height = 3.5, dpi=600)

###########
# figure 1c
tmp.plot3 <- ggplot( data %>% mutate(LogFC = DEPENDENCY_SCORE_LDO), aes( x = CNA, y = LogFC)) + geom_point(size=0.5) + facet_wrap(~SAMPLE_NAME) + geom_smooth(size = 2) + 
  theme_Publication(base_family="arial", base_size = 8) +
  ggtitle("Copy Number effect after LDO correction") + xlab("log2(CN)")
tmp.plot3

file.name <- "Figure1c.tiff"
ggsave(file.path(figure.dir, file.name), width = 3.5, height = 3.2, dpi=600)



###########
# figure 2a
pos.of.interest <- which(arranged.gs == "MET")
gs.of.interest <- arranged.gs[(pos.of.interest-6):(pos.of.interest+6)]
data.tmp <- data %>% filter( SAMPLE_NAME == "mkn45", GENESYMBOLS %in% gs.of.interest) %>%
  transmute(SAMPLE_NAME, SEQ, Original = DEPENDENCY_SCORE, LDO_corrected = DEPENDENCY_SCORE_LDO) %>%
  gather(SMethod, LogFC, -SAMPLE_NAME, -SEQ) %>% merge(data, by = c("SAMPLE_NAME", "SEQ"), all.x = T) %>%
  mutate(SMethod = factor(SMethod, levels = c("Original", "LDO_corrected")))
tmp.plot <- ggplot( data.tmp %>% dplyr::filter(SMethod %in% "Original"), aes( x = GENESYMBOLS)) + 
  geom_boxplot(aes(y = LogFC), outlier.alpha = 0) + 
  geom_jitter(aes(y = LogFC), size = 0.5) + 
  geom_line(aes(y = -0.18*CNA, group = 1), colour = "red", alpha = 0.3, size = 3) +
  # facet_wrap(~SMethod) + 
  theme_Publication(base_family="arial", base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_hline(aes(yintercept=0), colour = "dodgerblue") +
  ggtitle("MET amplicon in MKN45")+
  scale_y_continuous(sec.axis = sec_axis(~-./0.18, name = "log2(CN)")) +
  theme(axis.title.y.right = element_text(colour = "indianred1"))
tmp.plot

file.name <- "Figure2a.tiff"
ggsave(file.path(figure.dir, file.name), width = 3.5, height = 3, dpi=600)

###########
# figure 2b
tmp.plot <- ggplot(data_screen2 %>% mutate(LogFC = LogFC_org), aes( x = GENESYMBOLS, y = LogFC)) + geom_boxplot(aes(), outlier.alpha = 0) +
  theme_Publication(base_family="arial", base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_hline(aes(yintercept=0), colour = "dodgerblue") + geom_jitter(size = 0.5) + 
  geom_line(aes(y = -0.17*CNA, group = 1, colour = "log2(CNA)", alpha = 0.5), size = 3) + 
  guides(size = FALSE, alpha = FALSE) + scale_colour_discrete(name = "Labels") +
  ggtitle("MET amplicon in MKN45 [Screen 2]") + theme(legend.position = "none") + ylim(c(-1.9,0.6)) +
  scale_y_continuous(sec.axis = sec_axis(~-./0.17, name = "log2(CN)")) +
  theme(axis.title.y.right = element_text(colour = "indianred1"))
tmp.plot

file.name <- "Figure2b.tiff"
ggsave(file.path(figure.dir, file.name), width = 3.5, height = 3, dpi=600)


###########
# figure 3a

indir <- avana.data.dir
filename <- "genelevel_ldo.RDS"
data.gl <- readRDS( file.path(indir, filename))

data.gl[ which(data.gl$CN >= 10), "CN"] <- "10+"
data.gl[ is.na(data.gl$CN), "CN"] <- "NA"
data.gl$CN <- factor( data.gl$CN, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10+", "NA"))


data.tmp <- data.gl %>% 
  dplyr::filter( !(CN %in% "NA")) %>%
  transmute(Replicate, Gene.Symbol, CN, UncorrectedData = scaledScore, LDO = scaledScore_LDO) %>%
  gather(ScoreType, Score, -Replicate, -Gene.Symbol, -CN) %>%
  mutate(ScoreType = factor(ScoreType, levels = c("UncorrectedData","LDO")))

tmp.plot <- ggplot( data.tmp, aes( x = CN)) + 
  theme_Publication(base_family="arial", base_size = 10) +
  geom_boxplot(aes( y = Score, group = CN, colour = ScoreType), outlier.shape = NA) +
  geom_hline( yintercept = 0, linetype = "longdash") +
  geom_hline( yintercept = -1, linetype = "longdash") +
  ylim( c(-2, 1)) +
  xlab("Copy Number") +
  facet_grid(~ScoreType) + 
  guides( colour = FALSE)
tmp.plot

file.name <- "Figure3a.tiff"
ggsave(file.path(figure.dir, file.name), width = 12, height = 10, units = "cm", dpi=600)


###########
# figure 3b

indir <- avana.data.dir
filename <- "genelevel_ldo.RDS"
data.tmp <- readRDS( file.path( indir, filename)) %>% dplyr::filter( Replicate %in% "DAN.G.311Cas9.Rep.A.p6")

Score_recall_1 <- data.tmp %>%
  mutate(Essential = ifelse(Gene.Symbol %in% BAGEL_essential, "pos", "rest"),
         nonEssential = ifelse(Gene.Symbol %in% BAGEL_nonessential, "pos", "rest"),
         Amplified = ifelse(CN < 6 | is.na(CN), "neg", "pos")) %>%
  arrange(Replicate, Score) %>%
  group_by(Replicate) %>%
  mutate(Essential_pos_Sum = cumsum(Essential == "pos"),
         nonEssential_pos_Sum = cumsum(nonEssential == "pos"),
         Amplified_pos_Sum = cumsum(Amplified == "pos")) %>%
  mutate(Essential_Recall = Essential_pos_Sum / sum(Essential == "pos"),
         nonEssential_Recall = nonEssential_pos_Sum / sum(nonEssential == "pos"),
         Amplified_Recall = Amplified_pos_Sum / sum(Amplified == "pos")) %>%
  mutate(Dataset = "UncorrectedData", Score_percentile = 100*(1:length(Score))/length(Score), CN_burden = sum(CN > 5, na.rm = T)) %>%
  select(Replicate, Essential_Recall, nonEssential_Recall, Amplified_Recall, Dataset, Score_percentile, CN_burden) %>%
  ungroup %>%
  gather(Test, Recall, -Replicate, -Dataset, -Score_percentile, -CN_burden)

LDO_recall_1 <- data.tmp %>%
  mutate(Essential = ifelse(Gene.Symbol %in% BAGEL_essential, "pos", "rest"),
         nonEssential = ifelse(Gene.Symbol %in% BAGEL_nonessential, "pos", "rest"),
         Amplified = ifelse(CN < 6 | is.na(CN), "neg", "pos")) %>%
  arrange(Replicate, Score_LDO) %>%
  group_by(Replicate) %>%
  mutate(Essential_pos_Sum = cumsum(Essential == "pos"),
         nonEssential_pos_Sum = cumsum(nonEssential == "pos"),
         Amplified_pos_Sum = cumsum(Amplified == "pos")) %>%
  mutate(Essential_Recall = Essential_pos_Sum / sum(Essential == "pos"),
         nonEssential_Recall = nonEssential_pos_Sum / sum(nonEssential == "pos"),
         Amplified_Recall = Amplified_pos_Sum / sum(Amplified == "pos")) %>%
  mutate(Dataset = "LDO", Score_percentile = 100*(1:length(Score_LDO))/length(Score_LDO), CN_burden = sum(CN > 5, na.rm = T)) %>%
  select(Replicate, Essential_Recall, nonEssential_Recall, Amplified_Recall, Dataset, Score_percentile, CN_burden) %>%
  ungroup %>%
  gather(Test, Recall, -Replicate, -Dataset, -Score_percentile, -CN_burden)

merge_recall <- LDO_recall_1 %>% rbind(Score_recall_1)
merge_recall$Dataset <- factor(merge_recall$Dataset, levels = c("UncorrectedData", "LDO"))

tmp.plot <- ggplot( merge_recall %>% dplyr::filter(Replicate %in% "DAN.G.311Cas9.Rep.A.p6") %>% mutate(Test = Test %>% gsub("_Recall","",.)), aes(x = Score_percentile, y = Recall, colour = Test)) +
  theme_Publication(base_family="arial", base_size = 10) +
  geom_line() + facet_grid(Dataset~.) + geom_abline(intercept = 0, slope = 0.01, alpha = 0.5) +
  ggtitle("DAN-G pancreas") +
  theme(legend.title=element_blank())
tmp.plot

file.name <- "Figure3b.tiff"
ggsave(file.path(figure.dir, file.name), width = 6, height = 10, units = "cm", dpi=600)

###########
# figure 3c

indir <- avana.data.dir
filename <- "auc_ldo.RDS"

data3c <- readRDS( file.path( indir, filename))

tmp.plot <- ggplot( data3c %>% spread(Dataset, AUC) %>% mutate(CN_burden = ifelse(CN_burden > 25, "High", "Low") %>% factor(levels = c("Low","High"))), aes( x = UncorrectedData, y = LDO, colour = Test)) + 
  theme_Publication(base_family="arial", base_size = 10) +
  geom_point(aes(size = CN_burden), alpha = 0.5) + 
  facet_wrap(~Test) + 
  geom_abline(slope = 1, intercept = 0, alpha = 0.5) + geom_hline(yintercept = 50, linetype = "longdash", alpha = 0.5) + geom_vline(xintercept = 50, linetype = "longdash", alpha = 0.5) +
  ggtitle("Area under the recall curve") +
  guides(colour=FALSE, size=guide_legend(title = "Copy number burden")) +
  scale_size_discrete(range = c(1,5))
tmp.plot

file.name <- "Figure3c.tiff"
ggsave(file.path(figure.dir, file.name), width = 20, height = 10, units = "cm", dpi=600)



####################
### SUPP FIGURES ###
####################

###########
# supp Figure1

# MERGE SCREEN 1 AND SCREEN 2
jnk.amp.mkn1 <- data %>% dplyr::filter(SAMPLE_NAME %in% "mkn45") %>% transmute( SAMPLE_NAME, GENESYMBOLS, SEQ, POSITION, DEPENDENCY_SCORE, SCREEN = "Screen 1")
jnk.amp.mkn2 <- data_screen2 %>% transmute( SAMPLE_NAME = CLEANNAME, GENESYMBOLS, SEQ, POSITION, DEPENDENCY_SCORE = LogFC_org, SCREEN = "Screen 2")
screen.mkn45 <- jnk.amp.mkn1 %>% rbind( jnk.amp.mkn2)
COMMON.SEQ <- screen.mkn45$SEQ %>% table %>% .[is_greater_than(.,1)] %>% names
screen.mkn45 <- screen.mkn45 %>% mutate( COMMON = SEQ %in% COMMON.SEQ) 


tmp.plot <- ggplot( screen.mkn45 %>% dplyr::filter(GENESYMBOLS %in% "MET"), aes( x = POSITION, y = DEPENDENCY_SCORE)) + geom_point() + facet_wrap(~SCREEN) +
  theme_Publication(base_family="arial") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
tmp.plot

file.name <- "Supp_Figure1.tiff"
ggsave(file.path(figure.dir, file.name), width = 7, height = 5, dpi = 600)

###########
# supp Figure2
tmp.plot <- ggplot( data2 %>% mutate(LogFC = DEPENDENCY_SCORE_GAM), aes( x = CNA, y = LogFC)) + geom_point(size=0.5) + facet_wrap(~SAMPLE_NAME) + geom_smooth(size = 2) + 
  theme_Publication(base_family="arial", base_size = 8) +
  ggtitle("Copy Number effect after GAM correction") + xlab("log2(CN)")
tmp.plot

file.name <- "Supp_Figure2.tiff"
ggsave(file.path(figure.dir, file.name), width = 3.5, height = 3.2, dpi=600)


###########
# supp Figure3
data2 <- data %>% mutate(colour = "black")
idx.orange <- (data2$CNA > 1.8 & data2$LDO_fix > -0.2) %in% TRUE
idx.red    <- (data2$CNA < 2 & data2$LDO_fix < -0.25) %in% TRUE
data2[idx.orange,"colour"] <- "orange"
data2[idx.red,   "colour"]    <- "red"

tmp.plot <- ggplot( data2, aes( x = CNA, y = LDO_fix, colour = colour)) + geom_point() + facet_wrap(~SAMPLE_NAME) + 
  theme_Publication(base_family="arial") + scale_colour_manual(values = c("black", "orange", "red")) + guides(colour = FALSE) +
  ggtitle("Copy Number vs LDO correction") +
  xlab("log2(CN)")
tmp.plot

file.name <- "Supp_Figure3.tiff"
ggsave(file.path(figure.dir, file.name), width = 7, height = 5, dpi = 600)


###########
# supp Figure4
###########
# Supp figure 4a

indir <- avana.data.dir
filename <- "genelevel_gam.RDS"

data.gl <- readRDS( file.path(indir, filename))
data.gl[ which(data.gl$CN >= 10), "CN"] <- "10+"
data.gl[ is.na(data.gl$CN), "CN"] <- "NA"
data.gl$CN <- factor( data.gl$CN, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10+", "NA"))


data.tmp <- data.gl %>% 
  dplyr::filter( !(CN %in% "NA")) %>%
  transmute(Replicate, Gene.Symbol, CN, UncorrectedData = scaledScore, GAM = scaledScore_GAM) %>%
  gather(ScoreType, Score, -Replicate, -Gene.Symbol, -CN) %>%
  mutate(ScoreType = factor(ScoreType, levels = c("UncorrectedData","GAM")))

tmp.plot <- ggplot( data.tmp, aes( x = CN)) + 
  theme_Publication(base_family="arial", base_size = 10) +
  geom_boxplot(aes( y = Score, group = CN, colour = ScoreType), outlier.shape = NA) +
  geom_hline( yintercept = 0, linetype = "longdash") +
  geom_hline( yintercept = -1, linetype = "longdash") +
  ylim( c(-2, 1)) +
  xlab("Copy Number") +
  facet_grid(~ScoreType) + 
  guides( colour = FALSE)
tmp.plot

file.name <- "Supp_Figure4a.tiff"
ggsave(file.path(figure.dir, file.name), width = 10, height = 10, units = "cm", dpi = 600)


###########
# Supp figure 4b

indir <- avana.data.dir
filename <- "genelevel_gam.RDS"
data.tmp <- readRDS( file.path(indir, filename)) %>% dplyr::filter( Replicate %in% "DAN.G.311Cas9.Rep.A.p6")

Score_recall_1 <- data.tmp %>%
  mutate(Essential = ifelse(Gene.Symbol %in% BAGEL_essential, "pos", "rest"),
         nonEssential = ifelse(Gene.Symbol %in% BAGEL_nonessential, "pos", "rest"),
         Amplified = ifelse(CN < 6 | is.na(CN), "neg", "pos")) %>%
  arrange(Replicate, Score) %>%
  group_by(Replicate) %>%
  mutate(Essential_pos_Sum = cumsum(Essential == "pos"),
         nonEssential_pos_Sum = cumsum(nonEssential == "pos"),
         Amplified_pos_Sum = cumsum(Amplified == "pos")) %>%
  mutate(Essential_Recall = Essential_pos_Sum / sum(Essential == "pos"),
         nonEssential_Recall = nonEssential_pos_Sum / sum(nonEssential == "pos"),
         Amplified_Recall = Amplified_pos_Sum / sum(Amplified == "pos")) %>%
  mutate(Dataset = "UncorrectedData", Score_percentile = 100*(1:length(Score))/length(Score), CN_burden = sum(CN > 5, na.rm = T)) %>%
  select(Replicate, Essential_Recall, nonEssential_Recall, Amplified_Recall, Dataset, Score_percentile, CN_burden) %>%
  ungroup %>%
  gather(Test, Recall, -Replicate, -Dataset, -Score_percentile, -CN_burden)

GAM_recall_1 <- data.tmp %>%
  mutate(Essential = ifelse(Gene.Symbol %in% BAGEL_essential, "pos", "rest"),
         nonEssential = ifelse(Gene.Symbol %in% BAGEL_nonessential, "pos", "rest"),
         Amplified = ifelse(CN < 6 | is.na(CN), "neg", "pos")) %>%
  arrange(Replicate, Score_GAM) %>%
  group_by(Replicate) %>%
  mutate(Essential_pos_Sum = cumsum(Essential == "pos"),
         nonEssential_pos_Sum = cumsum(nonEssential == "pos"),
         Amplified_pos_Sum = cumsum(Amplified == "pos")) %>%
  mutate(Essential_Recall = Essential_pos_Sum / sum(Essential == "pos"),
         nonEssential_Recall = nonEssential_pos_Sum / sum(nonEssential == "pos"),
         Amplified_Recall = Amplified_pos_Sum / sum(Amplified == "pos")) %>%
  mutate(Dataset = "GAM", Score_percentile = 100*(1:length(Score_GAM))/length(Score_GAM), CN_burden = sum(CN > 5, na.rm = T)) %>%
  select(Replicate, Essential_Recall, nonEssential_Recall, Amplified_Recall, Dataset, Score_percentile, CN_burden) %>%
  ungroup %>%
  gather(Test, Recall, -Replicate, -Dataset, -Score_percentile, -CN_burden)

merge_recall <- GAM_recall_1 %>% rbind(Score_recall_1)
merge_recall$Dataset <- factor(merge_recall$Dataset, levels = c("UncorrectedData", "GAM"))

tmp.plot <- ggplot( merge_recall %>% dplyr::filter(Replicate %in% "DAN.G.311Cas9.Rep.A.p6") %>% mutate(Test = Test %>% gsub("_Recall","",.)), aes(x = Score_percentile, y = Recall, colour = Test)) +
  theme_Publication(base_family="arial", base_size = 10) +
  geom_line() + facet_grid(Dataset~.) + geom_abline(intercept = 0, slope = 0.01, alpha = 0.5) +
  ggtitle("DAN-G pancreas") +
  theme(legend.title=element_blank())
tmp.plot

file.name <- "Supp_Figure4b.tiff"
ggsave(file.path(figure.dir, file.name), width = 6, height = 10, units = "cm", dpi = 600)

###########
# Supp figure 4c

indir <- avana.data.dir
filename <- "auc_gam.RDS"
data <- readRDS( file.path(indir, filename))

tmp.plot <- ggplot( data %>% spread(Dataset, AUC) %>% mutate(CN_burden = ifelse(CN_burden > 25, "High", "Low") %>% factor(levels = c("Low","High"))), aes( x = UncorrectedData, y = GAM, colour = Test)) + 
  theme_Publication(base_family="arial", base_size = 10) +
  geom_point(aes(size = CN_burden), alpha = 0.5) + 
  facet_wrap(~Test) + 
  geom_abline(slope = 1, intercept = 0, alpha = 0.5) + geom_hline(yintercept = 50, linetype = "longdash", alpha = 0.5) + geom_vline(xintercept = 50, linetype = "longdash", alpha = 0.5) +
  ggtitle("Area under the recall curve") +
  guides(colour=FALSE, size=guide_legend(title = "Copy number burden")) +
  scale_size_discrete(range = c(1,5))
tmp.plot

file.name <- "Supp_Figure4c.tiff"
ggsave(file.path(figure.dir, file.name), width = 20, height = 10, units = "cm", dpi = 600)


