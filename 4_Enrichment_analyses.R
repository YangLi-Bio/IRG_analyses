##################################################################################
#                                                                                #
# KEGG pathway enrichment analyses for selected eGRNs, including:                #
# 1. Barplots to depict different enriched pathway                               #
# 2. heatmaps to showcase the pathway genes across different cell subpopulation  #
# 3. Others?                                                                     #
#                                                                                #
##################################################################################


# Research plan:
# 1. Heatmap : rows - eGRN genes associated to various phenotypes,
#    columns - cells divided into different time points or cell types, 
#    elements - mean expression values. If the adjusted p-value is insignificant, 
#    draw the negative logarithm p-value instead of adjusted p-values.
# 2. Barplots : negative logarithm FDR against KEGG, Hallmark, and Reactome


# Libraries
library(enrichR)
library(pbapply)
library(qs)


# Global parameters
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Rfile/"
table.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Tables/"
IgHV.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/1_classification_GSE113386/"
TP53.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/3_TP53_GSE11038/"
cox.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/2_survival_GSE22762/"


# Pathway enrichment analyses for MYCN
MYCN.df <- read.table(paste0(table.dir, "MYCN_nodes.txt"), header = T)
dim(MYCN.df)
head(MYCN.df)
# MYCN.enriched <- qs::qread(paste0(R.dir, "MYCN_enriched.qsave"))
# MYCN.pathways <- MYCN.enriched$KEGG_2019_Human
# dim(MYCN.pathways)
# MYCN.pathways <- MYCN.pathways[MYCN.pathways$Adjusted.P.value < 0.05,]
# dim(MYCN.pathways)
# colnames(MYCN.pathways)
# head(MYCN.pathways)
# geneSet.list <- Reduce("union", pbapply(MYCN.pathways, 1, function(x) {
#   strsplit(x[10], split = ";")[[1]]
# }))
