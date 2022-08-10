####################################################
#                                                  #
#      Tune parameters for three applications      #
#                                                  #
####################################################


# Libraries
library(Scissor)
library(qs)
library(Seurat)
library(dplyr)
library(data.table)
library(pbmcapply)
library(parallel)


# Global parameters
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Codes/"
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Rfile/Scissor_tuned/"
dir.create(R.dir) 
ighv.file <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/1_classification_GSE113386/bulk.qsave"


# Functions
hyper_test_enhancer_grns <- function(enhancer.grns, cells, total) {
  
  dt <- rbindlist(pbmclapply(enhancer.grns, function(x) {
    overlap <- length(intersect(x$cells, cells))
    group1 <- length(x$cells)
    group2 <- length(cells)
    pval <- phyper(overlap - 1, group2, total - group2, 
                   group1, lower.tail = FALSE)
    list(overlap, group1, group2, total, pval)
  }, mc.cores = detectCores()))
  colnames(dt) <- c("Overlap", "Group1", "Group2", "Total", "Pvalue")
  dt
}


# Load scRNA-seq data
sc_dataset <- readRDS("/fs/ess/PCON0022/liyang/STREAM/benchmarking/10x_lymph_node_lymphoma_14k.RDS") %>% 
  GetAssayData(slot = "data", assay = "RNA") # get the expression matrix
dim(sc_dataset) # 36601 14104
class(sc_dataset)
sc.obj <- Seurat_preprocessing(sc_dataset, verbose = F)
class(sc.obj)
qs::qsave(sc.obj, paste0(R.dir, "Scissor_processed_sc_obj.qsave"))


# Load the eGRNs
enhancer.grns <- qs::qread("/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Rfile/Extended_eGRNs.qsave")
names(enhancer.grns[[1]])
lapply(enhancer.grns, "[[", "genes") %>% sapply(., length)
lapply(enhancer.grns, "[[", "peaks") %>% sapply(., length)
lapply(enhancer.grns, "[[", "cells") %>% sapply(., length)
enhancer.grn.cells <- Reduce("union", lapply(enhancer.grns, "[[", "cells"))
length(enhancer.grn.cells) # 1,808


# Load IgHV bulk data
ighv.bulk <- qs::qread(ighv.file)
class(ighv.bulk)
names(ighv.bulk)
ighv.exp <- ighv.bulk$matrix # get the matrix
ighv.phenotype <- ighv.bulk$phenotype # get the phenotype
ighv.tag <- ighv.bulk$tag # get the list of phenotype
class(ighv.exp)
dim(ighv.exp)
class(ighv.phenotype)
length(ighv.phenotype)
head(ighv.phenotype)
table(ighv.phenotype)
class(ighv.tag)
table(ighv.tag)


# Run Scissor on IgHV
ighv.Scissor.0.5 <- Scissor(bulk_dataset = ighv.exp, sc_dataset = sc.obj, 
                            phenotype = ighv.phenotype, 
                            tag = ighv.tag, alpha = 0.5, 
                  family = "binomial", Save_file = paste0(R.dir, "IgHV.RData")) # run Scissor
qs::qsave(ighv.Scissor.0.5, paste0(R.dir, "Scissor_out_IgHV_0.5.qsave"))
ighv.plus <- ighv.Scissor.0.5$Scissor_pos
ighv.minus <- ighv.Scissor.0.5$Scissor_neg
length(ighv.plus)
length(ighv.minus)


# Select the Scissor+ cells in IgHV
ighv.plus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                      cells = ighv.plus, total = ncol(sc.obj))
ighv.plus.pval
apply(ighv.plus.pval, 2, range)
ighv.plus.ids <- which(ighv.plus.pval$Pvalue < 0.05)
length(ighv.plus.ids)


# Select the Scissor- cells in IgHV
ighv.minus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                           cells = ighv.minus, total = ncol(sc.obj))
ighv.minus.pval
apply(ighv.minus.pval, 2, range)
ighv.minus.ids <- which(ighv.minus.pval$Pvalue < 0.05)
length(ighv.minus.ids)


# Try another parameter
ighv.Scissor.0.1 <- Scissor(bulk_dataset = ighv.exp, sc_dataset = sc.obj, 
                            phenotype = ighv.phenotype, 
                            tag = ighv.tag, alpha = 0.1, 
                            family = "binomial", Save_file = paste0(R.dir, "IgHV.RData")) # run Scissor
qs::qsave(ighv.Scissor.0.1, paste0(R.dir, "Scissor_out_IgHV_0.1.qsave"))
ighv.plus.0.1 <- ighv.Scissor.0.1$Scissor_pos
ighv.minus.0.1 <- ighv.Scissor.0.1$Scissor_neg
length(ighv.plus.0.1)
length(ighv.minus.0.1)


# Select the Scissor+ cells in IgHV (alpha = 0.1)
ighv.plus.pval.0.1 <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                           cells = ighv.plus.0.1, total = ncol(sc.obj))
ighv.plus.pval.0.1
apply(ighv.plus.pval.0.1, 2, range)
ighv.plus.ids.0.1 <- which(ighv.plus.pval.0.1$Pvalue < 0.05)
length(ighv.plus.ids.0.1)


# Select the Scissor- cells in IgHV (alpha = 0.1)
ighv.minus.pval.0.1 <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                            cells = ighv.minus.0.1, total = ncol(sc.obj))
ighv.minus.pval.0.1
apply(ighv.minus.pval.0.1, 2, range)
ighv.minus.ids.0.1 <- which(ighv.minus.pval.0.1$Pvalue < 0.05)
length(ighv.minus.ids.0.1)
