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
setwd(R.dir)
# dir.create(R.dir) 
ighv.file <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/1_classification_GSE113386/bulk.qsave"
TP53.file <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/3_TP53_GSE11038/bulk.qsave"
cox.file <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/2_survival_GSE22762/bulk.qsave"


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
# sc_dataset <- readRDS("/fs/ess/PCON0022/liyang/STREAM/benchmarking/10x_lymph_node_lymphoma_14k.RDS") %>% 
#   GetAssayData(slot = "data", assay = "RNA") # get the expression matrix
# dim(sc_dataset) # 36601 14104
# class(sc_dataset)
# sc.obj <- Seurat_preprocessing(sc_dataset, verbose = F)
# class(sc.obj)
# qs::qsave(sc.obj, paste0(R.dir, "Scissor_processed_sc_obj.qsave"))
# sc.obj <- qs::qread(paste0(R.dir, "Scissor_processed_sc_obj.qsave"))


# Load the eGRNs
enhancer.grns <- qs::qread("/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Rfile/Extended_eGRNs.qsave")
names(enhancer.grns[[1]])
lapply(enhancer.grns, "[[", "genes") %>% sapply(., length)
lapply(enhancer.grns, "[[", "peaks") %>% sapply(., length)
lapply(enhancer.grns, "[[", "cells") %>% sapply(., length)
enhancer.grn.cells <- Reduce("union", lapply(enhancer.grns, "[[", "cells"))
length(enhancer.grn.cells) # 1,808





# Build Seurat objects for T and other cells
lymph.obj <- qs::qread("/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Rfile/Obj_celltyped.qsave")
dim(lymph.obj)
table(lymph.obj$cell.type)
t.celltypes <- levels(lymph.obj$cell.type)[c(2:5, 7:8)]
other.celltypes <- setdiff(levels(lymph.obj$cell.type), t.celltypes)
t.celltypes
other.celltypes
t.cells <- colnames(lymph.obj)[which(lymph.obj$cell.type %in% t.celltypes)]
dim(lymph.obj) # 70469 14104
length(t.cells) # 8900
other.cells <- setdiff(colnames(lymph.obj), t.cells)
length(other.cells) # 5204
sc.m <- GetAssayData(lymph.obj, slot = "data", assay = "RNA")
sc.m[1:5, 1:5]
t.m <- sc.m[, t.cells]
dim(t.m)
other.m <- sc.m[, other.cells]
dim(other.m)




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
ighv.tag


# Run Scissor on IgHV T cells
t.ighv.Scissor.0.5 <- Scissor(bulk_dataset = ighv.exp, sc_dataset = t.dataset, 
                            phenotype = ighv.phenotype, 
                            tag = ighv.tag, alpha = 0.5, 
                  family = "binomial", Save_file = paste0(R.dir, "IgHV_Scissor_T_cells_0.5.RData")) # run Scissor
t.ighv.plus <- t.ighv.Scissor.0.5$Scissor_pos
t.ighv.minus <- t.ighv.Scissor.0.5$Scissor_neg
length(t.ighv.plus)
length(t.ighv.minus)


# Run Scissor on IgHV other cells
other.ighv.Scissor.0.5 <- Scissor(bulk_dataset = ighv.exp, sc_dataset = other.dataset, 
                              phenotype = ighv.phenotype, 
                              tag = ighv.tag, alpha = 0.5, 
                              family = "binomial", Save_file = paste0(R.dir, "IgHV_Scissor_other_cells_0.5.RData")) # run Scissor
other.ighv.plus <- other.ighv.Scissor.0.5$Scissor_pos
other.ighv.minus <- other.ighv.Scissor.0.5$Scissor_neg
length(other.ighv.plus)
length(other.ighv.minus)


# Select the Scissor+ cells in IgHV T cells
t.ighv.plus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                             cells = t.ighv.plus, total = ncol(lymph.obj))
apply(t.ighv.plus.pval, 2, range)
t.ighv.plus.ids <- which(t.ighv.plus.pval$Pvalue < 0.05)
t.ighv.plus.ids # None


# Select the Scissor- cells in IgHV T cells
t.ighv.minus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                              cells = t.ighv.minus, total = ncol(lymph.obj))
t.ighv.minus.pval
apply(t.ighv.minus.pval, 2, range)
t.ighv.minus.ids <- which(t.ighv.minus.pval$Pvalue < 0.05)
t.ighv.minus.ids # None


# Select the Scissor+ cells in ighv other cells
other.ighv.plus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                                 cells = other.ighv.plus, total = ncol(lymph.obj))
apply(other.ighv.plus.pval, 2, range)
other.ighv.plus.ids <- which(other.ighv.plus.pval$Pvalue < 0.05)
other.ighv.plus.ids # 32


# Select the Scissor- cells in ighv other cells
other.ighv.minus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                                  cells = other.ighv.minus, total = ncol(lymph.obj))
other.ighv.minus.pval
apply(other.ighv.minus.pval, 2, range)
other.ighv.minus.ids <- which(other.ighv.minus.pval$Pvalue < 0.05)
other.ighv.minus.ids # None


# # Try another parameter
# ighv.Scissor.0.1 <- Scissor(bulk_dataset = ighv.exp, sc_dataset = sc.obj, 
#                             phenotype = ighv.phenotype, 
#                             tag = ighv.tag, alpha = 0.1, 
#                             family = "binomial", Save_file = paste0(R.dir, "IgHV.RData")) # run Scissor
# qs::qsave(ighv.Scissor.0.1, paste0(R.dir, "Scissor_out_IgHV_0.1.qsave"))
# ighv.plus.0.1 <- ighv.Scissor.0.1$Scissor_pos
# ighv.minus.0.1 <- ighv.Scissor.0.1$Scissor_neg
# length(ighv.plus.0.1)
# length(ighv.minus.0.1)
# 
# 
# # Select the Scissor+ cells in IgHV (alpha = 0.1)
# ighv.plus.pval.0.1 <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
#                                            cells = ighv.plus.0.1, total = ncol(sc.obj))
# ighv.plus.pval.0.1
# apply(ighv.plus.pval.0.1, 2, range)
# ighv.plus.ids.0.1 <- which(ighv.plus.pval.0.1$Pvalue < 0.05)
# length(ighv.plus.ids.0.1)
# 
# 
# # Select the Scissor- cells in IgHV (alpha = 0.1)
# ighv.minus.pval.0.1 <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
#                                             cells = ighv.minus.0.1, total = ncol(sc.obj))
# ighv.minus.pval.0.1
# apply(ighv.minus.pval.0.1, 2, range)
# ighv.minus.ids.0.1 <- which(ighv.minus.pval.0.1$Pvalue < 0.05)
# length(ighv.minus.ids.0.1)







# Load TP53 bulk data
TP53.bulk <- qs::qread(TP53.file)
class(TP53.bulk)
names(TP53.bulk)
TP53.exp <- TP53.bulk$matrix # get the matrix
TP53.phenotype <- TP53.bulk$phenotype # get the phenotype
TP53.tag <- TP53.bulk$tag # get the list of phenotype
class(TP53.exp)
dim(TP53.exp)
class(TP53.phenotype)
length(TP53.phenotype)
head(TP53.phenotype)
table(TP53.phenotype)
class(TP53.tag)
TP53.tag


# Run Scissor on TP53 T cells
t.dataset <- Seurat_preprocessing(t.m, verbose = F)
t.TP53.Scissor.0.5 <- Scissor(bulk_dataset = TP53.exp, sc_dataset = t.dataset, 
                            phenotype = TP53.phenotype, 
                            tag = TP53.tag, alpha = 0.5, 
                            family = "binomial", Save_file = paste0(R.dir, "TP53_T_cells_0.5.RData"))
# qs::qsave(TP53.Scissor.0.5, paste0(R.dir, "Scissor_out_TP53_0.5.qsave"))
t.TP53.plus <- t.TP53.Scissor.0.5$Scissor_pos
t.TP53.minus <- t.TP53.Scissor.0.5$Scissor_neg
length(t.TP53.plus)
length(t.TP53.minus)


# Run Scissor on TP53 other cells
other.dataset <- Seurat_preprocessing(other.m, verbose = F)
other.TP53.Scissor.0.5 <- Scissor(bulk_dataset = TP53.exp, sc_dataset = other.dataset, 
                              phenotype = TP53.phenotype, 
                              tag = TP53.tag, alpha = 0.5, 
                              family = "binomial", Save_file = paste0(R.dir, "TP53_other_cells_0.5.RData"))
other.TP53.plus <- other.TP53.Scissor.0.5$Scissor_pos
other.TP53.minus <- other.TP53.Scissor.0.5$Scissor_neg
length(other.TP53.plus)
length(other.TP53.minus)


# Select the Scissor+ cells in TP53 T cells
t.TP53.plus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                           cells = t.TP53.plus, total = ncol(lymph.obj))
apply(t.TP53.plus.pval, 2, range)
t.TP53.plus.ids <- which(t.TP53.plus.pval$Pvalue < 0.05)
t.TP53.plus.ids # none


# Select the Scissor- cells in TP53 T cells
t.TP53.minus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                            cells = t.TP53.minus, total = ncol(lymph.obj))
t.TP53.minus.pval
apply(t.TP53.minus.pval, 2, range)
t.TP53.minus.ids <- which(t.TP53.minus.pval$Pvalue < 0.05)
t.TP53.minus.ids # 14


# Select the Scissor+ cells in TP53 other cells
other.TP53.plus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                             cells = other.TP53.plus, total = ncol(lymph.obj))
apply(other.TP53.plus.pval, 2, range)
other.TP53.plus.ids <- which(other.TP53.plus.pval$Pvalue < 0.05)
other.TP53.plus.ids


# Select the Scissor- cells in TP53 other cells
other.TP53.minus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                              cells = other.TP53.minus, total = ncol(lymph.obj))
other.TP53.minus.pval
apply(other.TP53.minus.pval, 2, range)
other.TP53.minus.ids <- which(other.TP53.minus.pval$Pvalue < 0.05)
other.TP53.minus.ids


# Too many other plus & minus cells
intersect(other.TP53.plus.ids, other.TP53.minus.ids) %>% length # 20
other.tp53.common.ids <- intersect(other.TP53.plus.ids, other.TP53.minus.ids)
sapply(enhancer.grns[other.tp53.common.ids], "[[", "TF")
length(other.TP53.plus.ids) # 32
length(other.TP53.minus.ids) # 23
sapply(enhancer.grns[setdiff(other.TP53.plus.ids, 
                             other.TP53.minus.ids)], "[[", "TF")
sapply(enhancer.grns[setdiff(other.TP53.minus.ids, 
                             other.TP53.plus.ids)], "[[", "TF")






# Load cox bulk data
cox.bulk <- qs::qread(cox.file)
class(cox.bulk)
names(cox.bulk)
cox.exp <- cox.bulk$matrix # get the matrix
cox.surv <- cox.bulk$phenotype # get the phenotype
all(colnames(cox.exp) == cox.surv$Id)
cox.phenotype <- cox.surv[, 2:3]
colnames(cox.phenotype) <- c("time", "status")
head(cox.phenotype)
class(cox.exp)
dim(cox.exp)




# Run Scissor on cox T cells
t.cox.Scissor <- Scissor(bulk_dataset = cox.exp, sc_dataset = t.dataset, 
                              phenotype = cox.phenotype, 
                              alpha = 0.02, 
                              family = "cox", Save_file = paste0(R.dir, "cox_T_cells_0.5.RData"))
t.cox.plus <- t.cox.Scissor$Scissor_pos
t.cox.minus <- t.cox.Scissor$Scissor_neg
length(t.cox.plus)
length(t.cox.minus)


# Run Scissor on cox other cells
other.cox.Scissor <- Scissor(bulk_dataset = cox.exp, sc_dataset = other.dataset, 
                                  phenotype = cox.phenotype, 
                                  alpha = 0.01, 
                                  family = "cox", 
                                 Save_file = paste0(R.dir, "cox_other_cells_0.5.RData"))
other.cox.plus <- other.cox.Scissor$Scissor_pos
other.cox.minus <- other.cox.Scissor$Scissor_neg
length(other.cox.plus)
length(other.cox.minus)


# Select the Scissor+ cells in cox T cells
t.cox.plus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                             cells = t.cox.plus, total = ncol(lymph.obj))
apply(t.cox.plus.pval, 2, range)
t.cox.plus.pval
t.cox.plus.ids <- which(t.cox.plus.pval$Pvalue < 0.05)
t.cox.plus.ids # none


# Select the Scissor- cells in cox T cells
t.cox.minus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                              cells = t.cox.minus, total = ncol(lymph.obj))
t.cox.minus.pval
apply(t.cox.minus.pval, 2, range)
t.cox.minus.ids <- which(t.cox.minus.pval$Pvalue < 0.05)
t.cox.minus.ids # None


# Select the Scissor+ cells in cox other cells
other.cox.plus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                                 cells = other.cox.plus, total = ncol(lymph.obj))
apply(other.cox.plus.pval, 2, range)
other.cox.plus.ids <- which(other.cox.plus.pval$Pvalue < 0.05)
other.cox.plus.ids
length(other.cox.plus.ids) # 17


# Select the Scissor- cells in cox other cells
other.cox.minus.pval <- hyper_test_enhancer_grns(enhancer.grns = enhancer.grns, 
                                                  cells = other.cox.minus, total = ncol(lymph.obj))
other.cox.minus.pval
apply(other.cox.minus.pval, 2, range)
other.cox.minus.ids <- which(other.cox.minus.pval$Pvalue < 0.05)
other.cox.minus.ids
length(other.cox.minus.ids) # 44
