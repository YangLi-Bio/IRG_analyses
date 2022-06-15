##################################################################################
#                                                                                #
#                            Immune Related Genes (IRGs)                         #
#                                                                                #
##################################################################################


##################################################################################
#                                                                                #
#                                  Basic parameters                              #
#                                                                                #
##################################################################################


# Libraries
library(Seurat)
library(dplyr)
library(qs)
library(pbapply)


# Parameters
source("/fs/ess/PCON0022/liyang/r_utilities/functions/visual_tools.R")
work.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/"
egrn.dir <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/STREAM/10x_lymph_node_lymphoma_14k/"
phenotype.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/1_classification_GSE113386/"
param.file <- "hbc_15_var_3000_top_3000_cons_1.0.qsave"
boxplot.path <- paste0(work.dir, "genes_CREs_cells_IRG_ratios.eps")
rds.path <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/10x_lymph_node_lymphoma_14k.RDS"
setwd(work.dir)
getwd()


# Load R objects
irg.ratios <- qs::qread(paste0(work.dir, "irg_ratios.qsave")) # ranked ratios of overlapped IRGs
scissor.res <- qs::qread(paste0(phenotype.dir, "GSE113386_single_param.qsave"))
length(irg.ratios)
range(irg.ratios)


##################################################################################
#                                                                                #
#    1. Boxplots of the number of cells, genes, CREs, and IRG ratios in eGRNs    #
#                                                                                #
##################################################################################


egrns.ll <- qs::qread(paste0(egrn.dir, param.file))
library(dplyr)
genes.ll <- lapply(egrns.ll, "[[", "genes") %>% sapply(., length)
range(genes.ll)
p.genes <- get_boxplot(obj = genes.ll, path = NULL,
            y.lab = "Genes", color = "black", fill = "#BE2A3E")
print(p.genes)
peaks.ll <- lapply(egrns.ll, "[[", "peaks") %>% sapply(., length)
range(peaks.ll)
p.peaks <- get_boxplot(obj = peaks.ll, path = NULL,
                       y.lab = "CREs", color = "black", fill = "#BE2A3E")
print(p.peaks)
cells.ll <- lapply(egrns.ll, "[[", "cells") %>% sapply(., length)
range(peaks.ll)
p.cells <- get_boxplot(obj = cells.ll, path = NULL,
                       y.lab = "Cells", color = "black", fill = "#BE2A3E")
print(p.cells)
p.irg <- get_boxplot(obj = irg.ratios * 100, path = NULL,
            y.lab = "IRG ratio (%)", color = "black", fill = "#BE2A3E")
print(p.irg)
library(ggplotify)
p.boxplot <- as.ggplot(p.genes | p.peaks | p.cells | p.irg)
qs::qsave(p.boxplot, paste0(work.dir, "Genes_CREs_cells_IRGs.qsave"))


##################################################################################
#                                                                                #
#                 2. Ten eGRNs distinguished IgHV v.s. others                   #
#                                                                                #
##################################################################################


# Load the IgHV phenotype data
tumor.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/1_classification_GSE113386/"
tumor.scissor.file <- "GSE113386_single_param.qsave"
tumor.scissor.cells <- qs::qread(paste0(tumor.dir, tumor.scissor.file))
names(tumor.scissor.cells)
lymph.obj <- qs::qread("/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/lymph_obj.qsave")
ncol(lymph.obj)
tumor.plus.cells <- tumor.scissor.cells$Scissor_pos
tumor.minus.cells <- tumor.scissor.cells$Scissor_neg
length(tumor.plus.cells) / ncol(lymph.obj)
length(tumor.minus.cells) / ncol(lymph.obj)
library(data.table)
tumor.dt <- rbindlist(lapply(irg.egrns, function(x) {
  plus.ratio <- length(intersect(x$cells, tumor.plus.cells)) / length(x$cells)
  minus.ratio <- length(intersect(x$cells, tumor.minus.cells)) / length(x$cells)
  list(TF = x$TF, Plus = plus.ratio, Minus = minus.ratio)
}))
tumor.dt <- cbind(eGRN = seq_along(irg.egrns), tumor.dt)
head(tumor.dt)
tumor.dt
range(tumor.dt$Minus)
range(tumor.dt$Plus)
tumor.dt[tumor.dt$TF == "KLF5",]
max(tumor.dt$Plus)
which.max(tumor.dt$Plus)
max(tumor.dt$Minus)
which.max(tumor.dt$Minus)
tumor.phenotype <- c(rep("IgHV", length(tumor.plus.cells)), 
                    rep("unmutated", length(tumor.minus.cells)), 
                    rep("background", ncol(lymph.obj) - length(tumor.plus.cells) - 
                          length(tumor.minus.cells)))
length(tumor.phenotype)
names(tumor.phenotype) <- c(tumor.plus.cells, tumor.minus.cells, setdiff(colnames(lymph.obj), 
                                                                      c(tumor.plus.cells, 
                                                                        tumor.minus.cells)))
head(tumor.phenotype)
lymph.obj <- AddMetaData(lymph.obj, metadata = tumor.phenotype, 
                         col.name = "IgHV.phenotype")
lymph.obj$IgHV.phenotype <- factor(lymph.obj$IgHV.phenotype, 
                                   levels = c("IgHV", "background", 
                                              "unmutated"))
qs::qsave(lymph.obj, paste0(work.dir, "lymph_obj.qsave"))


# DEG analysis
length(tumor.plus.cells) / length(tumor.minus.cells)
library(Seurat)
Idents(lymph.obj) <- lymph.obj$IgHV.phenotype
tumor.DEGs <- FindMarkers(lymph.obj, ident.1 = "tumor")
dim(tumor.DEGs)
colnames(tumor.DEGs)
qs::qsave(tumor.DEGs, paste0(tumor.dir, "DEGs.qsave"))
# tumor.DEGs <- qs::qread(paste0(tumor.dir, "DEGs.qsave"))


# Overrepresented and underrepresented DEGs and eGRNs
tumor.up.DEGs <- rownames(tumor.DEGs[tumor.DEGs$p_val_adj < 0.05 & 
                                     tumor.DEGs$avg_log2FC > 0.25,])
length(tumor.up.DEGs)
tumor.down.DEGs <- rownames(tumor.DEGs[tumor.DEGs$p_val_adj < 0.05 & 
                                       tumor.DEGs$avg_log2FC < -0.25,])
length(tumor.down.DEGs)
tumor.eGRNs <- lapply(irg.egrns, function(x) {
  xx <- x
  xx$up <- intersect(tumor.up.DEGs, xx$genes)
  xx$down <- intersect(tumor.down.DEGs, xx$genes)
  xx
})
qs::qsave(tumor.eGRNs, paste0(tumor.dir, "DEG_eGRNs.qsave"))
# tumor.eGRNs <- qs::qread(paste0(tumor.dir, "DEG_eGRNs.qsave"))
tumor.up.DEG.ratio <- sapply(tumor.eGRNs, function(x) {
  length(x$up) / length(x$genes)
})
tumor.down.DEG.ratio <- sapply(tumor.eGRNs, function(x) {
  length(x$down) / length(x$genes)
})
which.max(tumor.up.DEG.ratio)
which.max(tumor.down.DEG.ratio)


# Obtain the top-five TFs regulating the overrepresented DEGs in eGRNs
head(sapply(tumor.eGRNs[order(tumor.up.DEG.ratio, decreasing = T)], "[[", "TF")) # overrepresented DEGs
# NFIC : important in lymphoma (https://www.proteinatlas.org/ENSG00000141905-NFIC/pathology/lymphoma)
# MYCN (x 2) : important in CLL (chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.pnas.org/doi/pdf/10.1073/pnas.0304717101)
# BHLHE22 : important in lymphoma (https://www.proteinatlas.org/ENSG00000180828-BHLHE22/pathology/lymphoma)
# ATF2 : important in lymphoma (https://www.proteinatlas.org/ENSG00000115966-ATF2/pathology/lymphoma)
# NR4A1 : important in lymphoma (https://www.nature.com/articles/s41598-018-32972-4)
# GABPA (x 3) : important in CLL (https://maayanlab.cloud/Harmonizome/gene_set/GABPA/ENCODE+Transcription+Factor+Targets)
# KLF5 (x 2) : important in lymphoma (https://www.proteinatlas.org/ENSG00000102554-KLF5/pathology/lymphoma)
# MYOG (x 3) : important in CLL (https://maayanlab.cloud/Harmonizome/gene_set/MYOG/ENCODE+Transcription+Factor+Targets)
# SP1 (x 3) : CLL (https://www.embopress.org/doi/full/10.15252/msb.20188339)
# STAT3 (x 3) : CLL (https://pubmed.ncbi.nlm.nih.gov/32075387/)
# TCF12 (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TCF12/ENCODE+Transcription+Factor+Targets)
# TCF3 (x 3) : CLL (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8096467/)
# TFAP2C (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TFAP2C/CHEA+Transcription+Factor+Targets)
# ZNF682 (x 2) : lymphoma (https://www.proteinatlas.org/ENSG00000197124-ZNF682/pathology/lymphoma)


head(sapply(tumor.eGRNs[order(tumor.down.DEG.ratio, decreasing = T)], "[[", "TF")) # underrepresented DEGs
# EGR1 (x 2) : lymphoma (https://www.karger.com/Article/Fulltext/494867)
# GABPA (x 3) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/GABPA/ENCODE+Transcription+Factor+Targets)
# KLF5 (x 2) : lymphoma (https://www.proteinatlas.org/ENSG00000102554-KLF5/pathology/lymphoma)
# MYCN (x 2) : CLL (chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.pnas.org/doi/pdf/10.1073/pnas.0304717101)
# MYOG (x 3) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/MYOG/ENCODE+Transcription+Factor+Targets)
# SP1 (x 3) : CLL (https://www.embopress.org/doi/full/10.15252/msb.20188339)
# STAT3 (x 3) : CLL (https://pubmed.ncbi.nlm.nih.gov/32075387/)
# TCF12 (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TCF12/ENCODE+Transcription+Factor+Targets)
# TCF3 (x 3) : CLL (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8096467/)
# TFAP2C (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TFAP2C/CHEA+Transcription+Factor+Targets)
# ZNF682 (x 2) : lymphoma (https://www.proteinatlas.org/ENSG00000197124-ZNF682/pathology/lymphoma)


# Co-factors
# MYCN : SP1
# SP1 : EGR1, GABPA, JUN, ESRRA, KLF4, MYCN, MYOG, STAT3, TP53, YY1
# EGR1 : SP1, TP53
# GABPA : SP1, YY1
# JUN : ATF2, FOS, KLF5, SP1, STAT3, TP53
# MYOG : SP1
# STAT3 : BATF3, JUN, SP1
# BATF3 : STAT3
# TFAP2C : TP53


# Build bulk expression table for GSVA analysis
source("/fs/ess/PCON0022/liyang/r_utilities/functions/transcriptome_tools.R")
tumor.gene.set <- lapply(tumor.eGRNs, function(x) {
  #x$genes
  x$up
  #x$down
})
names(tumor.gene.set) <- seq_along(tumor.gene.set)
tumor.bulk <- qs::qread("/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/1_classification_GSE113386/bulk.qsave")
tumor.bulk.expr <- tumor.bulk$matrix
tumor.bulk.meta <- tumor.bulk$phenotype
colnames(tumor.bulk.expr)
dim(tumor.bulk.expr)
tumor.gsva.m <- run_GSVA(X = tumor.bulk.expr, gs = tumor.gene.set)
dim(tumor.gsva.m)
tumor.gsva.m[1:5, 1:5]
rownames(tumor.gsva.m) <- names(tumor.gene.set)


# Generate the UMAP between IgHV v.s. unmutated
lymph.obj <- Seurat_reduce_dim(object = lymph.obj)
qs::qsave(lymph.obj, paste0(work.dir, "lymph_obj.qsave"))
Idents(lymph.obj) <- lymph.obj$IgHV.phenotype
p.IgHV.UMAP <- get_UMAP(object = lymph.obj, reduction.method = "wnn.umap",
                   txt = "IgHV")
qs::qsave(p.IgHV.UMAP, paste0(phenotype.dir, "UMAP.qsave"))


# Generate boxplots for the two MYCN eGRNs
tumor.MYCN.ids <- which(sapply(tumor.eGRNs, "[[", "TF") == "MYCN")
MYCN.tumor.pval <- lapply(tumor.MYCN.ids, function(i) {
  x <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 1))]
  y <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 0))]
  wilcox.test(x, y, alternative = "greater")$p.value
})
names(MYCN.tumor.pval) <- tumor.MYCN.ids


# Generate boxplot for the first eGRN of MYCN
i <- tumor.MYCN.ids[1]
x <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 1))]
y <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 0))]
obj <- data.frame(Sample = c(rep("IgHV", length(x)), rep("unmutated", length(y))), 
                  GSVA = c(x, y))
obj$Sample <- factor(obj$Sample, levels = c("IgHV", "unmutated"))
p.tumor.MYCN1.box <- get_boxplot(obj = obj, path = NULL, x.lab = "Sample",
                       y.lab = "GSVA", color = "black", title = "MYCN", 
                       fill = c("IgHV" = "red", "unmutated" = "blue"))


# Generate boxplot for the second eGRN of MYCN
i <- tumor.MYCN.ids[2]
x <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 1))]
y <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 0))]
obj <- data.frame(Sample = c(rep("IgHV", length(x)), rep("unmutated", length(y))), 
                  GSVA = c(x, y))
obj$Sample <- factor(obj$Sample, levels = c("IgHV", "unmutated"))
p.tumor.MYCN2.box <- get_boxplot(obj = obj, path = NULL, x.lab = "Sample",
                                 y.lab = "GSVA", color = "black", title = "MYCN", 
                                 fill = c("IgHV" = "red", "unmutated" = "blue"))


# Generate boxplot for JUN (co-factor of STAT3)
tumor.JUN.ids <- which(sapply(tumor.eGRNs, "[[", "TF") == "JUN")
i <- 1
tumor.JUN.pval <- lapply(tumor.JUN.ids, function(i) {
  x <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 1))]
  y <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 0))]
  wilcox.test(x, y, alternative = "greater")$p.value
})
names(tumor.JUN.pval) <- tumor.JUN.ids


# Generate boxplot for JUN
i <- tumor.JUN.ids[1]
x <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 1))]
y <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 0))]
obj <- data.frame(Sample = c(rep("IgHV", length(x)), rep("unmutated", length(y))), 
                  GSVA = c(x, y))
obj$Sample <- factor(obj$Sample, levels = c("IgHV", "unmutated"))
p.tumor.JUN.box <- get_boxplot(obj = obj, path = NULL, x.lab = "Sample",
                                 y.lab = "GSVA", color = "black", title = "JUN", 
                                 fill = c("IgHV" = "red", "unmutated" = "blue"))


# Generate boxplot for STAT3 (co-factor of JUN)
tumor.STAT3.ids <- which(sapply(tumor.eGRNs, "[[", "TF") == "STAT3")
tumor.STAT3.pval <- lapply(tumor.STAT3.ids, function(i) {
  x <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 1))]
  y <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 0))]
  wilcox.test(x, y, alternative = "greater")$p.value
})
names(tumor.STAT3.pval) <- tumor.STAT3.ids


# Generate boxplot for STAT3
i <- tumor.STAT3.ids[2]
x <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 1))]
y <- tumor.gsva.m[i, names(which(tumor.bulk.meta == 0))]
obj <- data.frame(Sample = c(rep("IgHV", length(x)), rep("unmutated", length(y))), 
                  GSVA = c(x, y))
obj$Sample <- factor(obj$Sample, levels = c("IgHV", "unmutated"))
p.tumor.STAT3.box <- get_boxplot(obj = obj, path = NULL, x.lab = "Sample",
                               y.lab = "GSVA", color = "black", title = "STAT3", 
                               fill = c("IgHV" = "red", "unmutated" = "blue"))


# Merge the boxplots
p.IgHV.box <- p.IgHV.UMAP | p.tumor.MYCN1.box | p.tumor.MYCN2.box | p.tumor.JUN.box | p.tumor.STAT3.box
qs::qsave(p.tumor.box, paste0(work.dir, "Boxplots_GSVA_IgHV.qsave"))


##################################################################################
#                                                                                #
#                 3. Ten eGRNs distinguished TP53 v.s. unmutated                 #
#                                                                                #
##################################################################################


# Load the TP53 phenotype data
TP53.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/3_TP53_GSE11038/"
TP53.scissor.file <- "3_TP53_GSE11038_single_param.qsave"
TP53.scissor.cells <- qs::qread(paste0(TP53.dir, TP53.scissor.file))
names(TP53.scissor.cells)
# lymph.obj <- readRDS(rds.path) # contains all the cells
ncol(lymph.obj)
TP53.plus.cells <- TP53.scissor.cells$Scissor_pos
TP53.minus.cells <- TP53.scissor.cells$Scissor_neg
length(TP53.plus.cells) / ncol(lymph.obj)
length(TP53.minus.cells) / ncol(lymph.obj)
library(data.table)
TP53.dt <- rbindlist(lapply(irg.egrns, function(x) {
  plus.ratio <- length(intersect(x$cells, TP53.plus.cells)) / length(x$cells)
  minus.ratio <- length(intersect(x$cells, TP53.minus.cells)) / length(x$cells)
  list(TF = x$TF, Plus = plus.ratio, Minus = minus.ratio)
}))
TP53.dt <- cbind(eGRN = seq_along(irg.egrns), TP53.dt)
head(TP53.dt)
TP53.dt
TP53.dt[TP53.dt$TF == "KLF5",]
max(TP53.dt$Plus)
which.max(TP53.dt$Plus)
max(TP53.dt$Minus)
which.max(TP53.dt$Minus)
TP53.phenotype <- c(rep("TP53", length(TP53.plus.cells)), 
                    rep("unmutated", length(TP53.minus.cells)), 
                    rep("background", ncol(lymph.obj) - length(TP53.plus.cells) - 
                          length(TP53.minus.cells)))
length(TP53.phenotype)
names(TP53.phenotype) <- c(TP53.plus.cells, TP53.minus.cells, setdiff(colnames(lymph.obj), 
                                                                      c(TP53.plus.cells, 
                                                                        TP53.minus.cells)))
head(TP53.phenotype)
lymph.obj <- AddMetaData(lymph.obj, metadata = TP53.phenotype, 
                         col.name = "TP53.phenotype")
lymph.obj$TP53.phenotype <- factor(lymph.obj$TP53.phenotype, 
                                   levels = c("TP53", "background", 
                                              "unmutated"))
qs::qsave(lymph.obj, paste0(work.dir, "lymph_obj.qsave"))
# lymph.obj <- qs::qread(paste0(work.dir, "lymph_obj.qsave"))


# DEG analysis
length(TP53.plus.cells) / length(TP53.minus.cells)
library(Seurat)
Idents(lymph.obj) <- lymph.obj$TP53.phenotype
TP53.DEGs <- FindMarkers(lymph.obj, ident.1 = "TP53", 
                         ident.2 = "unmutated")
dim(TP53.DEGs)
colnames(TP53.DEGs)
qs::qsave(TP53.DEGs, paste0(TP53.dir, "DEGs.qsave"))
# TP53.DEGs <- qs::qread(paste0(TP53.dir, "DEGs.qsave"))


# Overrepresented and underrepresented DEGs and eGRNs
TP53.up.DEGs <- rownames(TP53.DEGs[TP53.DEGs$p_val_adj < 0.05 & 
                                     TP53.DEGs$avg_log2FC > 0.25,])
length(TP53.up.DEGs)
TP53.down.DEGs <- rownames(TP53.DEGs[TP53.DEGs$p_val_adj < 0.05 & 
                                     TP53.DEGs$avg_log2FC < -0.25,])
length(TP53.down.DEGs)
TP53.eGRNs <- lapply(irg.egrns, function(x) {
  xx <- x
  xx$up <- intersect(TP53.up.DEGs, xx$genes)
  xx$down <- intersect(TP53.down.DEGs, xx$genes)
  xx
})
qs::qsave(TP53.eGRNs, paste0(TP53.dir, "DEG_eGRNs.qsave"))
TP53.up.DEG.ratio <- sapply(TP53.eGRNs, function(x) {
  length(x$up) / length(x$genes)
})
TP53.down.DEG.ratio <- sapply(TP53.eGRNs, function(x) {
  length(x$down) / length(x$genes)
})
which.max(TP53.up.DEG.ratio)
which.max(TP53.down.DEG.ratio)


# Obtain the top-five TFs regulating the overrepresented DEGs in eGRNs
head(sapply(TP53.eGRNs[order(TP53.up.DEG.ratio, decreasing = T)], "[[", "TF")) # overrepresented DEGs
# NFIC : important in lymphoma (https://www.proteinatlas.org/ENSG00000141905-NFIC/pathology/lymphoma)
# MYCN (x 2) : important in CLL (chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.pnas.org/doi/pdf/10.1073/pnas.0304717101)
# BHLHE22 : important in lymphoma (https://www.proteinatlas.org/ENSG00000180828-BHLHE22/pathology/lymphoma)
# ATF2 : important in lymphoma (https://www.proteinatlas.org/ENSG00000115966-ATF2/pathology/lymphoma)
# NR4A1 : important in lymohoma (https://www.nature.com/articles/s41598-018-32972-4)
# GABPA (x 3) : important in CLL (https://maayanlab.cloud/Harmonizome/gene_set/GABPA/ENCODE+Transcription+Factor+Targets)
# KLF5 (x 2) : important in lymphoma (https://www.proteinatlas.org/ENSG00000102554-KLF5/pathology/lymphoma)
# MYOG (x 3) : important in CLL (https://maayanlab.cloud/Harmonizome/gene_set/MYOG/ENCODE+Transcription+Factor+Targets)
# SP1 (x 3) : CLL (https://www.embopress.org/doi/full/10.15252/msb.20188339)
# STAT3 (x 3) : CLL (https://pubmed.ncbi.nlm.nih.gov/32075387/)
# TCF12 (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TCF12/ENCODE+Transcription+Factor+Targets)
# TCF3 (x 3) : CLL (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8096467/)
# TFAP2C (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TFAP2C/CHEA+Transcription+Factor+Targets)
# ZNF682 (x 2) : lymphoma (https://www.proteinatlas.org/ENSG00000197124-ZNF682/pathology/lymphoma)


head(sapply(TP53.eGRNs[order(TP53.down.DEG.ratio, decreasing = T)], "[[", "TF")) # underrepresented DEGs
# EGR1 (x 2) : lymphoma (https://www.karger.com/Article/Fulltext/494867)
# GABPA (x 3) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/GABPA/ENCODE+Transcription+Factor+Targets)
# KLF5 (x 2) : lymphoma (https://www.proteinatlas.org/ENSG00000102554-KLF5/pathology/lymphoma)
# MYCN (x 2) : CLL (chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.pnas.org/doi/pdf/10.1073/pnas.0304717101)
# MYOG (x 3) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/MYOG/ENCODE+Transcription+Factor+Targets)
# SP1 (x 3) : CLL (https://www.embopress.org/doi/full/10.15252/msb.20188339)
# STAT3 (x 3) : CLL (https://pubmed.ncbi.nlm.nih.gov/32075387/)
# TCF12 (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TCF12/ENCODE+Transcription+Factor+Targets)
# TCF3 (x 3) : CLL (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8096467/)
# TFAP2C (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TFAP2C/CHEA+Transcription+Factor+Targets)
# ZNF682 (x 2) : lymphoma (https://www.proteinatlas.org/ENSG00000197124-ZNF682/pathology/lymphoma)


# Co-factors
# MYCN : SP1
# SP1 : EGR1, GABPA, JUN, ESRRA, KLF4, MYCN, MYOG, STAT3, TP53, YY1
# EGR1 : SP1, TP53
# GABPA : SP1, YY1
# JUN : ATF2, FOS, KLF5, SP1, STAT3, TP53
# MYOG : SP1
# STAT3 : BATF3, JUN, SP1
# BATF3 : STAT3
# TFAP2C : TP53


# Generate the UMAP between TP53 v.s. unmutated
Idents(lymph.obj) <- lymph.obj$TP53.phenotype
p.TP53.UMAP <- get_UMAP(object = lymph.obj, reduction.method = "wnn.umap",
                        txt = "TP53")


# Build bulk expression table for GSVA analysis
TP53.gene.set <- lapply(TP53.eGRNs, function(x) {
  #x$genes
  x$up
  #x$down
})
names(TP53.gene.set) <- seq_along(TP53.gene.set)
TP53.bulk <- qs::qread("/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/3_TP53_GSE11038/bulk.qsave")
TP53.bulk.expr <- TP53.bulk$matrix
TP53.bulk.meta <- TP53.bulk$phenotype
colnames(TP53.bulk.expr)
dim(TP53.bulk.expr)
TP53.gsva.m <- run_GSVA(X = TP53.bulk.expr, gs = TP53.gene.set)
dim(TP53.gsva.m)
colnames(TP53.gsva.m) <- gsub(".CEL.gz", "", colnames(TP53.gsva.m))
# TP53.col.names <- gsub(".CEL.gz", "", colnames(TP53.gsva.m))
TP53.bulk.meta <- TP53.bulk.meta[gsub(".CEL.gz", "", colnames(TP53.gsva.m))]
TP53.gsva.m[1:5, 1:5]
rownames(TP53.gsva.m) <- names(TP53.gene.set)


# Generate boxplots for the three MYOG eGRNs
TP53.MYOG.ids <- which(sapply(TP53.eGRNs, "[[", "TF") == "MYOG")
length(TP53.MYOG.ids)


# TP53.pval <- unlist(lapply(seq_along(TP53.eGRNs), function(i) {
#   x <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 1))]
#   y <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 0))]
#   wilcox.test(x, y, alternative = "greater")$p.value
# }))
# names(TP53.pval) <- seq_along(TP53.eGRNs)
# table(sapply(TP53.eGRNs, "[[", "TF")[as.numeric(names(TP53.pval[which(TP53.pval < 0.05)]))])


MYOG.TP53.pval <- lapply(TP53.MYOG.ids, function(i) {
  x <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 1))]
  y <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 0))]
  wilcox.test(x, y, alternative = "greater")$p.value
})
names(MYOG.TP53.pval) <- TP53.MYOG.ids


# Generate boxplot for the first eGRN of MYOG
i <- TP53.MYOG.ids[1]
x <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 1))]
y <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 0))]
obj <- data.frame(Sample = c(rep("TP53", length(x)), rep("unmutated", length(y))), 
                  GSVA = c(x, y))
obj$Sample <- factor(obj$Sample, levels = c("TP53", "unmutated"))
p.TP53.MYOG1.box <- get_boxplot(obj = obj, path = NULL, x.lab = "Sample",
                                 y.lab = "GSVA", color = "black", title = "MYOG", 
                                 fill = c("TP53" = "red", "unmutated" = "blue"))


# Generate boxplot for the second eGRN of MYOG
i <- TP53.MYOG.ids[2]
x <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 1))]
y <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 0))]
obj <- data.frame(Sample = c(rep("TP53", length(x)), rep("unmutated", length(y))), 
                  GSVA = c(x, y))
obj$Sample <- factor(obj$Sample, levels = c("TP53", "unmutated"))
p.TP53.MYOG2.box <- get_boxplot(obj = obj, path = NULL, x.lab = "Sample",
                                 y.lab = "GSVA", color = "black", title = "MYOG", 
                                 fill = c("TP53" = "red", "unmutated" = "blue"))


# Generate boxplot for JUN (co-factor of STAT3)
TP53.JUN.ids <- which(sapply(TP53.eGRNs, "[[", "TF") == "JUN")
i <- 1
TP53.JUN.pval <- lapply(TP53.JUN.ids, function(i) {
  x <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 1))]
  y <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 0))]
  wilcox.test(x, y, alternative = "greater")$p.value
})
names(TP53.JUN.pval) <- TP53.JUN.ids


# Generate boxplot for JUN
i <- TP53.JUN.ids[1]
x <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 1))]
y <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 0))]
obj <- data.frame(Sample = c(rep("TP53", length(x)), rep("unmutated", length(y))), 
                  GSVA = c(x, y))
obj$Sample <- factor(obj$Sample, levels = c("TP53", "unmutated"))
p.TP53.JUN.box <- get_boxplot(obj = obj, path = NULL, x.lab = "Sample",
                               y.lab = "GSVA", color = "black", title = "JUN", 
                               fill = c("TP53" = "red", "unmutated" = "blue"))


# Generate boxplot for STAT3 (co-factor of JUN)
TP53.STAT3.ids <- which(sapply(TP53.eGRNs, "[[", "TF") == "STAT3")
TP53.STAT3.pval <- lapply(TP53.STAT3.ids, function(i) {
  x <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 1))]
  y <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 0))]
  wilcox.test(x, y, alternative = "greater")$p.value
})
names(TP53.STAT3.pval) <- TP53.STAT3.ids


# Generate boxplot for STAT3
i <- TP53.STAT3.ids[2]
x <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 1))]
y <- TP53.gsva.m[i, names(which(TP53.bulk.meta == 0))]
obj <- data.frame(Sample = c(rep("TP53", length(x)), rep("unmutated", length(y))), 
                  GSVA = c(x, y))
obj$Sample <- factor(obj$Sample, levels = c("TP53", "unmutated"))
p.TP53.STAT3.box <- get_boxplot(obj = obj, path = NULL, x.lab = "Sample",
                                 y.lab = "GSVA", color = "black", title = "STAT3", 
                                 fill = c("TP53" = "red", "unmutated" = "blue"))


# Merge the boxplots
p.TP53.box <- p.TP53.UMAP | p.TP53.MYOG1.box | p.TP53.MYOG2.box | p.TP53.JUN.box | p.TP53.STAT3.box
qs::qsave(phenotype.dir, "Boxplots_GSVA_TP53.qsave")


##################################################################################
#                                                                                #
#                   12. Ten eGRNs predict survival outcomes                      #
#                                                                                #
##################################################################################


# Obtain the Scissor+ and Scissor- cells
cox.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/2_survival_GSE22762/"
cox.scissor.file <- "2_survival_GSE22762_optimal_alpha_0.017.qsave"
cox.scissor.cells <- qs::qread(paste0(cox.dir, cox.scissor.file))
names(cox.scissor.cells)
cox.plus.cells <- cox.scissor.cells$Scissor_pos
length(cox.plus.cells)
cox.minus.cells <- cox.scissor.cells$Scissor_neg
length(cox.minus.cells)
ncol(lymph.obj)
length(cox.plus.cells) / ncol(lymph.obj)
length(cox.minus.cells) / ncol(lymph.obj)
cox.phenotype <- c(rep("good", length(cox.plus.cells)), 
                   rep("bad", length(cox.minus.cells)), 
                   rep("background", ncol(lymph.obj) - length(cox.plus.cells) - 
                         length(cox.minus.cells)))
head(cox.phenotype)
names(cox.phenotype) <- c(cox.plus.cells, cox.minus.cells, setdiff(colnames(lymph.obj), 
                                                                   c(cox.plus.cells, 
                                                                     cox.minus.cells)))
head(cox.phenotype)
cox.phenotype <- factor(cox.phenotype, levels = c("good", "background", "bad"))
cox.phenotype
lymph.obj <- AddMetaData(object = lymph.obj, metadata = cox.phenotype, col.name = "survival")
colnames(lymph.obj@meta.data)
qs::qsave(lymph.obj, paste0(work.dir, "lymph_obj.qsave"))


# DEG analysis
length(cox.plus.cells) / length(cox.minus.cells)
Idents(lymph.obj) <- lymph.obj$survival
Idents(lymph.obj)
library(Seurat)
cox.DEGs <- FindMarkers(lymph.obj, ident.1 = "good", ident.2 = "bad")
dim(cox.DEGs)
colnames(cox.DEGs)
qs::qsave(cox.DEGs, paste0(cox.dir, "DEGs.qsave"))
# cox.DEGs <- qs::qread(paste0(cox.dir, "DEGs.qsave"))


# Overrepresented and underrepresented DEGs and eGRNs
cox.up.DEGs <- rownames(cox.DEGs[cox.DEGs$p_val_adj < 0.05 & 
                                     cox.DEGs$avg_log2FC > 0.25,])
length(cox.up.DEGs)
cox.down.DEGs <- rownames(cox.DEGs[cox.DEGs$p_val_adj < 0.05 & 
                                       cox.DEGs$avg_log2FC < -0.25,])
length(cox.down.DEGs)
cox.eGRNs <- lapply(irg.egrns, function(x) {
  xx <- x
  xx$up <- intersect(cox.up.DEGs, xx$genes)
  xx$down <- intersect(cox.down.DEGs, xx$genes)
  xx
})
qs::qsave(cox.eGRNs, paste0(cox.dir, "DEG_eGRNs.qsave"))
cox.eGRNs <- qs::qread(paste0(cox.dir, "DEG_eGRNs.qsave"))
cox.up.DEG.ratio <- sapply(cox.eGRNs, function(x) {
  length(x$up) / length(x$genes)
})
cox.down.DEG.ratio <- sapply(cox.eGRNs, function(x) {
  length(x$down) / length(x$genes)
})
which.max(cox.up.DEG.ratio)
which.max(cox.down.DEG.ratio)


# # Select the most informative genes
# cox.overlap.genes <- sapply(cox.eGRNs, function(x) {
#   length(intersect(x$genes, dd))
# })


# Obtain the top-five TFs regulating the overrepresented DEGs in eGRNs
head(sapply(cox.eGRNs[order(cox.up.DEG.ratio, decreasing = T)], "[[", "TF")) # overrepresented DEGs
# NFIC : important in lymphoma (https://www.proteinatlas.org/ENSG00000141905-NFIC/pathology/lymphoma)
# MYCN (x 2) : important in CLL (chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.pnas.org/doi/pdf/10.1073/pnas.0304717101)
# BHLHE22 : important in lymphoma (https://www.proteinatlas.org/ENSG00000180828-BHLHE22/pathology/lymphoma)
# ATF2 : important in lymphoma (https://www.proteinatlas.org/ENSG00000115966-ATF2/pathology/lymphoma)
# NR4A1 : important in lymohoma (https://www.nature.com/articles/s41598-018-32972-4)
# GABPA (x 3) : importnt in CLL (https://maayanlab.cloud/Harmonizome/gene_set/GABPA/ENCODE+Transcription+Factor+Targets)
# KLF5 (x 2) : important in lymphoma (https://www.proteinatlas.org/ENSG00000102554-KLF5/pathology/lymphoma)
# MYOG (x 3) : important in CLL (https://maayanlab.cloud/Harmonizome/gene_set/MYOG/ENCODE+Transcription+Factor+Targets)
# SP1 (x 3) : CLL (https://www.embopress.org/doi/full/10.15252/msb.20188339)
# STAT3 (x 3) : CLL (https://pubmed.ncbi.nlm.nih.gov/32075387/)
# TCF12 (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TCF12/ENCODE+Transcription+Factor+Targets)
# TCF3 (x 3) : CLL (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8096467/)
# TFAP2C (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TFAP2C/CHEA+Transcription+Factor+Targets)
# ZNF682 (x 2) : lymphoma (https://www.proteinatlas.org/ENSG00000197124-ZNF682/pathology/lymphoma)
# EGR1 (x 2) : lymphoma (https://www.karger.com/Article/Fulltext/494867)
# SP1 (x 3) : CLL (https://www.embopress.org/doi/full/10.15252/msb.20188339)


head(sapply(cox.eGRNs[order(cox.down.DEG.ratio, decreasing = T)], "[[", "TF")) # underrepresented DEGs
# TFAP2C (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TFAP2C/ENCODE+Transcription+Factor+Targets)
# EGR1 (x 2) : lymphoma (https://www.karger.com/Article/Fulltext/494867)
# GABPA (x 3) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/GABPA/ENCODE+Transcription+Factor+Targets)
# KLF5 (x 2) : lymphoma (https://www.proteinatlas.org/ENSG00000102554-KLF5/pathology/lymphoma)
# MYCN (x 2) : CLL (chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.pnas.org/doi/pdf/10.1073/pnas.0304717101)
# MYOG (x 3) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/MYOG/ENCODE+Transcription+Factor+Targets)
# SP1 (x 3) : CLL (https://www.embopress.org/doi/full/10.15252/msb.20188339)
# STAT3 (x 3) : CLL (https://pubmed.ncbi.nlm.nih.gov/32075387/)
# TCF12 (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TCF12/ENCODE+Transcription+Factor+Targets)
# TCF3 (x 3) : CLL (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8096467/)
# TFAP2C (x 2) : CLL (https://maayanlab.cloud/Harmonizome/gene_set/TFAP2C/CHEA+Transcription+Factor+Targets)
# ZNF682 (x 2) : lymphoma (https://www.proteinatlas.org/ENSG00000197124-ZNF682/pathology/lymphoma)


# Co-factors
# MYCN : SP1
# SP1 : EGR1, GABPA, JUN, ESRRA, KLF4, MYCN, MYOG, STAT3, cox, YY1
# EGR1 : SP1, cox
# GABPA : SP1, YY1
# JUN : ATF2, FOS, KLF5, SP1, STAT3, cox
# MYOG : SP1
# STAT3 : BATF3, JUN, SP1
# BATF3 : STAT3
# TFAP2C : TP53


# Generate the UMAP associated with survival
Idents(lymph.obj) <- lymph.obj$survival
p.cox.UMAP <- get_UMAP(object = lymph.obj, reduction.method = "wnn.umap",
                   txt = "Survival")


# Build bulk expression table for GSVA analysis
source("/fs/ess/PCON0022/liyang/r_utilities/functions/transcriptome_tools.R")
cox.gene.set <- lapply(cox.eGRNs, function(x) {
  #x$genes
  x$up
  #x$down
})
names(cox.gene.set) <- seq_along(cox.gene.set)
cox.bulk <- qs::qread("/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/2_survival_GSE22762/bulk.qsave")
cox.bulk.expr <- cox.bulk$matrix
colnames(cox.bulk.expr)
dim(cox.bulk.expr)
cox.gsva.m <- run_GSVA(X = cox.bulk.expr, gs = cox.gene.set)
dim(cox.gsva.m)
cox.gsva.m[1:5, 1:5]
cox.quantiles <- get_quantiles(cox.gsva.m)
cox.bulk.meta <- cox.bulk$phenotype


# # Test which genes are overprepresented
# gene.diff <- apply(cox.gsva.m, 1, function(x) {
#   c(mean(x[long.samples]), mean(x[short.samples]))
# })


# # Make table to generate Kaplan-Meier curve for the first eGRN
# cox.KM.df <- make_DF_for_KM_curve(meta.df = cox.bulk.meta, col.names = c("time", "status"), 
#                                  quantiles = cox.quantiles[[1]])


# Generate Kaplan-Meier curve
source("/fs/ess/PCON0022/liyang/r_utilities/functions/visual_tools.R")
cox.pval <- lapply(cox.quantiles, function(quantiles) {
  cox.KM.df <- make_DF_for_KM_curve(meta.df = cox.bulk.meta, col.names = c("time", "status"),
                                   quantiles = quantiles)
  get_KM_curve(meta.df = cox.KM.df, category = "GSVA", path = NULL, draw = F, 
               legend.labs = c("high", "low"), pval = T, 
               risk.table.height = .15, 
               risk.table = T, palette = c("red", "blue"), 
               legend.title = "GSVA enrichment score")
}) %>% unlist
range(cox.pval)
cox.significant <- cox.pval[cox.pval < 0.05]
cox.significant.TFs <- sapply(irg.egrns[as.numeric(names(cox.significant))], 
                              "[[", "TF")
which.min(cox.significant)
table(cox.significant.TFs)


# Interpretations : 
# Two eGRNs : MYCN (x 2)
# Co-factors : JUN and KLF5


# MYCN 1
cox.KM.df <- make_DF_for_KM_curve(meta.df = cox.bulk.meta, col.names = c("time", "status"),
                                  quantiles = cox.quantiles[["12"]])
p.cox.MYCN1 <- get_KM_curve(meta.df = cox.KM.df, category = "GSVA", path = NULL, draw = T, 
             legend.labs = c("high", "low"), pval = F, 
             risk.table.height = .30, title = "MYCN", 
             risk.table = T, palette = c("red", "blue"), 
             legend.title = "GSVA enrichment score")


# MYCN 2
cox.KM.df <- make_DF_for_KM_curve(meta.df = cox.bulk.meta, col.names = c("time", "status"),
                                  quantiles = cox.quantiles[["42"]])
p.cox.MYCN2 <- get_KM_curve(meta.df = cox.KM.df, category = "GSVA", path = NULL, draw = T, 
                           legend.labs = c("high", "low"), pval = F, 
                           risk.table.height = .30, title = "MYCN", 
                           risk.table = T, palette = c("red", "blue"), 
                           legend.title = "GSVA enrichment score")


# JUN
cox.KM.df <- make_DF_for_KM_curve(meta.df = cox.bulk.meta, col.names = c("time", "status"),
                                  quantiles = cox.quantiles[["49"]])
p.cox.JUN <- get_KM_curve(meta.df = cox.KM.df, category = "GSVA", path = NULL, draw = T, 
                            legend.labs = c("high", "low"), pval = F, 
                            risk.table.height = .30, title = "JUN", 
                            risk.table = T, palette = c("red", "blue"), 
                            legend.title = "GSVA enrichment score")


# KLF5 : the co-factor of JUN
cox.KM.df <- make_DF_for_KM_curve(meta.df = cox.bulk.meta, col.names = c("time", "status"),
                                  quantiles = cox.quantiles[["50"]])
p.cox.KLF5 <- get_KM_curve(meta.df = cox.KM.df, category = "GSVA", path = NULL, draw = T, 
                          legend.labs = c("high", "low"), pval = F, 
                          risk.table.height = .30, title = "KLF5", 
                          risk.table = T, palette = c("red", "blue"), 
                          legend.title = "GSVA enrichment score")

p.cox.KM <- p.cox.UMAP | (p.cox.MYCN1$plot / p.cox.MYCN1$table) | 
  (p.cox.MYCN2$plot /  p.cox.MYCN2$table) | (p.cox.JUN$plot / p.cox.JUN$table) | 
  (p.cox.KLF5$plot / p.cox.KLF5$table)


##################################################################################
#                                                                                #
#                             13. Generate plots                                 #
#                                                                                #
##################################################################################


# Merge all panels
p.IRG <- ggpubr::ggarrange(p.boxplot, 
                           ggpubr::ggarrange(p.tumor, p.tumor.box, ncol = 2, labels = c("B", "C"), 
                                             widths = c(1, 2)), 
                           ggpubr::ggarrange(p.TP53, p.TP53.box, ncol = 2, labels = c("D", "E"), 
                                             widths = c(1, 2)), 
                           ggpubr::ggarrange(p.cox, p.cox.KM, ncol = 2, labels = c("F", "G"), 
                                             widths = c(1, 2)), 
          ncol = 1, heights = c(2, 1, 1, 2))
# save_image(p = p.IRG, path = paste0(work.dir, "IRG.eps"))
save_image(p = p.IRG, path = paste0(work.dir, "IRG.png"), 
           width = 6000, height = 5000)


# 1. Save the boxplots
save_image(p = p.boxplot, path = paste0(work.dir, "1_boxplots.eps"), 
           width = 3000, height = 500)


# 2. Save the UMAP plot of tumor v.s. healthy
save_image(p = p.tumor, path = paste0(work.dir, "1_UMAP_tumor.eps"))


# 3. Save the boxplots of four eGRNs between tumor v.s. healthy
save_image(p = p.tumor.box, path = paste0(work.dir, "1_boxplots_GSVA.eps"))


# 4. Save the UMAP plot of TP53 v.s. unmutated
save_image(p = p.TP53, path = paste0(work.dir, "2_UMAP_TP53.eps"))


# 5. Save the boxplots of four eGRNs between TP53 v.s. unmutated
save_image(p = p.TP53.box, path = paste0(work.dir, "2_boxplots_GSVA.eps"))


# 6. Save the UMAP plot of four eGRNs associated with survival
save_image(p = p.cox, path = paste0(work.dir, "3_UMAP_survival.eps"))


# 7. Save the Kaplan-Meier curve of four eGRNs associated with survival
save_image(p = p.cox.KM, path = paste0(work.dir, "3_KM_survival.eps"))


##################################################################################
#                                                                                #
#                                Merge all the panels                            #
#                                                                                #
##################################################################################


library(ggpubr)
# p <- (p.boxplot | p.heatmap | p.ct) / (p.phenotype.piechart | pos.bar | p.kegg.bars)
# print(p)
# save_image(p.boxplot | p.heatmap | p.ct, path = paste0(work.dir, "IRG_up.eps"))
# save_image(p.phenotype.piechart | pos.bar | p.kegg.bars, path = paste0(work.dir, "IRG_down.eps"))
# p <- ggarrange(p.boxplot, p.ct, p.phenotype, 
#           nrow = 1, labels = c("A", "B", "C"))
# p <- ggdraw() + draw_plot(p.boxplot, x = 0, y = .5, width = .5)
p.box.phenotype.ct <- egg::ggarrange(p.boxplot, p.ct, p.phenotype, heights = 1, 
                  widths = c(1, 1.5, 1.5), nrow = 1, labels = c("A", "B", "C"))
p.box.phenotype.ct
save_image(p = p.box.phenotype.ct, width = 1000, height = 500, 
           path = paste0(phenotype.dir, "boxplot_UMAP_Scissor_CT.eps"))


# Save the images of the piechart of phenotype, barplot of Scissor+ cell distribution across
# cell types, and barplot of KEGG pathway analysis
save_image(p = p.phenotype.piechart, width = 500, height = 300, 
           path = paste0(phenotype.dir, "piechart.eps"))
save_image(p = pos.bar, width = 500, height = 500, 
           path = paste0(phenotype.dir, "barplot_Scissor+_CTs.eps"))
save_image(p = p.kegg.bars, width = 500, height = 500, 
           path = paste0(phenotype.dir, "barplot_KEGG_pathways.eps"))
# p.heatmap.pie.bar.enrich <- egg::ggarrange(p.phenotype.piechart, pos.bar, p.kegg.bars, 
#                                            heights = 1, widths = c(1, 1, 1), nrow = 1,
#                                            labels = c("E", "F", "G"))
# p.heatmap.pie.bar.enrich
