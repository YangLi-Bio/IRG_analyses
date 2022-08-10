#######################################################
#                                                     #
#        Perform analyses regarding cell types        #
#                                                     #
#######################################################


# Libraries
library(Seurat)
library(Signac)
library(dplyr)
library(pbmcapply)
library(pbapply)
library(parallel)
library(data.table)
library(aricode)
library(CFTK)
library(Matrix)
library(reshape2)
library(igraph)
library(mstknnclust)
library(statGraph)
library(UpSetR)
library(ggplot2)
library(RColorBrewer)
library(scales)


# Global parameters
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Rfile/"
coherence.mpx <- 1.0 # multiple of ATAC ratio
distance <- 250000


# Source scripts
source(paste0(tool.dir, "multiome_tools.R"))


# Load files
obj.celltyped <- qs::qread(paste0(R.dir, "Obj_celltyped.qsave"))
exp.sum <- qs::qread(paste0(R.dir, "eGRN_exp_sum.qsave"))
acc.sum <- qs::qread(paste0(R.dir, "eGRN_acc_sum.qsave"))
eGRN.ct <- qs::qread(paste0(R.dir, "eGRN_active_cells_in_CT.qsave"))
occur.m <- qs::qread(paste0(R.dir, "eGRN_occur_matrix.qsave"))
eGRNs <- qs::qread(paste0(R.dir, "Extended_eGRNs.qsave"))
eGRN.ct <- qs::qread(paste0(R.dir, "eGRN_active_cells_in_CT.qsave"))


# Prepare matrix
range(occur.m)
dim(occur.m)
eGRN.cells <- Reduce("union", lapply(eGRNs, "[[", "cells"))
length(eGRN.cells)
obj.celltyped$cell.type
subset.cells <- obj.celltyped$cell.type[eGRN.cells]
length(subset.cells)
subset.m <- occur.m[, names(subset.cells), drop = F]
dim(subset.m)
min.df <- do.call("rbind", pbmclapply(2 : ncol(subset.m), function(i) {
  i.set <- which(subset.m[, i] > 0)
  Reduce("rbind", lapply(1:(i - 1), function(j) {
    j.set <- which(subset.m[, j] > 0)
    ij.jacc <- length(intersect(i.set, j.set)) / 
      min(length(i.set), length(j.set))
    c(i, j, ij.jacc)
  }))
}, mc.cores = detectCores()))
dim(min.df)


# Cell clustering
jacc.df <- min.df
jacc.total <- rbind(jacc.df, Reduce("rbind", pbmclapply(1:ncol(subset.m), function(i) {
  c(i, i, 1)
}, mc.cores = detectCores())))
colnames(jacc.total) <- c("i", "j", "x")
jacc.total <- as.data.frame(jacc.total)
head(jacc.total)
regularMatrix <- acast(jacc.total, i ~ j, value.var = "x")
class(regularMatrix)
dim(regularMatrix)
regularMatrix[1:3, 1:3]
regularMatrix[is.na(regularMatrix)] <- 1
regularMatrix[1:3, 1:3]
dim(regularMatrix)
rownames(regularMatrix) <- names(subset.cells)
colnames(regularMatrix) <- names(subset.cells)
regularMatrix[1:3, 1:3]
apply(regularMatrix, 1, sum) %>% range
apply(regularMatrix, 2, sum) %>% range
regularMatrix[1:3, 1:3]
disRegularMatrix <- 1 - regularMatrix
disRegularMatrix[1:3, 1:3]
disRegularMatrix <- disRegularMatrix / apply(disRegularMatrix, 1, sum)
distanceMatrix <- as.dist(disRegularMatrix)
length(distanceMatrix)
class(distanceMatrix)
distanceMatrix[1:3]
distanceMatrix[is.na(distanceMatrix)] <- 0
clusters <- hclust(distanceMatrix, method = "ward.D2")
plot(clusters)
length(unique(subset.cells))


# Calculate ARI
group <- cutree(clusters, k = 11)
group
head(group)
head(subset.cells)
jacc.ari <- ARI(group, subset.cells)
jacc.ari # 0.329393
# ARI will be mentioned in the main text without using plots.


# Showcase the distribution of eGRN-active cell subpopulations within
# cell types using Upset diagram and a barplot (logarithm p-values)
# weblink: https://upset.app/advanced/


# Build the nested list to generate upset diagram
class(eGRN.ct)
eGRN.ct[1:3, 1:3]
dim(eGRN.ct)
ct.sum <- apply(eGRN.ct, 2, sum)
pval.dt <- rbindlist(pbmclapply(1:nrow(eGRN.ct), function(i) {
  message(i, "\n")
  group1 <- length(eGRNs[[i]]$cells)
  lapply(seq_along(eGRN.ct[i,]), function(j) {
    overlap <- as.matrix(eGRN.ct)[i, j]
    group2 <- ct.sum[j]
    total <- ncol(obj.celltyped)
    pval <- phyper(overlap - 1, group2, total - group2, group1, lower.tail = FALSE)
    pval
  })
}, mc.cores = detectCores()))
dim(pval.dt)
pval.dt[1:3, 1:3]
colnames(pval.dt) <- colnames(eGRN.ct)
pval.dt[1:3, 1:3]
padj.dt <- pval.dt * ncol(eGRN.ct) # Bonferroni correction
range(pval.dt)
range(padj.dt)


# Select cell-type-specific eGRNs
ct.spec.dt <- padj.dt < 0.05
ct.spec.dt[1:3, 1:3]
apply(ct.spec.dt, 1, sum)
rownames(ct.spec.dt) <- rownames(padj.dt)
cts.eGRN.ids <- which(apply(ct.spec.dt, 1, sum) == 1)
cts.eGRN.ids # 16
general.eGRN.ids <- which(apply(ct.spec.dt, 1, sum) > 1)
general.eGRN.ids # 32
length(general.eGRN.ids) + length(cts.eGRN.ids)
save.image("Visual.RData")
dim(ct.spec.dt)


# Prepare the list to generate upset diagram
# http://vcg.github.io/upset/
ct.eGRN.list <- pbmclapply(1:ncol(ct.spec.dt), function(i) {
  rownames(ct.spec.dt)[which(ct.spec.dt[, i] > 0)]
}, mc.cores = detectCores())
length(ct.eGRN.list)
names(ct.eGRN.list) <- colnames(ct.spec.dt)
head(ct.eGRN.list)
upset(fromList(ct.eGRN.list), order.by = "freq", nsets = 11, 
      nintersects = NA)
# Upset graph is not suitable for this result. I will use dot plot instead.


# Generate UMAP plot
Idents(obj.celltyped) <- obj.celltyped$cell.type
levels(obj.celltyped$cell.type)
obj.celltyped[["abbreviations"]] <- obj.celltyped$cell.type
levels(obj.celltyped$abbreviations)
names(celltype.abbr)
celltype.abbr
levels(obj.celltyped$abbreviations) <- celltype.abbr
levels(obj.celltyped$abbreviations)
# col.celltypes <- brewer.pal(n = 11, name = "Spectral")
col.celltypes <- c(
  "#FF0000", 
  "#FEB139",
  "#FAD9A1",
  "#66BFBF",
  "#BD4291",
  "#187498",
  "#143F6B",
  "#357C3C",
  "#FF0075",
  "#30475E",
  "#FFB400"
  )
names(col.celltypes) <- levels(obj.celltyped$abbreviations)
col.celltypes
length(col.celltypes)
Idents(obj.celltyped) <- obj.celltyped$abbreviations
p.celltype.umap <- DimPlot(obj.celltyped, reduction = "wnn.umap", 
                           cols = col.celltypes)
p.celltype.umap


# Dotplot
# Rows (11): cell types
# Columns (50): eGRNs
# Dot color: p-value
# Dot size: overlaps


# Build the data frame to generate dot plot
ct.spec.dt[1:3, 1:3]
colnames(ct.spec.dt) <- as.character(levels(obj.celltyped$abbreviations))
as.data.frame(ct.spec.dt)[1:3, 1:3]
ct.spec.dt <- as.data.frame(ct.spec.dt)
pval.dt <- as.data.frame(pval.dt)
eGRN.ct <- as.data.frame(eGRN.ct)
eGRN.p.ratio.df <- do.call("rbind", pbmclapply(1:nrow(ct.spec.dt), function(i) {
  message (i, "\n")
  Reduce("rbind", lapply(1:ncol(ct.spec.dt), function(j) {
    if (ct.spec.dt[i, j]) {
      message (j, "\n")
      c(eGRN = rownames(ct.spec.dt)[i], 
        Cell.type = colnames(ct.spec.dt)[j], 
        `log p-value` = -log(pval.dt[i, j]), 
        `% Overlap` = 100 * eGRN.ct[i, j] / 
          length(eGRNs[[i]]$cells))
    }
  }))
}, mc.cores = detectCores()))
# , mc.cores = detectCores()
dim(eGRN.p.ratio.df)
head(eGRN.p.ratio.df)


# Generate the dot plot
class(eGRN.p.ratio.df)
eGRN.p.ratio.df[1:3, 1:3]
eGRN.p.ratio.df <- as.data.frame(eGRN.p.ratio.df)
rownames(eGRN.p.ratio.df) <- NULL
eGRN.p.ratio.df$`log p-value` <- as.numeric(eGRN.p.ratio.df$`log p-value`)
eGRN.p.ratio.df$`% Overlap` <- as.numeric(eGRN.p.ratio.df$`% Overlap`)
head(eGRN.p.ratio.df)
rownames(eGRN.p.ratio.df)
head(eGRN.p.ratio.df)
eGRN.p.ratio.df$eGRN <- factor(eGRN.p.ratio.df$eGRN, levels = unique(eGRN.p.ratio.df$eGRN))
p.dot.eGRN.CT.pval.ratio <- ggplot(eGRN.p.ratio.df, aes(x = eGRN, y = Cell.type, 
                            fill = `log p-value`, size = `% Overlap`)) + 
  geom_point(shape = 21) + 
  scale_fill_gradient(low = "#FFDEDE", high = "red", name = 'log p-value') + 
  cowplot::theme_cowplot() + 
  # theme(axis.line  = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
p.dot.eGRN.CT.pval.ratio
save_image(p = p.dot.eGRN.CT.pval.ratio, path = paste0(image.dir, "Dotplot_eGRN_celltype.eps"))


# Cell type distribuiton among IgHV v.s. unmutated
ighv.scissor.cells <- qs::qread("/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/1_classification_GSE113386/GSE113386_single_param.qsave")
sapply(ighv.scissor.cells[3:4], length)

# Scissor+ cells
cells.ighv.plus <- ighv.scissor.cells$Scissor_pos
length(cells.ighv.plus)
celltypes.ighv.plus <- obj.celltyped$cell.type[intersect(names(obj.celltyped$cell.type), 
                                                         cells.ighv.plus)]
table(celltypes.ighv.plus)
celltype.ratio.ighv.plus <- table(celltypes.ighv.plus) / sum(table(celltypes.ighv.plus))
celltype.ratio.ighv.plus <- sort(celltype.ratio.ighv.plus, decreasing = T)
head(celltype.ratio.ighv.plus)


# Scissor- cells
cells.ighv.minus <- ighv.scissor.cells$Scissor_neg
length(cells.ighv.minus)
celltypes.ighv.minus <- obj.celltyped$cell.type[intersect(names(obj.celltyped$cell.type), 
                                                         cells.ighv.minus)]
table(celltypes.ighv.minus)
celltype.ratio.ighv.minus <- table(celltypes.ighv.minus) / sum(table(celltypes.ighv.minus))
celltype.ratio.ighv.minus <- sort(celltype.ratio.ighv.minus, decreasing = T)
head(celltype.ratio.ighv.minus)


# Generate piechart and two barplots
length(cells.ighv.plus)
length(cells.ighv.minus)
ighv.pie.df <- data.frame(group = c("IgHV", "unmutated"), 
                     value = c(length(cells.ighv.plus), length(cells.ighv.minus)))
ighv.pie.df
ighv.pie.df <- cbind(ighv.pie.df, count = ighv.pie.df$value)
ighv.pie.df$value <- ighv.pie.df$value / sum(ighv.pie.df$value)
ighv.pie.df
ighv.pie.df$value <- round(ighv.pie.df$value, digits = 3)
ighv.pie.df
p.pie.ighv <- ggplot(ighv.pie.df, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 10, color = "white") +
  coord_polar("y", start = 0) + 
  # geom_label(aes(label = group),
  #             position = position_stack(vjust = 0.5)) +
  theme_void() + theme(legend.position = "none") + 
  theme(legend.position="none") +
  # geom_text(aes(y = value, label = paste0(group, "\n", count, "\n(", 100 * value, "%)")), 
  #           color = "white", size = 5, position = position_dodge(.9)) +
  scale_fill_manual(values = c("Black", "#413F42"))
p.pie.ighv
save_image(p = p.pie.ighv, path = paste0(image.dir, "Pie_IgHV.eps"))


# Barplots for Scissor+ cells of IgHV
col.celltypes
celltype.ratio.ighv.plus
names(celltype.ratio.ighv.plus) <- names(col.celltypes)
levels(celltypes.ighv.plus) <- names(col.celltypes)
names(col.celltypes)
celltype.abbr <- c(
  "Tumor B",
  "Exh CD4+ T", 
  "Exh CD8+ T",
  "Cent mem CD8+ T",
  "Effect CD4+ T 1", 
  "TAMs",
  "Effect CD4+ T 2",
  "Effect CD8+ T",
  "Normal B", 
  "Tumor B prolifer",
  "DCs"
)
names(celltype.abbr) <- names(col.celltypes)
sum(celltype.ratio.ighv.plus[1:3])
ighv.plus.cells <- table(celltypes.ighv.plus)
ighv.plus.cells <- sort(ighv.plus.cells, decreasing = T)
head(ighv.plus.cells)
ighv.plus.df <- as.data.frame(celltype.ratio.ighv.plus[1:3])
ighv.plus.df <- cbind(ighv.plus.df, ighv.plus.cells[1:3])
ighv.plus.df <- ighv.plus.df[, -3]
ighv.plus.df <- cbind(as.character(ighv.plus.df[, 1]), ighv.plus.df)
rownames(ighv.plus.df) <- NULL
colnames(ighv.plus.df) <- c("Abbreviation", "name", "value", "count")
ighv.plus.df$Abbreviation <- factor(ighv.plus.df$Abbreviation, 
                                    levels = celltype.abbr[names(celltype.ratio.ighv.plus[1:3])])
ighv.plus.df
ighv.plus.df$Abbreviation
p.bar.ighv.plus <- ggplot(data = ighv.plus.df, aes(x = Abbreviation,
                                                   y = count, 
                                                   fill = name)) + 
  geom_bar(position = 'dodge', stat = 'identity', color = "black") + 
  scale_fill_manual(values = col.celltypes[as.character(ighv.plus.df$name)]) + 
  geom_text(aes(label = paste0(count, "\n(", round(value * 100, digits = 1), "%)")), 
            position = position_dodge(width = 0.9), vjust = -0.50) + 
  theme_classic() + theme(legend.position = "none", axis.title.x = element_blank(), 
                          axis.text.x = element_text(angle = 45, colour = "black", vjust = 0.50, size = 10),
                          axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                          axis.text.y = element_text(colour = "black", size = 10)) + 
  labs(y = "Number of cells") + ylim(c(0, 800))
p.bar.ighv.plus
save_image(p = p.bar.ighv.plus, path = paste0(image.dir, "Barplot_IgHV_plus.eps"))


# Barplot for Scissor- cells in IgHV
ighv.minus.cells <- table(celltypes.ighv.minus)
ighv.minus.cells <- sort(ighv.minus.cells, decreasing = T)
head(ighv.minus.cells)
ighv.minus.df <- as.data.frame(celltype.ratio.ighv.minus[1:3])
ighv.minus.df <- cbind(ighv.minus.df, ighv.minus.cells[1:3])
ighv.minus.df <- ighv.minus.df[, -3]
ighv.minus.df <- cbind(celltype.abbr[as.character(ighv.minus.df[, 1])], ighv.minus.df)
rownames(ighv.minus.df) <- NULL
dim(ighv.minus.df)
colnames(ighv.minus.df) <- c("Abbreviation", "name", "value", "count")
ighv.minus.df$Abbreviation <- factor(ighv.minus.df$Abbreviation, 
                                    levels = celltype.abbr[names(celltype.ratio.ighv.minus[1:3])])
ighv.minus.df
ighv.minus.df$Abbreviation
p.bar.ighv.minus <- ggplot(data = ighv.minus.df, aes(x = Abbreviation,
                                                   y = count, 
                                                   fill = name)) + 
  geom_bar(position = 'dodge', stat = 'identity', color = "black") + 
  scale_fill_manual(values = col.celltypes[as.character(ighv.minus.df$name)]) + 
  geom_text(aes(label = paste0(count, "\n(", round(value * 100, digits = 1), "%)")), 
            position = position_dodge(width = 0.9), vjust = -0.50) + 
  theme_classic() + theme(legend.position = "none", axis.title.x = element_blank(), 
                          axis.text.x = element_text(angle = 45, colour = "black", vjust = 0.50),
                          axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                          axis.text.y = element_text(colour = "black")) + 
  labs(y = "Number of cells") + ylim(c(0, 60))
p.bar.ighv.minus
save_image(p = p.bar.ighv.minus, path = paste0(image.dir, "Barplot_IgHV_minus.eps"))


# Cell type distribution among TP53 v.s. unmutated
tp53.scissor.cells <- qs::qread("/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/1_classification_GSE113386/GSE113386_single_param.qsave")
sapply(tp53.scissor.cells[3:4], length)

# Scissor+ cells
cells.tp53.plus <- tp53.scissor.cells$Scissor_pos
length(cells.tp53.plus)
celltypes.tp53.plus <- obj.celltyped$cell.type[intersect(names(obj.celltyped$cell.type), 
                                                         cells.tp53.plus)]
table(celltypes.tp53.plus)
celltype.ratio.tp53.plus <- table(celltypes.tp53.plus) / sum(table(celltypes.tp53.plus))
celltype.ratio.tp53.plus <- sort(celltype.ratio.tp53.plus, decreasing = T)
head(celltype.ratio.tp53.plus)


# Scissor- cells
cells.tp53.minus <- tp53.scissor.cells$Scissor_neg
length(cells.tp53.minus)
celltypes.tp53.minus <- obj.celltyped$cell.type[intersect(names(obj.celltyped$cell.type), 
                                                          cells.tp53.minus)]
table(celltypes.tp53.minus)
celltype.ratio.tp53.minus <- table(celltypes.tp53.minus) / sum(table(celltypes.tp53.minus))
celltype.ratio.tp53.minus <- sort(celltype.ratio.tp53.minus, decreasing = T)
head(celltype.ratio.tp53.minus)


# Generate piechart and two barplots
length(cells.tp53.plus)
length(cells.tp53.minus)
tp53.pie.df <- data.frame(group = c("tp53", "unmutated"), 
                     value = c(length(cells.tp53.plus), length(cells.tp53.minus)))
tp53.pie.df
tp53.pie.df$value <- tp53.pie.df$value / sum(tp53.pie.df$value)
tp53.pie.df
p.pie.tp53 <- ggplot(tp53.pie.df, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) + 
  # geom_label(aes(label = group), 
  #             position = position_stack(vjust = 0.5)) + 
  theme_void() + theme(legend.position = "none") + 
  # theme(legend.position="none") +
  # geom_text(aes(y = value, label = value), color = "white", size = 10) +
  scale_fill_manual(values = c("Black", "#413F42"))
p.pie.tp53
save_image(p = p.pie.tp53, path = paste0(image.dir, "Pie_TP53.eps"))


# Barplots for Scissor+ cells of tp53
sum(celltype.ratio.tp53.plus[1:3])
tp53.plus.cells <- table(celltypes.tp53.plus)
tp53.plus.cells <- sort(tp53.plus.cells, decreasing = T)
head(tp53.plus.cells)
tp53.plus.df <- as.data.frame(celltype.ratio.tp53.plus[1:3])
tp53.plus.df <- cbind(tp53.plus.df, tp53.plus.cells[1:3])
tp53.plus.df <- tp53.plus.df[, -3]
tp53.plus.df <- cbind(celltype.abbr[as.character(tp53.plus.df[, 1])], tp53.plus.df)
rownames(tp53.plus.df) <- NULL
colnames(tp53.plus.df) <- c("Abbreviation", "name", "value", "count")
tp53.plus.df$Abbreviation <- factor(tp53.plus.df$Abbreviation, 
                                    levels = celltype.abbr[names(celltype.ratio.tp53.plus[1:3])])
tp53.plus.df
tp53.plus.df$Abbreviation
p.bar.tp53.plus <- ggplot(data = tp53.plus.df, aes(x = Abbreviation,
                                                   y = count, 
                                                   fill = name)) + 
  geom_bar(position = 'dodge', stat = 'identity', color = "black") + 
  scale_fill_manual(values = col.celltypes[as.character(tp53.plus.df$name)]) + 
  geom_text(aes(label = paste0(count, "\n(", round(value * 100, digits = 1), "%)")), 
            position = position_dodge(width = 0.9), vjust = -0.50) + 
  theme_classic() + theme(legend.position = "none", axis.title.x = element_blank(), 
                          axis.text.x = element_text(angle = 45, colour = "black", vjust = 0.50),
                          axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                          axis.text.y = element_text(colour = "black")) + 
  labs(y = "Number of cells") + ylim(c(0, 700))
p.bar.tp53.plus
save_image(p = p.bar.tp53.plus, path = paste0(image.dir, "Barplot_TP53_plus.eps"))


# Barplot for Scissor- cells in tp53
tp53.minus.cells <- table(celltypes.tp53.minus)
tp53.minus.cells <- sort(tp53.minus.cells, decreasing = T)
head(tp53.minus.cells)
tp53.minus.df <- as.data.frame(celltype.ratio.tp53.minus[1:4])
sum(celltype.ratio.tp53.minus[1:4])
tp53.minus.df <- cbind(tp53.minus.df, tp53.minus.cells[1:4])
tp53.minus.df <- tp53.minus.df[, -3]
tp53.minus.df <- cbind(celltype.abbr[as.character(tp53.minus.df[, 1])], tp53.minus.df)
rownames(tp53.minus.df) <- NULL
dim(tp53.minus.df)
colnames(tp53.minus.df) <- c("Abbreviation", "name", "value", "count")
tp53.minus.df$Abbreviation <- factor(tp53.minus.df$Abbreviation, 
                                     levels = celltype.abbr[names(celltype.ratio.tp53.minus[1:4])])
tp53.minus.df
tp53.minus.df$Abbreviation
p.bar.tp53.minus <- ggplot(data = tp53.minus.df, aes(x = Abbreviation,
                                                     y = count, 
                                                     fill = name)) + 
  geom_bar(position = 'dodge', stat = 'identity', color = "black") + 
  scale_fill_manual(values = col.celltypes[as.character(tp53.minus.df$name)]) + 
  geom_text(aes(label = paste0(count, "\n(", round(value * 100, digits = 1), "%)")), 
            position = position_dodge(width = 0.9), vjust = -0.50) + 
  theme_classic() + theme(legend.position = "none", axis.title.x = element_blank(), 
                          axis.text.x = element_text(angle = 45, colour = "black", vjust = 0.50),
                          axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                          axis.text.y = element_text(colour = "black")) + 
  labs(y = "Number of cells") + ylim(c(0, 300))
p.bar.tp53.minus
save_image(p = p.bar.tp53.minus, path = paste0(image.dir, "Barplot_TP53_minus.eps"))


# Cell type distribution among survival- worse v.s. good
cox.scissor.cells <- qs::qread("/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/2_survival_GSE22762/2_survival_GSE22762_optimal_alpha_0.017.qsave")
sapply(cox.scissor.cells[3:4], length)

# Scissor+ cells
cells.cox.plus <- cox.scissor.cells$Scissor_pos
length(cells.cox.plus)
celltypes.cox.plus <- obj.celltyped$cell.type[intersect(names(obj.celltyped$cell.type), 
                                                         cells.cox.plus)]
table(celltypes.cox.plus)
celltype.ratio.cox.plus <- table(celltypes.cox.plus) / sum(table(celltypes.cox.plus))
celltype.ratio.cox.plus <- sort(celltype.ratio.cox.plus, decreasing = T)
head(celltype.ratio.cox.plus)


# Scissor- cells
cells.cox.minus <- cox.scissor.cells$Scissor_neg
length(cells.cox.minus)
celltypes.cox.minus <- obj.celltyped$cell.type[intersect(names(obj.celltyped$cell.type), 
                                                          cells.cox.minus)]
table(celltypes.cox.minus)
celltype.ratio.cox.minus <- table(celltypes.cox.minus) / sum(table(celltypes.cox.minus))
celltype.ratio.cox.minus <- sort(celltype.ratio.cox.minus, decreasing = T)
head(celltype.ratio.cox.minus)


# Generate piechart and two barplots
length(cells.cox.plus)
length(cells.cox.minus)
cox.pie.df <- data.frame(group = c("cox", "unmutated"), 
                          value = c(length(cells.cox.plus), length(cells.cox.minus)))
cox.pie.df
cox.pie.df$value <- cox.pie.df$value / sum(cox.pie.df$value)
cox.pie.df
p.pie.cox <- ggplot(cox.pie.df, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) + 
  # geom_label(aes(label = group), 
  #             position = position_stack(vjust = 0.5)) + 
  theme_void() + theme(legend.position = "none") + 
  # theme(legend.position="none") +
  # geom_text(aes(y = value, label = value), color = "white", size = 10) +
  scale_fill_manual(values = c("Black", "#413F42"))
p.pie.cox
save_image(p = p.pie.cox, path = paste0(image.dir, "Pie_Cox.eps"))


# Barplots for Scissor+ cells of cox
sum(celltype.ratio.cox.plus[1:3])
cox.plus.cells <- table(celltypes.cox.plus)
cox.plus.cells <- sort(cox.plus.cells, decreasing = T)
head(cox.plus.cells)
cox.plus.df <- as.data.frame(celltype.ratio.cox.plus[1:3])
cox.plus.df <- cbind(cox.plus.df, cox.plus.cells[1:3])
cox.plus.df <- cox.plus.df[, -3]
cox.plus.df <- cbind(celltype.abbr[as.character(cox.plus.df[, 1])], cox.plus.df)
rownames(cox.plus.df) <- NULL
colnames(cox.plus.df) <- c("Abbreviation", "name", "value", "count")
cox.plus.df$Abbreviation <- factor(cox.plus.df$Abbreviation, 
                                    levels = celltype.abbr[names(celltype.ratio.cox.plus[1:3])])
cox.plus.df
cox.plus.df$Abbreviation
range(cox.plus.df$count)
p.bar.cox.plus <- ggplot(data = cox.plus.df, aes(x = Abbreviation,
                                                   y = count, 
                                                   fill = name)) + 
  geom_bar(position = 'dodge', stat = 'identity', color = "black") + 
  scale_fill_manual(values = col.celltypes[as.character(cox.plus.df$name)]) + 
  geom_text(aes(label = paste0(count, "\n(", round(value * 100, digits = 1), "%)")), 
            position = position_dodge(width = 0.9), vjust = -0.50) + 
  theme_classic() + theme(legend.position = "none", axis.title.x = element_blank(), 
                          axis.text.x = element_text(angle = 45, colour = "black", vjust = 0.50),
                          axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                          axis.text.y = element_text(colour = "black")) + 
  labs(y = "Number of cells") + ylim(c(0, 620))
p.bar.cox.plus
save_image(p = p.bar.cox.plus, path = paste0(image.dir, "Barplot_cox_worse.eps"))


# Barplot for Scissor- cells in cox
cox.minus.cells <- table(celltypes.cox.minus)
cox.minus.cells <- sort(cox.minus.cells, decreasing = T)
head(cox.minus.cells)
sum(celltype.ratio.cox.minus[1:3])
cox.minus.df <- as.data.frame(celltype.ratio.cox.minus[1:3])
sum(celltype.ratio.cox.minus[1:3])
cox.minus.df <- cbind(cox.minus.df, cox.minus.cells[1:3])
cox.minus.df <- cox.minus.df[, -3]
cox.minus.df <- cbind(celltype.abbr[as.character(cox.minus.df[, 1])], cox.minus.df)
rownames(cox.minus.df) <- NULL
dim(cox.minus.df)
colnames(cox.minus.df) <- c("Abbreviation", "name", "value", "count")
cox.minus.df$Abbreviation <- factor(cox.minus.df$Abbreviation, 
                                     levels = celltype.abbr[names(celltype.ratio.cox.minus[1:3])])
cox.minus.df
cox.minus.df$Abbreviation
range(cox.minus.df$count)
p.bar.cox.minus <- ggplot(data = cox.minus.df, aes(x = Abbreviation,
                                                     y = count, 
                                                     fill = name)) + 
  geom_bar(position = 'dodge', stat = 'identity', color = "black") + 
  scale_fill_manual(values = col.celltypes[as.character(cox.minus.df$name)]) + 
  geom_text(aes(label = paste0(count, "\n(", round(value * 100, digits = 1), "%)")), 
            position = position_dodge(width = 0.9), vjust = -0.50) + 
  theme_classic() + theme(legend.position = "none", axis.title.x = element_blank(), 
                          axis.text.x = element_text(angle = 45, colour = "black", vjust = 0.50),
                          axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                          axis.text.y = element_text(colour = "black")) + 
  labs(y = "Number of cells") + ylim(c(0, 1050))
p.bar.cox.minus
save_image(p = p.bar.cox.minus, path = paste0(image.dir, "Barplot_cox_good.eps"))


# Determine the panels to showcase in Figure 4
# UMAP plot
# Three phenotypes: piechart, two barplots
# One dotplot
