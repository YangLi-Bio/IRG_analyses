#######################################################
#                                                     #
#   Annotate cell types based on Tong's annotation    #
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


# Global parameters
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Rfile/"
coherence.mpx <- 1.0 # multiple of ATAC ratio
distance <- 250000


# Source scripts
source(paste0(tool.dir, "multiome_tools.R"))


# # Load the object
# obj <- qs::qread(paste0(R.dir, "lymph_obj.qsave"))
# dim(obj)
# names(obj@assays)
# dim(obj[['ATAC']])
# dim(obj[['RNA']])
# dim(obj[['SCT']])
eGRNs <- qs::qread(paste0(R.dir, "Extended_eGRNs.qsave"))
# 
# 
# # Load the cell types from Tong Xiao
# ct.vec <- readRDS(paste0(R.dir, "lymphoma_CT.RDS"))
# length(ct.vec)
# table(ct.vec)
# length(table(ct.vec))
# ncol(obj)
# length(intersect(names(ct.vec), colnames(obj))) == length(ct.vec)
# obj.clustered <- cluster_multiome(pbmc = obj)
# qs::qsave(obj.clustered, paste0(R.dir, "obj_clustered.qsave"))
# 
# 
# # Calculate the cell type maximally overlapped with annotated cell type
# Idents(obj.clustered) <- obj.clustered$seurat_clusters
# cluster.ct.dt <- rbindlist(pbmclapply(levels(obj.clustered$seurat_clusters), function(i) {
#   cells.cluster <- colnames(obj.clustered)[which(Idents(obj.clustered) == i)]
#   overlapped <- intersect(cells.cluster, names(ct.vec))
#   as.list(table(ct.vec[overlapped]))
# }, mc.cores = detectCores()))
# rownames(cluster.ct.dt) <- levels(obj.clustered$seurat_clusters)
# dim(cluster.ct.dt)
# head(cluster.ct.dt)
# cluster.ct.ratio <- cluster.ct.dt / apply(cluster.ct.dt, 1, sum)
# rownames(cluster.ct.ratio) <- levels(obj.clustered$seurat_clusters)
# dim(cluster.ct.ratio)
# head(cluster.ct.ratio)
# qs::qsave(cluster.ct.ratio, paste0(R.dir, "Cluster_cell_type_ratios.qsave"))
# 
# 
# # Check the cell clusters
# apply(cluster.ct.ratio, 1, max) %>% range
# apply(cluster.ct.ratio, 1, max)
# which(apply(cluster.ct.ratio, 1, max) <= 0.50) %>% length
# 
# 
# # Perform cell subclustering for cluster IDs 2, 3, 5, 13 (not levels)
# ids.subclustered <- levels(obj.clustered$seurat_clusters)[c(2, 3, 5, 13)]
# ids.subclustered
# obj.subclustered <- FindSubCluster(
#   obj.clustered,
#   graph.name = "wsnn",
#   cluster = ids.subclustered, 
#   resolution = 0.1
# )
# table(obj.subclustered$sub.cluster)
# 
# 
# # Cell cubclusters
# Idents(obj.subclustered) <- obj.subclustered$sub.cluster
# subcluster.ct.dt <- rbindlist(pbmclapply(unique(obj.subclustered$sub.cluster), function(i) {
#   cells.cluster <- colnames(obj.subclustered)[which(Idents(obj.subclustered) == i)]
#   overlapped <- intersect(cells.cluster, names(ct.vec))
#   as.list(table(ct.vec[overlapped]))
# }, mc.cores = detectCores()))
# rownames(subcluster.ct.dt) <- unique(obj.subclustered$sub.cluster)
# rownames(subcluster.ct.dt)
# dim(subcluster.ct.dt)
# head(subcluster.ct.dt)
# subcluster.ct.ratio <- subcluster.ct.dt / apply(subcluster.ct.dt, 1, sum)
# rownames(subcluster.ct.ratio) <- rownames(subcluster.ct.dt)
# dim(subcluster.ct.ratio)
# head(subcluster.ct.ratio)
# qs::qsave(subcluster.ct.ratio, paste0(R.dir, "Subcluster_cell_type_ratios.qsave"))
# 
# 
# # Check the cell subclusters
# apply(subcluster.ct.ratio, 1, max) %>% range
# apply(subcluster.ct.ratio, 1, max)
# rownames(subcluster.ct.ratio)
# which(apply(subcluster.ct.ratio, 1, max) <= 0.50) %>% length
# 
# 
# # Annotate cell types according to cell clusters
# cluster.ct.ratio
# apply(cluster.ct.ratio, 1, which.max) %>% unique %>% length
# ct.list <- names(table(ct.vec))
# ct.list
# levels(obj.subclustered$seurat_clusters)
# cluster.ct.map <- colnames(cluster.ct.ratio)[apply(cluster.ct.ratio, 1, which.max)]
# names(cluster.ct.map) <- rownames(cluster.ct.ratio)
# head(cluster.ct.map)
# obj.celltyped <- obj.subclustered
# obj.celltyped[['cell.type']] <- obj.celltyped$seurat_clusters
# obj.celltyped$cell.type
# levels(obj.celltyped$cell.type) <- cluster.ct.map
# levels(obj.celltyped$cell.type)
# qs::qsave(obj.celltyped, paste0(R.dir, "Obj_celltyped.qsave"))
# 
# 
# # Compare the eGRN-active cell subpopulation with cell types
# sapply(eGRNs, "[[", "cells") %>% sapply(., length) %>% range
# eGRN.ct <- rbindlist(pbmclapply(eGRNs, function(x) {
#   x.ct <- intersect(x$cells, names(obj.celltyped$cell.type))
#   table(obj.celltyped$cell.type[x.ct]) %>% as.list
# }, mc.cores = detectCores()))
# qs::qsave(eGRN.ct, paste0(R.dir, "eGRN_active_cells_in_CT.qsave"))
# ct.cover <- apply(eGRN.ct > 0, 1, sum)
# table(ct.cover)




# Compare eGRN-active cell subpopulations with cell types
eGRN.ct <- qs::qread(paste0(R.dir, "eGRN_active_cells_in_CT.qsave"))
obj.celltyped <- qs::qread(paste0(R.dir, "Obj_celltyped.qsave"))


eGRN.ct[1:3, 1:3]
eGRN.ratio <- eGRN.ct / apply(eGRN.ct, 1, sum)
eGRN.ratio[1:3, 1:3]
apply(eGRN.ratio, 1, max) %>% range
overlap.vec <- apply(eGRN.ct, 1, max)
ct.sum <- apply(eGRN.ct, 2, sum)
# To-do: generate a stacked barplots


# Calculate the hypergeometric p-value for overlaps between eGRN-active cell subpopulations
# and cell types
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
cts.eGRN.ids <- which(apply(ct.spec.dt, 1, sum) == 1)
cts.eGRN.ids
general.eGRN.ids <- setdiff(seq_along(eGRNs), cts.eGRN.ids)
general.eGRN.ids
length(general.eGRN.ids) + length(cts.eGRN.ids)


# Build matrices for cell clustering based on eGRNs
rna.m <- GetAssayData(obj.celltyped, slot = "data", assay = "RNA")
atac.m <- GetAssayData(obj.celltyped, slot = "data", assay = "ATAC")
exp.m <- rbindlist(pblapply(eGRNs, function(x) {
  mclapply(colnames(obj.celltyped), function(y) {
    mean(rna.m[x$genes, y])
  }, mc.cores = detectCores())
}))
acc.m <- rbindlist(pblapply(eGRNs, function(x) {
  mclapply(colnames(obj.celltyped), function(y) {
    mean(atac.m[x$peaks, y])
  }, mc.cores = detectCores())
}))

# exp.sum <- rbindlist(pblapply(eGRNs, function(x) {
#   mclapply(colnames(obj.celltyped), function(y) {
#     sum(rna.m[x$genes, y])
#   }, mc.cores = detectCores())
# }))
# acc.sum <- rbindlist(pblapply(eGRNs, function(x) {
#   mclapply(colnames(obj.celltyped), function(y) {
#     sum(atac.m[x$peaks, y])
#   }, mc.cores = detectCores())
# }))
# 


qs::qsave(exp.m, paste0(R.dir, "eGRN_expression.qsave"))
qs::qsave(acc.m, paste0(R.dir, "eGRN_accessibility.qsave"))


exp.sum <- Reduce("rbind", pbmclapply(eGRNs, function(x) {
  apply(rna.m[x$genes,], 2, sum)
}, mc.cores = detectCores()))
acc.sum <- Reduce("rbind", pbmclapply(eGRNs, function(x) {
  apply(atac.m[x$peaks,], 2, sum)
}, mc.cores = detectCores()))
rownames(exp.sum) <- 1:nrow(exp.sum)
rownames(acc.sum) <- 1:nrow(acc.sum)
qs::qsave(exp.sum, paste0(R.dir, "eGRN_exp_sum.qsave"))
qs::qsave(acc.sum, paste0(R.dir, "eGRN_acc_sum.qsave"))
dim(exp.sum)
dim(acc.sum)


# Normalize matrices
exp.sum.normed <- exp.sum / apply(exp.sum, 1, sum)
apply(exp.sum.normed, 1, range)
apply(exp.sum.normed, 2, range)
acc.sum.normed <- acc.sum / apply(acc.sum, 1, sum)
qs::qsave(exp.sum.normed, paste0(R.dir, "eGRN_sum_normed.qsave"))
qs::qsave(acc.sum.normed, paste0(R.dir, "eGRN_sum_normed.qsave"))


# Perform cell clustering based on expression summed matrix
# ARI of exp.sum: 0.07245299
exp.sum.d <- dist(t(exp.sum), method = "euclidean")
length(exp.sum.d)
class(exp.sum.d)
exp.sum.fit <- hclust(exp.sum.d, method = "ward")
plot(exp.sum.fit)
exp.sum.clusters <- cutree(exp.sum.fit, k = ncol(eGRN.ct))
identical(names(exp.sum.clusters), names(obj.celltyped$cell.type))
head(exp.sum.clusters)
head(obj.celltyped$cell.type)
exp.sum.ari <- ARI(exp.sum.clusters, obj.celltyped$cell.type)
exp.sum.ari


# Normalized expression sum matrix
# ARI of exp.sum.normed euclidean: 0.0726653
# ARI of exp.sum.normed maximum: 0.09276864 (bingo)
# ARI of exp.sum.normed manhattan: 0.07544453
# ARI of exp.sum.normed canberra: 0.0763576
# ARI of exp.sum.normed binary: 0.000368554
# ARI of exp.sum.normed minkowski: 0.0726653
exp.sum.normed[1:3, 1:3]
apply(exp.sum.normed, 1, sum)
exp.sum.normed.d <- dist(t(exp.sum.normed), method = "maximum")
length(exp.sum.normed.d)
class(exp.sum.normed.d)
exp.sum.normed.fit <- hclust(exp.sum.normed.d, method = "ward.D")
plot(exp.sum.normed.fit)
exp.sum.normed.clusters <- cutree(exp.sum.normed.fit, k = ncol(eGRN.ct))
exp.sum.normed.C6 <- cutree(exp.sum.normed.fit, k = 6)
identical(names(exp.sum.normed.clusters), names(obj.celltyped$cell.type))
identical(names(exp.sum.normed.C6), names(obj.celltyped$cell.type))
head(exp.sum.normed.clusters)
head(obj.celltyped$cell.type)
exp.sum.normed.ari <- ARI(exp.sum.normed.clusters, obj.celltyped$cell.type)
exp.sum.normed.ari
# ARI(exp.sum.normed.C6, obj.celltyped$cell.type)


# Accessibility matrix
# ARI: 0.0816956
acc.sum.normed[1:3, 1:3]
apply(acc.sum.normed, 1, sum)
acc.sum.normed.d <- dist(t(acc.sum.normed), method = "maximum")
length(acc.sum.normed.d)
class(acc.sum.normed.d)
acc.sum.normed.fit <- hclust(acc.sum.normed.d, method = "ward.D")
acc.sum.normed.clusters <- cutree(acc.sum.normed.fit, k = ncol(eGRN.ct))
identical(names(acc.sum.normed.clusters), names(obj.celltyped$cell.type))
head(acc.sum.normed.clusters)
head(obj.celltyped$cell.type)
acc.sum.normed.ari <- ARI(acc.sum.normed.clusters, obj.celltyped$cell.type)
acc.sum.normed.ari


# Try different clustering methods using maximum
# ward.D: 0.09276864 (bingo)
# ward.D2: 0.06859666
# single: 8.245143e-05
# complete: 0.001398377
# average: 0.001989182
# mcquitty: 0.001493673
# median: -2.509646e-05
# centroid: 0.0005179491
exp.sum.normed.d <- dist(t(exp.sum.normed), method = "maximum")
length(exp.sum.normed.d)
class(exp.sum.normed.d)


exp.sum.normed.fit <- hclust(exp.sum.normed.d, method = "ward.D")
exp.sum.normed.clusters <- cutree(exp.sum.normed.fit, k = ncol(eGRN.ct))
identical(names(exp.sum.normed.clusters), names(obj.celltyped$cell.type))
head(exp.sum.normed.clusters)
head(obj.celltyped$cell.type)
exp.sum.normed.ari <- ARI(exp.sum.normed.clusters, obj.celltyped$cell.type)
exp.sum.normed.ari


# Build eGRN occurrence matrix
occur.m <- Reduce("rbind", pbmclapply(eGRNs, function(x) {
  as.numeric(colnames(obj.celltyped) %in% x$cells)
}, mc.cores = detectCores()))
rownames(occur.m) <- seq_along(eGRNs)
colnames(occur.m) <- colnames(obj.celltyped)
qs::qsave(occur.m, paste0(R.dir, "eGRN_occur_matrix.qsave"))
dim(occur.m)


# Perform clustering and calculate ARI
# ARI: 0.008237966
exp.sum.occur <- exp.sum.normed * occur.m
exp.sum.occur.d <- dist(t(exp.sum.occur), method = "maximum")
length(exp.sum.occur.d)
class(exp.sum.occur.d)
exp.sum.occur.fit <- hclust(exp.sum.occur.d, method = "ward.D")
exp.sum.occur.clusters <- cutree(exp.sum.occur.fit, k = ncol(eGRN.ct))
identical(names(exp.sum.occur.clusters), names(obj.celltyped$cell.type))
head(exp.sum.occur.clusters)
head(obj.celltyped$cell.type)
exp.sum.occur.ari <- ARI(exp.sum.occur.clusters, obj.celltyped$cell.type)
exp.sum.occur.ari


# Perform PCA and then ward.D
dim(exp.sum.normed)
exp.sum.pca <- prcomp(x = t(exp.sum.normed), scale = F)
class(exp.sum.pca)
str(exp.sum.pca)
class(exp.sum.pca$x)
dim(exp.sum.pca$x)
colnames(exp.sum.pca$x)
rownames(exp.sum.pca$x)
exp.sum.pca$x[1:3, 1:3]


# Clustering after PCA
# ARI: 0.04950886
exp.sum.pca.d <- dist(exp.sum.pca$x, method = "maximum")
length(exp.sum.pca.d)
class(exp.sum.pca.d)
exp.sum.pca.fit <- hclust(exp.sum.pca.d, method = "ward.D")
exp.sum.pca.clusters <- cutree(exp.sum.pca.fit, k = ncol(eGRN.ct))
identical(names(exp.sum.pca.clusters), names(obj.celltyped$cell.type))
head(exp.sum.pca.clusters)
head(obj.celltyped$cell.type)
exp.sum.pca.ari <- ARI(exp.sum.pca.clusters, obj.celltyped$cell.type)
exp.sum.pca.ari


# Perform PCA on occurrence matrix and then ward.D
dim(exp.sum.occur)
exp.occur.pca <- prcomp(x = t(exp.sum.occur), scale = T)
class(exp.occur.pca)
str(exp.occur.pca)
class(exp.occur.pca$x)
dim(exp.occur.pca$x)
colnames(exp.occur.pca$x)
rownames(exp.occur.pca$x)
exp.occur.pca$x[1:3, 1:3]


# Clustering after PCA and occurrence
# PCA did not perform scaling
# maximum: 0.008520937
# euclidean: 0.00851466

# PCA performed scaling
# ARI: 0.009480027
# PC = 2: 0.009971001
exp.occur.pca.d <- dist(exp.occur.pca$x, method = "maximum")
length(exp.occur.pca.d)
class(exp.occur.pca.d)
exp.occur.pca.fit <- hclust(exp.occur.pca.d, method = "ward.D")
exp.occur.pca.clusters <- cutree(exp.occur.pca.fit, k = ncol(eGRN.ct))
identical(names(exp.occur.pca.clusters), names(obj.celltyped$cell.type))
head(exp.occur.pca.clusters)
head(obj.celltyped$cell.type)
exp.occur.pca.ari <- ARI(exp.occur.pca.clusters, obj.celltyped$cell.type)
exp.occur.pca.ari


# PCA and kmeans
# Default: 0.001791755
# PC = 10: 0.001430359
dim(exp.occur.pca$x)
exp.occur.pca.norm.kmeans <- kmeans(exp.occur.pca$x[, 1:10], centers = ncol(eGRN.ct))
exp.occur.pca.norm.kmeans$cluster
unique(exp.occur.pca.norm.kmeans$cluster)
identical(names(exp.occur.pca.norm.kmeans$cluster), names(obj.celltyped$cell.type))
head(exp.occur.pca.norm.kmeans$cluster)
head(obj.celltyped$cell.type)
exp.occur.pca.kmeans.ari <- ARI(exp.occur.pca.norm.kmeans$cluster, obj.celltyped$cell.type)
exp.occur.pca.kmeans.ari


# PCA is not necessary
# Use normalized and maximum
# ARI: 0.05870446
dim(exp.sum.normed)
exp.sum.norm.kmeans <- kmeans(t(exp.sum.normed), centers = ncol(eGRN.ct))
exp.sum.norm.kmeans$cluster
unique(exp.sum.norm.kmeans$cluster)
identical(names(exp.sum.norm.kmeans$cluster), names(obj.celltyped$cell.type))
head(exp.sum.norm.kmeans$cluster)
head(obj.celltyped$cell.type)
exp.sum.kmeans.ari <- ARI(exp.sum.norm.kmeans$cluster, obj.celltyped$cell.type)
exp.sum.kmeans.ari


# Kmeans without normalization
# ARI: 0.05202223
exp.sum
exp.sum.kmeans <- kmeans(t(exp.sum), centers = ncol(eGRN.ct))
exp.sum.kmeans$cluster
unique(exp.sum.kmeans$cluster)
identical(names(exp.sum.kmeans$cluster), names(obj.celltyped$cell.type))
head(exp.sum.kmeans$cluster)
head(obj.celltyped$cell.type)
exp.sum.noNorm.kmeans.ari <- ARI(exp.sum.kmeans$cluster, obj.celltyped$cell.type)
exp.sum.noNorm.kmeans.ari


# Kmeans without normalization; occurrence
# ARI: 0.0009066792
exp.sum
exp.sum.occur2 <- exp.sum * occur.m
exp.sum.occur2.kmeans <- kmeans(t(exp.sum.occur2), centers = ncol(eGRN.ct))
exp.sum.occur2.kmeans$cluster
unique(exp.sum.occur2.kmeans$cluster)
identical(names(exp.sum.occur2.kmeans$cluster), names(obj.celltyped$cell.type))
head(exp.sum.occur2.kmeans$cluster)
head(obj.celltyped$cell.type)
exp.sum.occur2.noNorm.kmeans.ari <- ARI(exp.sum.occur2.kmeans$cluster, obj.celltyped$cell.type)
exp.sum.occur2.noNorm.kmeans.ari


# Inspect hierarchical clustering after normalization, occurrence, maximum, and ward.D
exp.sum.normed.clusters
obj.celltyped$cell.type
clust.ct.m <- Reduce("rbind", pbmclapply(unique(exp.sum.normed.clusters), function(x) {
  overlap <- intersect(names(obj.celltyped$cell.type), names(which(exp.sum.normed.clusters == x)))
  table(obj.celltyped$cell.type[overlap])
}, mc.cores = detectCores()))
dim(clust.ct.m)
rownames(clust.ct.m) <- unique(exp.sum.normed.clusters)
clust.ct.m[1:3, 1:3]
clust.ct.m
clust.ct.ratio <- clust.ct.m / apply(clust.ct.m, 1, sum)
clust.ct.ratio
apply(clust.ct.ratio, 1, max)
# Clusters 2, 3, 5, and 10 are confused!!!


# Focus on occurrence matrix normalized by rows
# ARI: 0.002402324
dim(occur.m)
apply(occur.m, 1, sum)
occur.row.normed <- occur.m / apply(occur.m, 1, sum)
occur.row.normed[1:3, 1:3]
occur.row.normed.kmeans <- kmeans(t(occur.row.normed), centers = ncol(eGRN.ct))
occur.row.normed.kmeans$cluster
unique(occur.row.normed.kmeans$cluster)
identical(names(occur.row.normed.kmeans$cluster), names(obj.celltyped$cell.type))
head(occur.row.normed.kmeans$cluster)
head(obj.celltyped$cell.type)
occur.row.normed.kmeans.ari <- ARI(occur.row.normed.kmeans$cluster, obj.celltyped$cell.type)
occur.row.normed.kmeans.ari


# maximum: 0.02769357
# euclidean: 0.03587221
dim(exp.sum)
exp.col.normed <- exp.sum / apply(exp.sum, 2, sum)
exp.col.d <- dist(t(exp.col.normed), method = "euclidean")
length(exp.col.d)
class(exp.col.d)
exp.col.fit <- hclust(exp.col.d, method = "ward.D")
exp.col.clusters <- cutree(exp.col.fit, k = ncol(eGRN.ct))
identical(names(exp.col.clusters), names(obj.celltyped$cell.type))
head(exp.col.clusters)
head(obj.celltyped$cell.type)
exp.col.ari <- ARI(exp.col.clusters, obj.celltyped$cell.type)
exp.col.ari


# Row normalized
# euclidean: 0.06332401
# maximum: 0.06323437
# multiply 100: 0.06280073
# 1000: 0.00791905
# 10: 0.0883072
dim(exp.sum)
exp.max.normed <- exp.sum / apply(exp.sum, 1, max)
exp.max.normed <- exp.max.normed + exp.max.normed * occur.m * 10
exp.max.normed <- exp.max.normed / apply(exp.max.normed, 1, sum)
exp.max.d <- dist(t(exp.max.normed), method = "maximum")
length(exp.max.d)
class(exp.max.d)
exp.max.fit <- hclust(exp.max.d, method = "ward.D")
exp.max.clusters <- cutree(exp.max.fit, k = nrow(eGRN.ct))
identical(names(exp.max.clusters), names(obj.celltyped$cell.type))
head(exp.max.clusters)
head(obj.celltyped$cell.type)
exp.max.ari <- ARI(exp.max.clusters, obj.celltyped$cell.type)
exp.max.ari
