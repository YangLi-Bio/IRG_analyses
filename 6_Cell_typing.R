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


# Global parameters
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Rfile/"
coherence.mpx <- 1.0 # multiple of ATAC ratio
distance <- 250000


# Source scripts
source(paste0(tool.dir, "multiome_tools.R"))


# Load the object
obj <- qs::qread(paste0(R.dir, "lymph_obj.qsave"))
dim(obj)
names(obj@assays)
dim(obj[['ATAC']])
dim(obj[['RNA']])
dim(obj[['SCT']])
eGRNs <- qs::qread(paste0(R.dir, "Extended_eGRNs.qsave"))


# Load the cell types from Tong Xiao
ct.vec <- readRDS(paste0(R.dir, "lymphoma_CT.RDS"))
length(ct.vec)
table(ct.vec)
length(table(ct.vec))
ncol(obj)
length(intersect(names(ct.vec), colnames(obj))) == length(ct.vec)
obj.clustered <- cluster_multiome(pbmc = obj)
qs::qsave(obj.clustered, paste0(R.dir, "obj_clustered.qsave"))


# Calculate the cell type maximally overlapped with annotated cell type
Idents(obj.clustered) <- obj.clustered$seurat_clusters
cluster.ct.dt <- rbindlist(pbmclapply(levels(obj.clustered$seurat_clusters), function(i) {
  cells.cluster <- colnames(obj.clustered)[which(Idents(obj.clustered) == i)]
  overlapped <- intersect(cells.cluster, names(ct.vec))
  as.list(table(ct.vec[overlapped]))
}, mc.cores = detectCores()))
rownames(cluster.ct.dt) <- levels(obj.clustered$seurat_clusters)
dim(cluster.ct.dt)
head(cluster.ct.dt)
cluster.ct.ratio <- cluster.ct.dt / apply(cluster.ct.dt, 1, sum)
rownames(cluster.ct.ratio) <- levels(obj.clustered$seurat_clusters)
dim(cluster.ct.ratio)
head(cluster.ct.ratio)
qs::qsave(cluster.ct.ratio, paste0(R.dir, "Cluster_cell_type_ratios.qsave"))


# Check the cell clusters
apply(cluster.ct.ratio, 1, max) %>% range
apply(cluster.ct.ratio, 1, max)
which(apply(cluster.ct.ratio, 1, max) <= 0.50) %>% length


# Perform cell subclustering for cluster IDs 2, 3, 5, 13 (not levels)
ids.subclustered <- levels(obj.clustered$seurat_clusters)[c(2, 3, 5, 13)]
ids.subclustered
obj.subclustered <- FindSubCluster(
  obj.clustered,
  graph.name = "wsnn",
  cluster = ids.subclustered, 
  resolution = 0.1
)
table(obj.subclustered$sub.cluster)


# Cell cubclusters
Idents(obj.subclustered) <- obj.subclustered$sub.cluster
subcluster.ct.dt <- rbindlist(pbmclapply(unique(obj.subclustered$sub.cluster), function(i) {
  cells.cluster <- colnames(obj.subclustered)[which(Idents(obj.subclustered) == i)]
  overlapped <- intersect(cells.cluster, names(ct.vec))
  as.list(table(ct.vec[overlapped]))
}, mc.cores = detectCores()))
rownames(subcluster.ct.dt) <- unique(obj.subclustered$sub.cluster)
rownames(subcluster.ct.dt)
dim(subcluster.ct.dt)
head(subcluster.ct.dt)
subcluster.ct.ratio <- subcluster.ct.dt / apply(subcluster.ct.dt, 1, sum)
rownames(subcluster.ct.ratio) <- rownames(subcluster.ct.dt)
dim(subcluster.ct.ratio)
head(subcluster.ct.ratio)
qs::qsave(subcluster.ct.ratio, paste0(R.dir, "Subcluster_cell_type_ratios.qsave"))


# Check the cell subclusters
apply(subcluster.ct.ratio, 1, max) %>% range
apply(subcluster.ct.ratio, 1, max)
rownames(subcluster.ct.ratio)
which(apply(subcluster.ct.ratio, 1, max) <= 0.50) %>% length


# Annotate cell types according to cell clusters
cluster.ct.ratio
apply(cluster.ct.ratio, 1, which.max) %>% unique %>% length
ct.list <- names(table(ct.vec))
ct.list
levels(obj.subclustered$seurat_clusters)
cluster.ct.map <- colnames(cluster.ct.ratio)[apply(cluster.ct.ratio, 1, which.max)]
names(cluster.ct.map) <- rownames(cluster.ct.ratio)
head(cluster.ct.map)
obj.celltyped <- obj.subclustered
obj.celltyped[['cell.type']] <- obj.celltyped$seurat_clusters
obj.celltyped$cell.type
levels(obj.celltyped$cell.type) <- cluster.ct.map
levels(obj.celltyped$cell.type)
qs::qsave(obj.celltyped, paste0(R.dir, "Obj_celltyped.qsave"))


# Compare the eGRN-active cell subpopulation with cell types
sapply(eGRNs, "[[", "cells") %>% sapply(., length) %>% range
eGRN.ct <- rbindlist(pbmclapply(eGRNs, function(x) {
  x.ct <- intersect(x$cells, names(obj.celltyped$cell.type))
  table(obj.celltyped$cell.type[x.ct]) %>% as.list
}, mc.cores = detectCores()))
qs::qsave(eGRN.ct, paste0(R.dir, "eGRN_active_cells_in_CT.qsave"))
