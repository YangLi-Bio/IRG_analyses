#######################################################
#                                                     #
#                  Obtain a good ARI                  #
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


# Activity scores in references
# SCENIC: AUC
# IRIS3: Regulon activity
# Cicero: gene activity + PCA
# cisTopic: topic probability
# sc-compReg & coupled NMF: decomposed matrix


# Try calculating Jaccard index between cells
range(occur.m)
dim(occur.m)
eGRN.cells <- Reduce("union", lapply(eGRNs, "[[", "cells"))
length(eGRN.cells)
obj.celltyped$cell.type
subset.cells <- obj.celltyped$cell.type[eGRN.cells]
length(subset.cells)
subset.m <- occur.m[, names(subset.cells), drop = F]
dim(subset.m)
union.df <- do.call("rbind", pbmclapply(2 : ncol(subset.m), function(i) {
  i.set <- which(subset.m[, i] > 0)
  Reduce("rbind", lapply(1:(i - 1), function(j) {
    j.set <- which(subset.m[, j] > 0)
    ij.jacc <- length(intersect(i.set, j.set)) / 
      length(union(i.set, j.set))
    c(i, j, ij.jacc)
  }))
}, mc.cores = detectCores()))
max.df <- do.call("rbind", pbmclapply(2 : ncol(subset.m), function(i) {
  i.set <- which(subset.m[, i] > 0)
  Reduce("rbind", lapply(1:(i - 1), function(j) {
    j.set <- which(subset.m[, j] > 0)
    ij.jacc <- length(intersect(i.set, j.set)) / 
      max(length(i.set), length(j.set))
    c(i, j, ij.jacc)
  }))
}, mc.cores = detectCores()))
min.df <- do.call("rbind", pbmclapply(2 : ncol(subset.m), function(i) {
  i.set <- which(subset.m[, i] > 0)
  Reduce("rbind", lapply(1:(i - 1), function(j) {
    j.set <- which(subset.m[, j] > 0)
    ij.jacc <- length(intersect(i.set, j.set)) / 
      min(length(i.set), length(j.set))
    c(i, j, ij.jacc)
  }))
}, mc.cores = detectCores()))
dim(jacc.df)
dim(max.df)
dim(min.df)
head(jacc.df)
head(max.df)
head(min.df)


# union: 0.2491488
# max: 0.08249952
# min: 0.329393
jacc.df <- min.df


jacc.total <- rbind(jacc.df, Reduce("rbind", pbmclapply(1:ncol(subset.m), function(i) {
  c(i, i, 1)
}, mc.cores = detectCores())))
colnames(jacc.total) <- c("i", "j", "x")
# jacc.m <- sparseMatrix(i = jacc.total[, 1], j = jacc.total[, 2], x = jacc.total[, 3])
# dim(jacc.m)
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
# regularMatrix <- regularMatrix / apply(regularMatrix, 1, max)
regularMatrix[1:3, 1:3]


# Try different values
# 1: 0.2491488 (Bingo)
# 2: 0.07530222
# 4: 0.06349345
# 8: 0.09066146
# 16: 0.126822
# 32: 0.06485978
disRegularMatrix <- 1 - regularMatrix
regularMatrix[1:3, 1:3]
disRegularMatrix[1:3, 1:3]
disRegularMatrix <- disRegularMatrix / apply(disRegularMatrix, 1, sum)


# Clustering
# ward.D2: 0.2491488 (Bingo)
# ward.D: 0.0785303
# single: 0.05874954
# complete: 0.09687345
# average: 0.006195763
# mcquitty: 0.006901614
# median: 0.003949436
# centroid: 0.00707892
distanceMatrix <- as.dist(disRegularMatrix)
length(distanceMatrix)
class(distanceMatrix)
distanceMatrix[1:3]
distanceMatrix[is.na(distanceMatrix)] <- 0
clusters <- hclust(distanceMatrix, method = "ward.D2")
plot(clusters)
length(unique(subset.cells))


# Try different k
# 22: 0.2014389
# 11: 0.329393 (Bingo)
# 6: 0.330807
# 3: 0.09064958
group <- cutree(clusters, k = 11)
group
head(group)
head(subset.cells)
jacc.ari <- ARI(group, subset.cells)
jacc.ari


# Try graph clustering
# Only one cluster was returned!!!
regularMatrix[1:3, 1:3]

# g <- graph_from_data_frame(jacc.total, directed = F)
# V(g)
# E(g)
# E(g)$weight <- jacc.total[, 3]

results <- mst.knn(disRegularMatrix, suggested.k = 11)
results
str(results)
names(results)
unique(results$cluster)
plot(results$network)


# Graph clustering
# Each single vertex is a cluster!!!
# unweighted: 0.143473
# weighted: 0.1434438


# 0.10: 0.1434744
# 0.20: 0.1443595
# 0.30: 0.1487188
# 0.40: 0.3147129
# 0.50: 0.318827
# 0.75: 0.3286345
# 0.90: 0.3302682
# 0.95: 0.3199338


# 0.95:
# resolution:
# 0.10: 0.3302682
# 0.05: 0.1594658


# n_iterations
# 2: 0.3302682
# 10: 0.3292472
# 20: 0.3298457


# beta
# 0.01: 0.3302682
# 0.1: 0.3241771


# objective_function
# modularity: 
# CPM: 0.3257425


# Cluster_louvain
# Default: 0.2556639
# weights: 0.2556639


# edge
# 0.95: 0.2556639
# 0.90: 0.2623397
# 0.50: 0.2445445
# 0.10: 0.2923139
jacc.df <- as.data.frame(jacc.df)
colnames(jacc.df) <- c("i", "j", "x")
range(jacc.df$x) 
dim(jacc.df)
boxplot(jacc.df$x)
g.df <- jacc.df[jacc.df$x > 0.90,]
g.df <- jacc.df
dim(g.df)
g <- graph_from_data_frame(g.df, directed = F)
E(g)$weight <- g.df$x


# optimal: 0.2623397
# 
# g.clusters <- cluster_leiden(g, resolution_parameter = 0.10)
# g.clusters <- cluster_louvain(g, weights = E(g)$weight)
# g.clusters <- cluster_optimal(g, weights = E(g)$weight)
g.clusters <- cluster_spinglass(g, weights = E(g)$weight)
unique(g.clusters$membership)
g.clusters$membership
length(g.clusters$membership)
cluster.v <- g.clusters$membership
names(cluster.v) <- colnames(subset.m)
head(cluster.v)
ARI(cluster.v, subset.cells)


# Conclusion: try other clustering algorithms in addition to graph clustering
# If no better result exists, use hierarchical clustering!!!
# Showcase the clustering results; remove several general eGRNs, and perform 
# clustering again and check whether better results can be obtained.