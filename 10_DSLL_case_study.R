###########################################################
#                                                         #
#                Perform case study on DSLL               #
#                                                         #
###########################################################


# Analysis plan:
# 1. Filter out the CLL and lymphoma related TFs and corresponding eGRNs
# 2. Screen the eGRNs by retaining the ones significantly enriched in cell types, 
#    i.e., cell-type-specific eGRNs
# 3. Run CellChat to discover cell-cell communications (CCCs) among cell types
# 4. Retain TFs showing difference between normal and diseased cell types
# 5. Pathway enrichment analyses of these cell-type-specific eGRNs and extract the
#    pathways related with CLL
# 6. Find the ligand-receptor pairs corresponding to eGRN genes
# 7. Follow the sequential order defined by CCCs
# 8. The three phenotypes should be analyzed one by one
# If possible, we had better exploit this case in a sequential order


# Relation among IgHV, TP53, and survival

# Weblink: https://news.mayocliniclabs.com/2019/08/12/ighv-and-tp53-sequencing-clinical-utility-in-chronic-lymphocytic-leukemia-cll/

# IgHV:
# The IgVH gene mutation status is one of the discriminators of clinical outcome in patients with CLL. The mutational status of the immunoglobulin genes expressed by CLL cells can be used to segregate patients into two subsets that have significantly different tendencies for disease progression.
# IgHV is better for patients.

# Somatic hypermutation is a process that allows B cells to mutate the genes that they use to produce antibodies. This enables the B cells to produce antibodies that are better able to bind to bacteria, viruses and other infections.
# Somatic hyper-mutation of IGHV is a stable marker in CLL, with unmutated IGHV predicting for an unfavorable response to chemoimmunotherapy.
# Mutational status doesn’t seem to make a difference and both groups of patients continue to benefit equally as they stay on therapy.
# Therefore, IgHV may not be predictive for better patient survival.

# TP53:
# The TP53 gene provides instructions for making a protein called tumor protein p53 (or p53). This protein acts as a tumor suppressor, which means that it regulates cell division by keeping cells from growing and dividing (proliferating) too fast or in an uncontrolled way.


# Analysis pipeline:
# 1. TP53 mutation leads to worse patient survival.
# 2. Calculate enrichment of TP53-enriched and worse survival-enriched eGRNs against 
#    normal B cells and DSLL B states, respectively. We obtain cell-type-specific eGRNs.
# 3. Run pathway enrichment for eGRNs specific to both phenotype and cell type.
# 4. Run CellChat to identify CCCs between various B cell states.
# 5. Identify eGRNs including ligands in normal or tumor B cells that are enriched 
#    in pathways incorporating TP53.
# 6. Discover eGRNs including receptors in tumor B cells that are enriched in pathways
#    that may lead to worse survival, i.e., BAFF signaling pathways.
# 7. IgHV mutation results in better patient survival, which is under debate.
# 8. Follow steps 2-6 between IgHV mutated cells v.s. better survival-associated cells.




# Libraries
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)




# Global parameters
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Rfile/"
load(paste0(R.dir, "Scissor_tuned/Scissor_tuned.RData"))
table.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Tables/TP53_to_survival/"
setwd(R.dir)
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
image.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Images/TP53_to_survival/"
dir.create(image.dir)



# Source codes
source(paste0(tool.dir, "transcriptome_tools.R"))




# 2. Calculate enrichment of TP53-enriched and worse survival-enriched eGRNs against 
#    normal B cells and DSLL B states, respectively. We obtain cell-type-specific eGRNs.


# TP53-specific eGRN prediction
t.TP53.plus.ids # none
length(other.TP53.plus.ids) # 32 eGRNs significantly enriched in other cell types
other.TP53.plus.pval # plot: a boxplot
levels(lymph.obj$cell.type)
t.celltypes
other.celltypes
padj.df <- as.data.frame(padj.dt)
colnames(padj.df)
dim(padj.df[, as.character(other.celltypes)])
other.padj <- padj.df[, as.character(other.celltypes)]
dim(other.padj)
other.cts.ids <- which(apply(other.padj, 1, min) < 0.05)
length(other.cts.ids) # 40
tp53.cts.ids <- intersect(other.TP53.plus.ids, other.cts.ids)
tp53.cts.ids
length(tp53.cts.ids) # 32 cell-type-specific eGRNs enriched in TP53 mutated cells
sapply(enhancer.grns[tp53.cts.ids], "[[", "TF") %>% unique %>% length # 24 TFs
enhancer.grn.genes <- lapply(enhancer.grns, "[[", "genes")
"TP53" %in% Reduce("union", enhancer.grn.genes[tp53.cts.ids]) # TP53 is not regulated
table(sapply(enhancer.grns[tp53.cts.ids], "[[", "TF"))


# # Pathway enrichment analysis (may bee unnecessary)
# pathways <- run_GO_and_KEGG(genes.ll = enhancer.grn.genes, dbs = "KEGG")$KEGG_2019_Human
# dim(pathways) # 7058
# length(unique(pathways$Id)) # 50
# qs::qsave(R.dir, "Enriched_KEGG_pathways.qsave")
# pathways <- pathways[pathways$Adjusted.P.value < 0.05,] # 1318
# which(grepl("ignaling", pathways$Term)) %>% length # 338 signaling pathways
# which(grepl("ctivat", pathways$Term)) %>% length # 21 activation pathways
# which(grepl("B cell", pathways$Term)) %>% length # 33 pathways related with B cells


tp53.enriched.celltypes <- apply(padj.df[tp53.cts.ids,], 1, function(x) {
  names(which(x < 0.05))
})
tp53.enriched.celltypes
sapply(enhancer.grns[tp53.cts.ids], "[[", "TF")
sapply(enhancer.grns[tp53.cts.ids], "[[", "TF") %>% table



# cox-specific eGRN prediction
t.cox.plus.ids # none
length(other.cox.plus.ids) # 17 eGRNs significantly enriched in other cell types
other.cox.plus.pval # plot: a boxplot
length(other.cts.ids) # 40
cox.cts.ids <- intersect(other.cox.plus.ids, other.cts.ids)
cox.cts.ids
length(cox.cts.ids) # 17 cell-type-specific eGRNs enriched in cox mutated cells
sapply(enhancer.grns[cox.cts.ids], "[[", "TF") %>% unique %>% length # 14 TFs
table(sapply(enhancer.grns[cox.cts.ids], "[[", "TF"))


# # Pathway enrichment analysis (may bee unnecessary)
# pathways <- run_GO_and_KEGG(genes.ll = enhancer.grn.genes, dbs = "KEGG")$KEGG_2019_Human
# dim(pathways) # 7058
# length(unique(pathways$Id)) # 50
# qs::qsave(R.dir, "Enriched_KEGG_pathways.qsave")
# pathways <- pathways[pathways$Adjusted.P.value < 0.05,] # 1318
# which(grepl("ignaling", pathways$Term)) %>% length # 338 signaling pathways
# which(grepl("ctivat", pathways$Term)) %>% length # 21 activation pathways
# which(grepl("B cell", pathways$Term)) %>% length # 33 pathways related with B cells


cox.enriched.celltypes <- apply(padj.df[cox.cts.ids,], 1, function(x) {
  names(which(x < 0.05))
})
cox.enriched.celltypes
sapply(enhancer.grns[cox.cts.ids], "[[", "TF")
sapply(enhancer.grns[cox.cts.ids], "[[", "TF") %>% table


# ighv-specific eGRN prediction
t.ighv.plus.ids # none
t.ighv.minus.ids # none
length(other.ighv.plus.ids) # 32 eGRNs significantly enriched in other cell types
other.ighv.plus.ids
other.ighv.plus.pval # plot: a boxplot
other.ighv.minus.ids # none
length(other.cts.ids) # 40
ighv.cts.ids <- intersect(other.ighv.plus.ids, other.cts.ids)
ighv.cts.ids
length(ighv.cts.ids) # 32 cell-type-specific eGRNs enriched in ighv mutated cells
sapply(enhancer.grns[ighv.cts.ids], "[[", "TF") %>% unique %>% length # 25 TFs
table(sapply(enhancer.grns[ighv.cts.ids], "[[", "TF"))


# # Pathway enrichment analysis (may bee unnecessary)
# pathways <- run_GO_and_KEGG(genes.ll = enhancer.grn.genes, dbs = "KEGG")$KEGG_2019_Human
# dim(pathways) # 7058
# length(unique(pathways$Id)) # 50
# qs::qsave(R.dir, "Enriched_KEGG_pathways.qsave")
# pathways <- pathways[pathways$Adjusted.P.value < 0.05,] # 1318
# which(grepl("ignaling", pathways$Term)) %>% length # 338 signaling pathways
# which(grepl("ctivat", pathways$Term)) %>% length # 21 activation pathways
# which(grepl("B cell", pathways$Term)) %>% length # 33 pathways related with B cells


ighv.enriched.celltypes <- apply(padj.df[ighv.cts.ids,], 1, function(x) {
  names(which(x < 0.05))
})
ighv.enriched.celltypes
sapply(enhancer.grns[ighv.cts.ids], "[[", "TF")
sapply(enhancer.grns[ighv.cts.ids], "[[", "TF") %>% table


# Comparison between TP53 and cox
tp53.cox.TFs <- intersect(sapply(enhancer.grns[cox.cts.ids], "[[", "TF") %>% unique, 
                          sapply(enhancer.grns[tp53.cts.ids], "[[", "TF") %>% unique)
tp53.cox.TFs
length(tp53.cox.TFs) # 14




# 3. Run pathway enrichment for eGRNs specific to both phenotype and cell type.
# Herein, we use CellChat to identify the signaling pathways.
levels(lymph.obj$cell.type)
short.celltypes <- lymph.obj$cell.type
levels(short.celltypes) <- c("Tumor B", "Exh CD4+ T", 
                             "Exh CD8+ T", "CM CD8+ T", 
                             "Eff CD4+ T_1", "TAMs", 
                             "Eff CD4+ T_2", "Eff CD8+ T", 
                             "Normal B", "Tumor B prol", 
                             "DCs")
head(short.celltypes)
lymph.obj <- AddMetaData(lymph.obj, metadata = short.celltypes, 
                         col.name = "short.celltypes")
qs::qsave(lymph.obj, paste0(R.dir, "lymph_obj.qsave"))
cellchat.meta <- as.data.frame(lymph.obj$short.celltypes)
colnames(cellchat.meta) <- "labels"
dim(cellchat.meta)
head(cellchat.meta)
cellchat.exp <- normalizeData(GetAssayData(lymph.obj, slot = "data", assay = "RNA"))
dim(cellchat.exp)
cellchat.exp[1:5, 1:5]
cellchat.obj <- createCellChat(object = cellchat.exp, meta = cellchat.meta, 
                           group.by = "labels")
qs::qsave(cellchat.obj, paste0(R.dir, "CellChat_object.qsave"))


# Show the structure of the database
showDatabaseCategory(CellChatDB.human)
dplyr::glimpse(CellChatDB.human$interaction)


# Subset the expression data of signaling genes for saving computation cost
cellchat.obj@DB <- CellChatDB.human
cellchat.obj <- subsetData(cellchat.obj) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat.obj <- identifyOverExpressedGenes(cellchat.obj)
cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)


# Compute the communication probability and infer cellular communication network
cellchat.obj <- computeCommunProb(cellchat.obj)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.obj <- filterCommunication(cellchat.obj, min.cells = 10)


# Infer the cell-cell communication at a pathway level
cellchat.obj <- computeCommunProbPathway(cellchat.obj)


# Calculate the aggregated cell-cell communication network
cellchat.obj <- aggregateNet(cellchat.obj)
groupSize <- as.numeric(table(cellchat.obj@idents))
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat.obj@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.obj@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Interaction weights/strength")


# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
ccc.mat <- cellchat.obj@net$weight
tiff(filename = paste0(image.dir, "DSLL_cell_cell_commun.tiff"), 
     width = 1500, height = 2000, res = 150)
par(mfrow = c(3, 4), mar = c(1, 1, 1, 1), xpd = T)
for (i in 1 : nrow(ccc.mat)) {
  ccc.mat2 <- matrix(0, nrow = nrow(ccc.mat), ncol = ncol(ccc.mat), 
                 dimnames = dimnames(ccc.mat))
  ccc.mat2[i, ] <- ccc.mat[i, ]
  netVisual_circle(ccc.mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(ccc.mat), title.name = rownames(ccc.mat)[i])
}
dev.off()


# # Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
# pathways.show <- c("TP53") 
# vertex.receiver = seq(1, 4) # a numeric vector. 
# netVisual_aggregate(cellchat.obj, signaling = pathways.show,  
#                     vertex.receiver = vertex.receiver)
# # Circle plot
# par(mfrow = c(1, 1))
# netVisual_aggregate(cellchat.obj, signaling = pathways.show, layout = "circle")


# Check the CCC pathways connecting TP53 with survival
lig.rec.pairs <- subsetCommunication(cellchat.obj)
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
dim(lig.rec.pairs)
head(lig.rec.pairs)
colnames(lig.rec.pairs)
lig.genes <- unique(lig.rec.pairs$ligand)
rec.genes <- unique(lig.rec.pairs$receptor)
length(lig.genes) # 60
length(rec.genes) # 51
t.TP53.plus.ids # none
length(tp53.cts.ids) # 32
sapply(enhancer.grns[tp53.cts.ids], "[[", "TF") %>% table


# Co-factors of TP53 (10): SP1, EGR1, JUN, TFAP2C, FOS, KLF5, MYOG, NR4A1, TFAP2C, YY2
tp53.cof <- c("SP1", "EGR1", "JUN", "TFAP2C", "FOS", 
              "KLF5", "MYOG", "NR4A1", "TFAP2C", "YY2")
tp53.cts.cof.ids <- tp53.cts.ids[sapply(enhancer.grns[tp53.cts.ids], "[[", "TF") %in% tp53.cof]
tp53.cts.cof.ids # 10
tp53.cof.celltypes <- Reduce("union", tp53.enriched.celltypes[as.character(tp53.cts.cof.ids)])
tp53.cof.celltypes # 9


# Find target genes overlapped with ligands regulated by co-factors of TP53
# and receptors leading to worse survival
tp53.cts.cof.targets <- Reduce("union", lapply(enhancer.grns[tp53.cts.cof.ids], 
                                               "[[", "genes"))
length(tp53.cts.cof.targets) # 408
tp53.cts.cof.ligands <- intersect(tp53.cts.cof.targets, lig.genes)
length(tp53.cts.cof.ligands) # 10
tp53.cts.cof.ligands
# tp53.cts.cof.receptors <- intersect(tp53.cts.cof.targets, rec.genes)
# length(tp53.cts.cof.receptors) # 6
# tp53.cts.cof.receptors
cox.cts.targets <- Reduce("union", lapply(enhancer.grns[cox.cts.ids], 
                                               "[[", "genes"))
length(cox.cts.targets) # 424
cox.cts.receptors <- intersect(cox.cts.targets, rec.genes)
length(cox.cts.receptors) # 4


# Locate the ligand-receptor interactions and asssociated cell types
lig.rec.pairs[lig.rec.pairs$ligand %in% tp53.cts.cof.ligands & 
                lig.rec.pairs$receptor %in% cox.cts.receptors, c(1:4, 7, 9)]
# To-do: Tumor B prol and Tumor B cells
