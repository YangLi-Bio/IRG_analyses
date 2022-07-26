#######################################################
#                                                     #
#         Extend the set of CREs in an eGRN           #
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


# Source codes
source(paste0(code.dir, "reguome_tools.R"))
source(paste0(tool.dir, "cistrome_tools.R"))


# Load the Seurat object
obj <- qs::qread(paste0(R.dir, "lymph_obj.qsave"))
names(obj@assays)
m.rna <- GetAssayData(object = obj, slot = "data", assay = 'RNA')
m.atac <- GetAssayData(object = obj, slot = "data", assay = 'ATAC')
dim(m.rna)
dim(m.atac)


# Load the IRG-eGRNs
IRG.eGRNs <- qs::qread(paste0(R.dir, "irg_egrns.qsave"))
length(IRG.eGRNs)
sapply(IRG.eGRNs, "[[", "genes") %>% sapply(., length) %>% range
sapply(IRG.eGRNs, "[[", "peaks") %>% sapply(., length) %>% range
sapply(IRG.eGRNs, "[[", "cells") %>% sapply(., length) %>% range


# Annotate TFBSs
load("/fs/ess/PCON0022/liyang/STREAM/Codes/stream/data/TFBS_list.rda")
TF.CRE.pairs <- find_TFBS(m.atac, TFBS.list = TFBS.list, org = "hg38")
length(TF.CRE.pairs)
names(TF.CRE.pairs)
bound.TFs <- TF.CRE.pairs$CRE
binding.CREs <- TF.CRE.pairs$TF
rm(TF.CRE.pairs)
length(bound.TFs)
length(binding.CREs)


# Add putative CREs to each IRG-eGRN
extended.eGRNs <- pblapply(IRG.eGRNs, function(x) {
  length(x$genes)
  length(x$peaks)
  length(x$cells)
  coherence.cutoff <- x$atac.ratio * coherence.mpx * length(x$cells)
  peak.pool <- binding.CREs[[x$TF]]
  length(peak.pool)
  m.pool <- (m.atac > 0)[peak.pool, x$cells, drop = F]
  dim(m.pool)
  m.pool[1:10, 1:10]
  colnames(m.pool)
  add.peaks <- which(apply(m.pool, 1, sum) >= coherence.cutoff) %>% names
  gene.annotations <- obj[['ATAC']]@annotation %>% CollapseToLongestTranscript
  gene.annotations <- gene.annotations[gene.annotations$gene_name %in% x$genes]
  head(gene.annotations)
  peak.GR <- StringToGRanges(add.peaks)
  head(peak.GR)
  summ <- DistanceToTSS(peaks = peak.GR, genes = gene.annotations, 
                        distance = distance, sep = c("-", "-")) %>% summary
  summ
  cis.gr <- peak.GR[summ$i]
  mcols(cis.gr)$gene <- gene.annotations$gene_name[summ$j]
  cis.gr
  x[['peaks']] <- union(x[['peaks']], unique(GRangesToString(cis.gr)))
  x[["links"]] <- cis.gr
  x
})
qs::qsave(extended.eGRNs, paste0(R.dir, "Extended_eGRNs.qsave"))


# Comparison between IRG-eGRNs and extended eGRNs
sapply(IRG.eGRNs, "[[", "genes") %>% sapply(., length) %>% range
sapply(extended.eGRNs, "[[", "genes") %>% sapply(., length) %>% range

sapply(IRG.eGRNs, "[[", "peaks") %>% sapply(., length) %>% range
sapply(extended.eGRNs, "[[", "peaks") %>% sapply(., length) %>% range

sapply(IRG.eGRNs, "[[", "cells") %>% sapply(., length) %>% range
sapply(extended.eGRNs, "[[", "cells") %>% sapply(., length) %>% range
