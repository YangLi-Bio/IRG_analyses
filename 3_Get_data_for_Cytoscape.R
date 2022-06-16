##################################################################################
#                                                                                #
#                            Prepare data for Cytoscape                          #
#                                                                                #
##################################################################################


# Parameters
work.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/"
IgHV.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/1_classification_GSE113386/"
TP53.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/3_TP53_GSE11038/"
cox.dir <- "/fs/ess/PCON0022/liyang/STREAM/case_1_IRG/2_survival_GSE22762/"


# Load eGRNs
IgHV.eGRNs <- qs::qread(paste0(IgHV.dir, "DEG_eGRNs.qsave"))
TP53.eGRNs <- qs::qread(paste0(TP53.dir, "DEG_eGRNs.qsave"))
cox.eGRNs <- qs::qread(paste0(cox.dir, "DEG_eGRNs.qsave"))


# Important eGRN IDs
IgHV.IDs <- c(12, 42, 49, 44) # MYCN (x 2), JUN, STAT3
TP53.IDs <- c(24, 45, 49, 44) # MYOG (x 2), JUN, STAT3
cox.IDs <- c(12, 42, 49, 50) # MYCN (x 2), JUN, KLF5


# # Selected eGRNs for IgHV, TP53, and cox
# IgHV.selected <- IgHV.eGRNs[IgHV.IDs]
# TP53.selected <- TP53.eGRNs[TP53.IDs]
# cox.selected <- cox.eGRNs[cox.IDs]


# The eGRNs composed of all IRGs
IRG.eGRNs <- qs::qread(paste0(work.dir, "irg_egrns.qsave"))
length(IRG.eGRNs)


##################################################################################
#                                                                                #
#            Determine the TFs for Cytoscape to generate networks                #
#                                                                                #
##################################################################################


# If necessary, only retain the CREs overlapped with enhancers
# Perform functional genomics analysis using both KEGG and Cistrome-GO


# MYCN (x 2) in IgHV and cox
MYCN_1.G_union <- union(IgHV.eGRNs[[12]]$up, cox.eGRNs[[12]]$up)
MYCN_2.G_union <- union(IgHV.eGRNs[[42]]$up, cox.eGRNs[[42]]$up)
length(MYCN_1.G_union)
length(MYCN_2.G_union)
MYCN_1 <- IRG.eGRNs[[12]]
MYCN_2 <- IRG.eGRNs[[42]]
MYCN_1$TF
MYCN_2$TF

library(dplyr)
MYCN.G_overlap <- intersect(MYCN_1$genes, MYCN_2$genes) %>% intersect(., MYCN_1.G_union) %>% 
  intersect(., MYCN_2.G_union)
length(MYCN.G_overlap)
MYCN.R_overlap <- intersect(MYCN_1$peaks, MYCN_2$peaks)
length(MYCN.R_overlap)

MYCN_1.G_diff <- setdiff(MYCN_1$genes, MYCN_2$genes) %>% intersect(., MYCN_1.G_union)
length(MYCN_1.G_diff)
MYCN_1.R_diff <- setdiff(MYCN_1$peaks, MYCN_2$peaks)
length(MYCN_1.R_diff)
MYCN_2.G_diff <- setdiff(MYCN_2$genes, MYCN_1$genes) %>% intersect(., MYCN_2.G_union)
length(MYCN_2.G_diff)
MYCN_2.R_diff <- setdiff(MYCN_2$peaks, MYCN_1$peaks)
length(MYCN_2.R_diff)


# JUN and STAT3 in IgHV and TP53
JUN.G_union <- union(IgHV.eGRNs[[49]]$up, TP53.eGRNs[[49]]$up)
STAT3.G_union <- union(IgHV.eGRNs[[44]]$up, TP53.eGRNs[[44]]$up)
length(JUN.G_union)
length(STAT3.G_union)
JUN <- IRG.eGRNs[[49]]
STAT3 <- IRG.eGRNs[[44]]
JUN$TF
STAT3$TF

library(dplyr)
JUN_STAT3.G_overlap <- intersect(JUN$genes, STAT3$genes) %>% intersect(., JUN.G_union) %>% 
  intersect(., STAT3.G_union)
length(JUN_STAT3.G_overlap)
JUN_STAT3.R_overlap <- intersect(JUN$peaks, STAT3$peaks)
length(JUN_STAT3.R_overlap)

JUN.G_diff <- setdiff(JUN$genes, STAT3$genes) %>% intersect(., JUN.G_union)
length(JUN.G_diff)
JUN.R_diff <- setdiff(JUN$peaks, STAT3$peaks)
length(JUN.R_diff)
STAT3.G_diff <- setdiff(STAT3$genes, JUN$genes) %>% intersect(., STAT3.G_union)
length(STAT3.G_diff)
STAT3.R_diff <- setdiff(STAT3$peaks, JUN$peaks)
length(STAT3.R_diff)


# JUN in IgHV, TP53, and cox
library(dplyr)
JUN.G_overlap <- intersect(IgHV.eGRNs[[49]]$up, TP53.eGRNs[[49]]$up) %>% 
  intersect(., cox.eGRNs[[49]]$up)
length(JUN.G_overlap)
JUN.R_overlap <- intersect(IgHV.eGRNs[[49]]$peaks, TP53.eGRNs[[49]]$peaks) %>% 
  intersect(., cox.eGRNs[[49]]$peaks)
length(JUN.R_overlap)
JUN <- IRG.eGRNs[[49]]
JUN$TF

# JUN.IgHV_G <- setdiff(IgHV.eGRNs[[49]]$up, TP53.eGRNs[[49]]$up) %>% setdiff(., cox.eGRNs[[49]]$up)
# length(JUN.IgHV_G)
# JUN.IgHV_R <- setdiff(IgHV.eGRNs[[49]]$peaks, TP53.eGRNs[[49]]$peaks) %>% setdiff(., cox.eGRNs[[49]]$peaks)
# length(JUN.IgHV_R)
# 
# JUN.TP53_G <- setdiff(TP53.eGRNs[[49]]$up, IgHV.eGRNs[[49]]$up) %>% setdiff(., cox.eGRNs[[49]]$up)
# length(JUN.TP53_G)
# JUN.TP53_R <- setdiff(TP53.eGRNs[[49]]$peaks, IgHV.eGRNs[[49]]$peaks) %>% setdiff(., cox.eGRNs[[49]]$peaks)
# length(JUN.TP53_R)
# 
# JUN.cox_G <- setdiff(cox.eGRNs[[49]]$up, TP53.eGRNs[[49]]$up) %>% setdiff(., IgHV.eGRNs[[49]]$up)
# length(JUN.cox_G)
# JUN.cox_R <- setdiff(cox.eGRNs[[49]]$peaks, TP53.eGRNs[[49]]$peaks) %>% setdiff(., IgHV.eGRNs[[49]]$peaks)
# length(JUN.cox_R)


##################################################################################
#                                                                                #
#            Prepare csv files for Cytoscape to generate networks                #
#                                                                                #
##################################################################################


# Selected eGRNs : MYCN_1, JUN and STAT3, and JUN : 12, (49, 44), 49
selected.IDs <- c(12, 44, 49)
library(pbapply)
source("/fs/ess/PCON0022/liyang/r_utilities/functions/cistrome_tools.R")
cistrome.list <- pblapply(selected.IDs, function(i) {
  eGRN <- list(peaks = IRG.eGRNs[[i]]$peaks, genes = IRG.eGRNs[[i]]$genes, 
              TF = IRG.eGRNs[[i]]$TF)
  link_peaks_to_genes(peak.obj = eGRN$peaks, gene.obj = eGRN$genes, 
                      distance = "gene")
})
qs::qsave(cistrome.list, paste0(work.dir, "Cistrome_list.qsave"))


# Prepare txt files for the first MYCN (ID = 12)
prepare_Cytoscape(gr = cistrome.list[[1]], node.table = T, 
                  path = paste0(work.dir, "MYCN"))
MYCN.df <- read.table(paste0(work.dir, "MYCN_nodes.txt"), header = T)
head(MYCN.df)
rownames(MYCN.df) <- MYCN.df$Nodes
head(MYCN.df)
MYCN.blue <- setdiff(IgHV.eGRNs[[12]]$up, cox.eGRNs[[12]]$up) # IgHV
MYCN.yellow <- setdiff(cox.eGRNs[[12]]$up, IgHV.eGRNs[[12]]$up) # Cox
# MYCN.grey <- setdiff(MYCN.df$Nodes, MYCN.blue) %>% setdiff(., MYCN.yellow)
MYCN.df <- cbind(MYCN.df, Color = "grey")
MYCN.df[MYCN.blue, "Color"] <- "blue"
MYCN.df[MYCN.yellow, "Color"] <- "yellow"
MYCN.df$Color <- factor(MYCN.df$Color, levels = c("red", "blue", "green", 
                                                  "yellow", "grey"))
head(MYCN.df)
tail(MYCN.df$Color)
write.table(MYCN.df, paste0(work.dir, "MYCN_nodes.txt"), 
            quote = F, row.names = F, sep = "\t")


# Prepare txt files for JUN and STAT3 (IDs = 49, 44)
JUN.list <- cistrome.list[[3]]
STAT3.list <- cistrome.list[[2]]
JUN_STAT3.list <- JUN.list[JUN.list$gene %in% JUN_STAT3.G_overlap]
# JUN.peaks <- GRangesToString(JUN.list)
# JUN.list <- JUN.list[JUN.peaks %in% JUN_STAT3.R_overlap]
prepare_Cytoscape(gr = JUN_STAT3.list, node.table = T, 
                  path = paste0(work.dir, "JUN_and_STAT3"))
JUN.STAT3.df <- read.table(paste0(work.dir, "JUN_and_STAT3_nodes.txt"), header = T)
head(JUN.STAT3.df)
rownames(JUN.STAT3.df) <- JUN.STAT3.df$Nodes
head(JUN.STAT3.df)
JUN.STAT3.red <- intersect(unique(JUN.STAT3.df$Nodes), IgHV.eGRNs[[49]]$up) %>% 
  intersect(., IgHV.eGRNs[[44]]$up) %>% intersect(., TP53.eGRNs[[49]]$up) %>% 
  intersect(., TP53.eGRNs[[44]]$up)
JUN.STAT3.blue <- intersect(JUN.STAT3.df$Nodes, IgHV.eGRNs[[49]]$up) %>% 
  intersect(., IgHV.eGRNs[[44]]$up) %>% setdiff(., TP53.eGRNs[[44]]$up) %>% 
  setdiff(., TP53.eGRNs[[49]]$up) # IgHV
JUN.STAT3.green <- intersect(JUN.STAT3.df$Nodes, TP53.eGRNs[[49]]$up) %>% 
  intersect(., TP53.eGRNs[[44]]$up) %>% setdiff(., IgHV.eGRNs[[44]]$up) %>% 
  setdiff(., IgHV.eGRNs[[49]]$up) # TP53
# MYCN.grey <- setdiff(MYCN.df$Nodes, MYCN.blue) %>% setdiff(., MYCN.yellow)
JUN.STAT3.df <- cbind(JUN.STAT3.df, Color = "grey")
JUN.STAT3.df[JUN.STAT3.red,]$Color <- "red"
JUN.STAT3.df[JUN.STAT3.blue,]$Color <- "blue"
JUN.STAT3.df[JUN.STAT3.green,]$Color <- "green"
JUN.STAT3.df$Color <- factor(JUN.STAT3.df$Color, levels = c("red", "blue", "green", 
                                                  "yellow", "grey"))
head(JUN.STAT3.df)
tail(JUN.STAT3.df$Color)
write.table(JUN.STAT3.df, paste0(work.dir, "JUN_STAT3_nodes.txt"), 
            quote = F, row.names = F, sep = "\t")


# Prepare txt files for JUN (ID = 49)
prepare_Cytoscape(gr = cistrome.list[[3]], node.table = T, 
                  path = paste0(work.dir, "JUN"))
JUN.df <- read.table(paste0(work.dir, "JUN_nodes.txt"), header = T)
head(JUN.df)
rownames(JUN.df) <- JUN.df$Nodes
head(JUN.df)
JUN.red <- intersect(unique(JUN.df$Nodes), IgHV.eGRNs[[49]]$up) %>% 
  intersect(., TP53.eGRNs[[49]]$up) %>% intersect(., cox.eGRNs[[49]]$up)
JUN.blue <- intersect(unique(JUN.df$Nodes), IgHV.eGRNs[[49]]$up) %>% 
  setdiff(., cox.eGRNs[[49]]$up) %>% setdiff(., TP53.eGRNs[[49]]$up) # IgHV
JUN.yellow <- intersect(unique(JUN.df$Nodes), cox.eGRNs[[49]]$up) %>% 
  setdiff(., IgHV.eGRNs[[49]]$up) %>% setdiff(., TP53.eGRNs[[49]]$up) # Cox
JUN.green <- intersect(unique(JUN.df$Nodes), TP53.eGRNs[[49]]$up) %>% 
  setdiff(., IgHV.eGRNs[[49]]$up) %>% setdiff(., cox.eGRNs[[49]]$up) # TP53
JUN.df <- cbind(JUN.df, Color = "grey")
JUN.df[JUN.red,]$Color <- "red"
JUN.df[JUN.blue,]$Color <- "blue"
JUN.df[JUN.yellow,]$Color <- "yellow"
JUN.df$Color <- factor(JUN.df$Color, levels = c("red", "blue", "green", 
                                                  "yellow", "grey"))
head(JUN.df)
tail(JUN.df$Color)
write.table(JUN.df, paste0(work.dir, "JUN_nodes.txt"), 
            quote = F, row.names = F, sep = "\t")
