########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
lnames = load(file = "./data/MCL_alignment.Rda")
lnames
table(MCL@ident)
idents <- as.data.frame(table(MCL@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("T cells",
                     "CD14 Monocytes",
                     "B cells",
                     "CD8 T cells",
                     "B cells",
                     "NK T cells",
                     "CD16 Monocytes",
                     "CD8 T cells",
                     "Dendritic Cells",
                     "Myeloid cells")

MCL@ident <- plyr::mapvalues(x = MCL@ident,
                             from = old.ident.ids,
                             to = new.cluster.ids)
MCL@meta.data$conditions <- sub("_",".",MCL@meta.data$conditions)
#===========answer to 6/1's email=====
Treg.markers.to.plot <- c("FOXP3","CTLA4","PDCD1","ENTPD1","CD38","ICOS",
                     "TNFSF9","TNFRSF9")
markers.to.plot <- HumanGenes(MCL,markers.to.plot,unique=T)
SplitDotPlotGG.1(object=MCL, genes.plot = rev(Treg.markers.to.plot),
                 cols.use = c("blue","red"), x.lab.rot = T, 
                 plot.legend = T, dot.scale = 8, do.return = T, 
                 grouping.var = "conditions")
#===========answer to 6/5's email=====
markers.to.plot <- c("CD79A","CD3G","CD8A","CD40","CXCR4")
markers.to.plot <- HumanGenes(MCL,markers.to.plot,unique=T)
SplitDotPlotGG.1(object=MCL, genes.plot = rev(markers.to.plot),
                 cols.use = c("blue","red"), x.lab.rot = T, 
                 plot.legend = T, dot.scale = 8, do.return = T, 
                 grouping.var = "conditions")
cellcycle.to.plot <- c("CCND1","CCND2","CCND3","CDK4","CDK6",
                       "PCNA","SOX11")
SplitDotPlotGG.1(object=MCL, genes.plot = rev(cellcycle.to.plot),
                 cols.use = c("blue","red"), x.lab.rot = T, 
                 plot.legend = T, dot.scale = 8, do.return = T, 
                 grouping.var = "conditions")

gde.all <- FindAllMarkersAcrossConditions(MCL)
write.csv(gde.all,"./doc/MCL_patient_vs_normal.csv")
#===========answer to 6/17 & 6/21's email=====
Featureplot <- function(x,object = object,...){
        p <- FeaturePlot(object = object, 
                         reduction.use = "tsne",
                         features.plot = x, min.cutoff = NA, 
                         cols.use = c("lightgrey","blue"), pt.size = 0.5,...)
        return(p)
}
markers.to.plot <- HumanGenes(MCL,c("CD79A","CD3G","CD8A","CD40","CXCR4","CEACAM1"))
MCL.subsets <- SplitCells(MCL)
MCL.subsets[[3]]
MCL.patient <- MCL.subsets[[1]]
MCL.normal <- MCL.subsets[[2]]
# Featureplot
Featureplot(cellcycle[-5],MCL.patient) # cellcycle
Featureplot(cellcycle[-5],MCL.normal) # cellcycle
Featureplot(T_Cell[1:6],MCL.patient) 
Featureplot(T_Cell[1:6],MCL.normal) 
Featureplot(Treg[1:6],MCL.patient) # Treg
Featureplot(Treg[1:6],MCL.normal) # Treg
Featureplot(Treg[5:10],MCL.patient) # Treg
Featureplot(Treg[5:10],MCL.normal) # Treg
Featureplot(c(CD4_Naive_T,"CD8A"),MCL.patient) 
Featureplot(c(CD4_Naive_T,"CD8A"),MCL.normal) 
Featureplot(NK,MCL.patient) # Treg
Featureplot(NK,MCL.normal) # Treg

Featureplot(B_Cell[1:6],MCL.patient) 
Featureplot(B_Cell[1:6],MCL.normal) 
Featureplot(c("CD1C","CD5","CD27","CD38","IL4R","CD40"),MCL.patient) 
Featureplot(c("CD1C","CD5","CD27","CD38","IL4R","CD40"),MCL.normal) 
Featureplot(c(Neutrophil,B_Cell[1:2]),MCL.patient) 
Featureplot(c(Neutrophil,B_Cell[1:2]),MCL.normal)

Featureplot(markers.to.plot,MCL.patient) 
Featureplot(markers.to.plot,MCL.normal)

Featureplot(c(CD14_Monocytes[-4],CD16_Monocytes),MCL.patient) 
Featureplot(c(CD14_Monocytes[-4],CD16_Monocytes),MCL.normal)

Featureplot(Stem_cell,MCL.patient) 
Featureplot(Stem_cell,MCL.normal)
