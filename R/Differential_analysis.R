########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(SingleR)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
lname1 = load(file = "./data/MCL_alignment20180525.Rda")
lnames
table(MCL@ident)
.FindIdentLabel(MCL)
B_cell_id <- grepl("B_cell",MCL@meta.data$singler1sub)
B_cell_id <- MCL@cell.names[B_cell_id]
B_cell <- SubsetData(MCL, cells.use = B_cell_id)
B_cell <- SetAllIdent(B_cell, id = "orig.ident")
B_cell_markers <- FindAllMarkers.UMI(B_cell, test.use = "MAST")
write.csv(B_cell_markers,paste0(path,"B_cell_DE.csv"))
jpeg(paste0(path,"Bcell_heatmap.jpeg"), units="in", width=10, height=7,
     res=600)
DoHeatmap.1(B_cell,B_cell_markers,Top_n = 50,ident.use = "B cells",
            group.label.rot =F,cex.row = 6)
dev.off()
