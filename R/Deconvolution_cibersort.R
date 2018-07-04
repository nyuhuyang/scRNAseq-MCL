library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")
# 5.1.1 load data
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
TSNEPlot(object = MCL, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 5)+
        ggtitle("TSNE plot of major cell types")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
MCL.subsets <- SplitCells(MCL)
MCL.subsets[[3]]
MCL.patient <- as.matrix(MCL.subsets[[1]]@data)
MCL.normal <- as.matrix(MCL.subsets[[2]]@data)
MCL.patient <- rowSums(MCL.patient)
MCL.normal <- rowSums(MCL.normal)
MCL_expression <- data.frame("patient"=MCL.patient,
                             "normal"=MCL.normal)
write.table(MCL.expression,file = "./output/MCL_expression.txt",sep = "\t")
# add Gene symbol to the first line.
# remove all "

CIBERSORT <- read.csv("./output/CIBERSORT.csv")
CIBERSORT <- 
table(MCL.subsets[[1]]@ident)/ncol(MCL.subsets[[1]]@scale.data)
table(MCL.subsets[[2]]@ident)/ncol(MCL.subsets[[2]]@scale.data)
