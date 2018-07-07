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

# 5.2 CIBERSORT
CIBERSORT <- read.csv("./output/CIBERSORT.csv")
CIBERSORT <- 
table(MCL.subsets[[1]]@ident)/ncol(MCL.subsets[[1]]@scale.data)
table(MCL.subsets[[2]]@ident)/ncol(MCL.subsets[[2]]@scale.data)

# EpiDISH
library(EpiDISH)
MCL.expression <- read.delim2("./output/MCL_expression.txt",row.names = 1)
LM22 <- read.delim2("./data/LM22.txt",row.names = 1)
MCL.expression <- data.matrix(MCL.expression) # convert dataframe to matrix
LM22 <- data.matrix(LM22) # convert dataframe to matrix
head(sapply(LM22,class))
head(sapply(MCL.expression,class))

out.l <- epidish(MCL.expression, LM22, method = "RPC") 
results_RPC <- out.l$estF
results_RPC

out.l <- epidish(MCL.expression, LM22, method = "CBS") 
results_CBS <- out.l$estF
results_CBS
#out.l <- epidish(MCL.expression, LM22, method = "CP",constraint = "inequality") 
#results_CP <- out.l$estF

#DeconRNASeq
library(DeconRNASeq)
MCL.expression <- read.delim2("./output/MCL_expression.txt",row.names = 1)
LM22 <- read.delim2("./data/LM22.txt",row.names = 1)
MCL.expression <- as.data.frame(data.matrix(MCL.expression)) # convert to dataframe 
LM22 <-  as.data.frame(data.matrix(LM22)) # convert to dataframe 
head(sapply(LM22,class)) 
head(sapply(MCL.expression,class))

results_Decon <- DeconRNASeq(MCL.expression, LM22, NULL, checksig=FALSE, 
            known.prop = FALSE, use.scale = FALSE, fig = FALSE)
rownames(results_Decon$out.all) <- colnames(MCL.expression)
results_Decon <- results_Decon$out.all
results_Decon

# xCell
library(xCell)
MCL.expression <- read.delim2("./output/MCL_expression.txt",row.names = 1)
results_xCell <- xCellAnalysis(MCL.expression)
results_xCell

# DeMixT
library("DeMixT")
MCL.expression <- read.delim2("./output/MCL_expression.txt",row.names = 1)
LM22 <- read.delim2("./data/LM22.txt",row.names = 1)
MCL.expression <- data.matrix(MCL.expression) # convert dataframe to matrix
LM22 <- data.matrix(LM22) # convert dataframe to matrix
head(sapply(LM22,class))
head(sapply(MCL.expression,class))

#Error in cbind(data.comp1, data.Y) : number of rows of matrices must match (see arg 2)
CommonGenes <- intersect(rownames(MCL.expression),rownames(LM22))
MCL.expression <- MCL.expression[CommonGenes,]
LM22 <- LM22[CommonGenes,]
dim(MCL.expression)
dim(LM22)

results_DeMixT <- DeMixT(data.Y = MCL.expression, data.comp1 = LM22, if.filter = FALSE) 
results_DeMixT <- DeMixT.S1(data.Y = MCL.expression, data.comp1 = LM22, if.filter = FALSE) 
