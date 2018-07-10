library(Seurat)
library(tibble)
library(dplyr)
source("./R/Seurat_functions.R")
library(reshape2)
# 5.1.1 load data
# Rename ident
lnames = load(file = "./data/MCL_alignment.Rda")
lnames
table(MCL@ident)
idents <- as.data.frame(table(MCL@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("T.cells",
                     "Monocytes.CD14",
                     "B.cells",
                     "T.cells.CD8",
                     "B.cells",
                     "NK.T.cells",
                     "Monocytes.CD16",
                     "T.cells.CD8",
                     "Dendritic.cells",
                     "Macrophages")
MCL@ident <- plyr::mapvalues(x = MCL@ident,
                             from = old.ident.ids,
                             to = new.cluster.ids)
MCL@meta.data$conditions <- sub("_",".",MCL@meta.data$conditions)
TSNEPlot(object = MCL, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 5)+
        ggtitle("TSNE plot of major cell types")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
freq_table <- prop.table(x = table(MCL@ident, MCL@meta.data[, "conditions"]),
                         margin = 2)
freq_table
Seurat_freq <- matrix(freq_table,ncol=2,
                      dimnames=list(rownames(freq_table),
                                    c("patient.Seurat","normal.Seurat")))
Seurat_freq <- data.frame(Seurat_freq)

# generate MCL_expression
MCL.subsets <- SplitCells(MCL)
MCL.subsets[[3]]
MCL.patient <- as.matrix(MCL.subsets[[1]]@raw.data)
MCL.normal <- as.matrix(MCL.subsets[[2]]@raw.data)
MCL.patient<- MCL.patient[,MCL.subsets[[1]]@cell.names]
MCL.normal<- MCL.normal[,MCL.subsets[[2]]@cell.names]
MCL_patient <- rowSums(MCL.patient)
MCL_normal <- rowSums(MCL.normal)
MCL_expression <- data.frame("patient"=MCL_patient,
                             "normal"=MCL_normal)
write.table(MCL_expression,file = "./output/MCL_expression.txt",sep = "\t")
# add "Gene symbol" to the first line.
# remove all "


# CIBERSORT
CIBERSORT <- read.csv("./output/CIBERSORT.csv",row.names = 1)
CIBERSORT <- CIBERSORT[,-((ncol(CIBERSORT)-2):ncol(CIBERSORT))]
t(CIBERSORT)[1:5,]
# 5.1.2 EpiDISH
library(EpiDISH)
MCL.expression <- read.delim2("./output/MCL_expression.txt",row.names = 1)
LM22 <- read.delim2("./data/LM22.txt",row.names = 1)
MCL.expression <- data.matrix(MCL.expression) # convert dataframe to matrix
LM22 <- data.matrix(LM22) # convert dataframe to matrix
head(sapply(LM22,class))
head(sapply(MCL.expression,class))

out.l <- epidish(MCL.expression, LM22, method = "RPC", maxit = 100) 
results_RPC <- out.l$estF
t(results_RPC)[1:5,]

out.l <- epidish(MCL.expression, LM22, method = "CBS") 
results_CBS <- out.l$estF
t(results_CBS)[1:5,]

#out.l <- epidish(avdata.m=MCL.expression, ref.m=LM22, method = "CP",constraint = "equality") 
#results_CP <- out.l$estF

# 5.1.3 DeconRNASeq
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
t(results_Decon)[1:5,]

# 5.1.4 xCell
library(xCell)
par(mar=c(14,3,3,3))
MCL.expression <- read.delim2("./output/MCL_expression.txt",row.names = 1)
results_xCell <- xCellAnalysis(expr=MCL.expression)
results_xCell1 <- results_xCell[-((nrow(results_xCell)-2):nrow(results_xCell)),]
results_xCell1 <- melt(results_xCell1)
colnames(results_xCell1) <- c("Cell_Types","patient_vs_normal","p_value")
results_xCell1 <- results_xCell1[results_xCell1$p_value<0.01,]
ggplot(results_xCell1,aes(x = Cell_Types, y = p_value,
                          group=patient_vs_normal))+
        geom_col(aes(fill = Cell_Types))+
        facet_grid(.~patient_vs_normal)+
        scale_colour_gradient(limits=c(0, 0.1),low="black", high="black") +
        scale_y_reverse()+
        #geom_text(aes(label = p_value), vjust = -.5)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position="none") #remove legend
         
# 5.1.5 DeMixT
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

results_DeMixT <- DeMixT(data.Y = MCL.expression, data.comp1 = LM22,
                         if.filter = FALSE) 
results_DeMixT <- DeMixT.S1(data.Y = MCL.expression, data.comp1 = LM22, if.filter = FALSE) 

# 5.2 Summarize all deconvolution results
rownames(CIBERSORT) <- paste0(rownames(CIBERSORT),".CIBERSORT")
rownames(results_RPC) <- paste0(rownames(results_RPC),".epidish_RPC")
rownames(results_CBS) <- paste0(rownames(results_CBS),".epidish_CBS")
rownames(results_Decon) <- paste0(rownames(results_Decon),".DeconRNASeq")

results <- data.frame(t(CIBERSORT),t(results_RPC),t(results_CBS),t(results_Decon))
results_all <- full_join(rownames_to_column(results), 
                         rownames_to_column(Seurat_freq), by = "rowname")
results_all$Cell_Types <- gsub('\\..*', '',results_all$rowname)
results_all <- results_all[,-1]
results_all[is.na(results_all)] <- 0
results_short <- aggregate(. ~ Cell_Types, data = results_all, FUN=sum)
rownames(results_short) <- results_short$Cell_Types
head(results_short)
#results_short <- results_short[,-1]
#boxplot(t(results_short))

results_short2 <- melt(results_short, value.name = "percentage")
results_short2$variable <- as.character(results_short2$variable)

Catalog <- matrix(unlist(strsplit(results_short2$variable,".",fixed = T)),
               ncol=2,byrow = T)
colnames(Catalog) <- c("patient_vs_normal","software")
results_short2 <- cbind.data.frame(results_short2,Catalog)
head(results_short2)
results_short2$Cell_Types <- as.factor(results_short2$Cell_Types)
ggplot(results_short2,aes(x = Cell_Types, y = percentage,
                                  group=software))+
        geom_col(aes(fill = Cell_Types))+
        facet_grid(software~patient_vs_normal)+
        scale_y_continuous(labels = scales::percent_format())+
        geom_text(aes(label = scales::percent(percentage)), vjust = -.5)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
