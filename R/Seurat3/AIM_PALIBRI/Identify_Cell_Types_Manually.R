library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# ======== 2.1 =========== test with known markers==================
(load(file = "data/MCL_AIM_58_20201009.Rda"))
DefaultAssay(object) <- "SCT"
Idents()
# global
features <- c("CD19","CCND1","SOX11",
              "CD3D","CD4","CD8A",
              "MS4A7","CD14","FCGR1A",
              "GNLY","KLRC1","NCAM1")
FeaturePlot.1(object,features = features, pt.size = 0.005,
              cols = c("lightgrey", "red"),
              alpha = 1,reduction = "tsne",
              unique.name = "cell.types",label = F,
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 3,
              units = "in",width=9, height=12, no.legend = T)
QC <- c("percent.mt","nCount_SCT","nFeature_SCT")
FeaturePlot.1(object,features = QC, pt.size = 0.005,
              cols = c("lightgrey", "red"),
              alpha = 1,reduction = "tsne",
              unique.name = F,label = F,
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 3,
              units = "in",width=9, height=4, no.legend = T)
cc <- c("CCND1","CDK4","RB1","E2F1","MCM7","CCNB2")
FeaturePlot.1(object,features = cc, pt.size = 0.005,
              cols = c("lightgrey", "red"),
              alpha = 1,reduction = "tsne",
              unique.name = "cell.types",label = F,
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 2,
              units = "in",width=9, height=9, no.legend = T)

write.csv(as.data.frame.table(table(object$label.fine)),file = paste0(path,"label.fine.csv"))
object$label = object$label.fine
object$label %<>% gsub("T cells, CD4+.*","T cells, CD4+",.)
object$label %<>% gsub("T cells, CD8+.*","T cells, CD8+",.)
object$label %<>% gsub("Monocytes.*","Monocytes",.)
object$label %<>% gsub("B cells,.*","B cells",.)
write.csv(as.data.frame.table(table(object$label)),file = paste0(path,"label.csv"))

labels = c("T cells, CD4+","T cells, CD8+","NK cells","B cells","MCL","Monocytes")
Idents(object) = "label"

exp_list <- vector("list", length = length(labels))
names(exp_list) = labels
for(i in seq_along(labels)){
    sub_object <- subset(object, idents = labels[i])
    Idents(sub_object) = "orig.ident"
    cell.number.df <- as.data.frame.table(table(sub_object$orig.ident))
    cell.number = cell.number.df$Freq
    names(cell.number) = cell.number.df$Var1
    exp = AverageExpression(sub_object, assays = "SCT")
    exp_list[[labels[i]]] = rbind("cell.number" = cell.number, exp$SCT)
    svMisc::progress(i/length(exp_list)*100)

}

openxlsx::write.xlsx(exp_list,
                     file = paste0(path,"Expression_Cell.types.xlsx"),
                     colNames = TRUE, rowNames = TRUE,
                     borders = "surrounding",colWidths = c(NA, "auto", "auto"))

for(i in seq_along(exp_list)){
    write.table(exp_list[[i]],file = paste0(path,"Expression_",names(exp_list)[i],".txt"),
                sep = "\t",row.names = TRUE)
    svMisc::progress(i/length(exp_list)*100)
}
