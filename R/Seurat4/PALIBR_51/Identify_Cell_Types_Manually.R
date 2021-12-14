library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# ======== 2.1 =========== test with known markers==================
object = readRDS(file = "data/MCL_SCT_51_20210724.rds")
reductions = readRDS("output/20210723/reductions_npcs=100_perplexity=175.rds")
barcode1 = gsub("-.*","",rownames(reductions$tsne@cell.embeddings))
barcode2 = gsub(".*-","",rownames(object@reductions$tsne@cell.embeddings))
table(barcode1 == barcode2)
rownames(reductions$tsne@cell.embeddings) = rownames(object@reductions$tsne@cell.embeddings)
object@reductions$tsne = reductions$tsne

# global
features <- c("CD19","CCND1","PCNA",
              "CD3D","CD4","CD8A",
              "MS4A7","CD14","FCGR3A",
              "GNLY","KLRC1","NCAM1")
FeaturePlot.1(object,features = features, pt.size = 0.005,
              cols = c("lightgrey", "red"),
              alpha = 1,reduction = "tsne",
              unique.name = "cell.types",label = F,
              threshold = 1, text.size = 20, border = T,
              file.name="TSNE_51_markers.jpeg",do.print = T, do.return = F,ncol = 3,
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

TSNEPlot.1(object, group.by = "seurat_clusters",label =T,do.print = T)
meta.data = object@meta.data
meta.data$cell.types %<>% gsub("Monocytes", "Monocytes:CD14+",.)
meta.data[meta.data$cell.types %in% "Monocytes:CD14+" & meta.data$SCT_snn_res.0.8 %in% 8,
          "cell.types"] = "Monocytes:CD16+"
meta.data[meta.data$SCT_snn_res.0.8 %in% 24 & meta.data$tSNE_2 > 20,
          "cell.types"] = "Monocytes:CD16+"

meta.data$cell.types.colors = meta.data$cell.types
meta.data$cell.types.colors %<>% plyr::mapvalues(from = c("B_cells","Erythrocytes","HSC/progenitors",
                                                          "MCL","Monocytes:CD14+","Monocytes:CD16+",
                                                          "NK cells","Nonhematopoietic cells","other Myeloid cells",
                                                          "Plasma cells","T_cells:CD4+","T_cells:CD8+",
                                                          "T_cells:regs","unknown"),
                                                 to = c("#E6AB02", "#ff0000", "#6A3D9A",
                                                        "#2055da", "#ADDFEE","#FB9A99",
                                                        "#A65628","#B3B3B3","#FDDAEC",
                                                        "#1B9E77","#B3DE69","#F0027F",
                                                        "#7570B3","#F2F2F2"))
df_samples <- readxl::read_excel("doc/20210715_scRNAseq_info.xlsx",sheet = "fastq")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples %<>% filter(sequence %in% "GEX") %>% filter(phase %in% "PALIBR_I") %>%
    filter(sample != "Pt11_31")
meta.data$orig.ident %<>% factor(levels = df_samples$`sample`)
object@meta.data = meta.data

saveRDS(object, file = "data/MCL_SCT_51_20210724.rds")
object[['RNA']] <- NULL
object[['integrated']] <- NULL
format(object.size(object),unit = "GB")


saveRDS(object, file = "data/MCL_SCT_51_20210724.rds")

MCL = subset(object, subset =  UMAP_1 > 0
                        & tSNE_2 < tSNE_1 +1
                        & UMAP_2 <10
                        & Doublets == "Singlet"
                        & cell.types  %in% c("B_cells","MCL")
)

MCL %<>% FindNeighbors(reduction = "tsne",dims = 1:2,k.param = 20)

save.path <- paste0(path,"serial_resolutions_tsne/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

resolutions = c(seq(0.001,0.009, by = 0.001),seq(0.01,0.09, by = 0.01))
for(i in 1:length(resolutions)){
    MCL %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    TSNEPlot.1(MCL, group.by=paste0("SCT_snn_res.",resolutions[i]),
               split.by = paste0("SCT_snn_res.",resolutions[i]),pt.size = 0.3,label = T,
               label.repel = T,alpha = 0.9,
               do.return = F,
               no.legend = T,label.size = 4, repel = T,
               title = paste("res =",resolutions[i]),
               do.print = T, save.path = save.path)
    Progress(i,length(resolutions))
}

TSNEPlot.1(MCL, group.by=paste0("SCT_snn_res.",0.006),pt.size = 0.3,label = T,
           label.repel = T,alpha = 0.9,
           do.return = F,
           no.legend = T,label.size = 4, repel = T,
           title = paste("res =",0.006),
           do.print = T, save.path = save.path)

MCL$X4cluster = as.integer(MCL$SCT_snn_res.0.006)
MCL = subset(MCL, subset =  X4cluster != 5)

MCL$X4cluster %<>% gsub("6|7|8|9|10","1",.)
MCL$SCT_snn_res.0.006 %<>% as.integer() %>% as.factor()

object$X4cluster = object$cell.types
object@meta.data[rownames(MCL@meta.data),"X4cluster"] = MCL$X4cluster
saveRDS(object, file = "data/MCL_SCT_51_20210724.rds")



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

# Extend Data
Idents(object) = "cell.types"
object %<>% sortIdent()

cell_Freq <- table(Idents(object)) %>% as.data.frame
cell_Freq$Percent <- round(prop.table(cell_Freq$Freq),digits = 3) %>%
    scales::percent()

cell_Freq = cell_Freq[order(cell_Freq$Var1),]
df_samples <- readxl::read_excel("doc/singler.colors.xlsx")
cell_Freq$col = plyr::mapvalues(cell_Freq$Var1,
                                from = na.omit(df_samples$Cell.types),
                                to = na.omit(df_samples$singler.color2))
cell_Freq$col %<>% as.character
cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")
cell_Freq$Cell_Type %<>% gsub("_"," ",.)


jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=6, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,x.text.angle = 45,
          ylab = "Cell Number",
          label = cell_Freq$Percent,
          sort.val = "desc",
          width = 1, size = 0.5,
          title = "Numbers of major cell types in total 74 samples")+NoLegend()+
    theme(plot.title = element_text(hjust = 0.5,size=15))+
    scale_y_continuous(expand = c(0, 0), limits = c(0,max(cell_Freq$Cell_Number)+4000))
dev.off()
