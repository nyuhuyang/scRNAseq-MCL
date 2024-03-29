# conda activate r4.0.3
library(Seurat)
library(magrittr)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame
###############################
pred = readRDS("output/MCL_SCT_51_20210724_singleR_pred.rds")
object = readRDS(file = "data/MCL_SCT_51_20210724.rds")

singlerDF = data.frame("label.fine" = pred$pruned.labels,
                       row.names = rownames(pred))
table(rownames(pred) == rownames(object@meta.data))
table(is.na(singlerDF$label.fine))
singlerDF$label.fine[is.na(singlerDF$label.fine)]= "unknown"

##############################
# adjust cell label
##############################
# combine cell types
singlerDF$cell.types = singlerDF$label.fine
grep("B-cells|B cells, PB",singlerDF$cell.types, value = T) %>% table
singlerDF[grep("B-cells|B cells",singlerDF$cell.types),"cell.types"] ="B_cells"
singlerDF[grep("CD4+",singlerDF$cell.types),"cell.types"] ="T_cells:CD4+"
singlerDF[grep("CD8+",singlerDF$cell.types),"cell.types"] ="T_cells:CD8+"
singlerDF$cell.types %<>% gsub("Tregs","T_cells:regs",.)
singlerDF$cell.types %<>% gsub("MEP|CLP|HSC|CMP|GMP|MPP","HSC/progenitors",.)
singlerDF$cell.types %<>% gsub("DC|Macrophages|Macrophages M1","other Myeloid cells",.)
singlerDF$cell.types %<>% gsub("Eosinophils|Megakaryocytes","other Myeloid cells",.)
singlerDF$cell.types %<>% gsub("Adipocytes|Fibroblasts|mv Endothelial cells|Keratinocytes","Nonhematopoietic cells",.)

# reduce false positive results (B cells are labeled as MCL in normal samples)
# and false negative results (MCL cells are labeled as B cells in MCL samples)
# singler1sub false negative results  =========
CCND1 = FetchData(object,"CCND1")
singlerDF$CCND1 = CCND1$CCND1
singlerDF[(singlerDF$CCND1 >0 & singlerDF$cell.types %in% "B_cells"),"cell.types"] = "MCL"
# cell.types false positive results  ========
table(singlerDF$cell.types, object@meta.data$orig.ident) %>% kable %>% kable_styling()
normal_cells <- object$orig.ident %in% c("N01","N02","N03","N04") %>% rownames(singlerDF)[.]
singlerDF[normal_cells,"cell.types"] %<>% gsub("MCL","B_cells",.)
singlerDF$CCND1 = NULL

##############################
# process color scheme
##############################
singlerDF$cell.types.colors = singlerDF$cell.types
singlerDF$cell.types.colors %<>% plyr::mapvalues(from = c("B_cells","Erythrocytes","HSC/progenitors",
                                                          "Macrophages","MCL","Monocytes",
                                                          "NK cells","Nonhematopoietic cells","other Myeloid cells",
                                                          "Plasma cells","T_cells:CD4+","T_cells:CD8+",
                                                          "T_cells:regs","unknown"),
                                                 to = c("#E6AB02", "#ff0000", "#6A3D9A",
                                                        "#FB9A99", "#2055da", "#ADDFEE",
                                                        "#A65628","#B3B3B3","#FDDAEC",
                                                        "#1B9E77","#B3DE69","#F0027F",
                                                        "#7570B3","#F2F2F2"))
table(colnames(object) == rownames(singlerDF))
object@meta.data %<>% cbind(singlerDF)


lapply(c(TSNEPlot.1,UMAPPlot.1), function(fun)
    fun(object = object, label = T, label.repel = T,group.by = "cell.types",
        no.legend = T,
        pt.size = 0.1,label.size = 3,
        do.print = T,do.return = F,
        title ="labeling by blue_encode and MCL RNA-seq"))

saveRDS(object, file = "data/MCL_SCT_51_20210724.rds")


# by barplot
cell_Freq <- table(object$label.fine) %>% as.data.frame
cell_Freq$Percent <- prop.table(cell_Freq$Freq) %>% round(digits = 2) %>% scales::percent()
cols = ExtractMetaColor(object)
cell_Freq$cols = cols[cell_Freq$Var1]
cell_Freq = cell_Freq[order(cell_Freq$Var1),]

cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=6, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,x.text.angle = 90,
          ylab = "Cell Number",
          label = cell_Freq$Percent,
          lab.size = 3,
          sort.val = "desc",
          width = 1, size = 0.5,
          title = "Numbers of cell types in total 6 samples")+NoLegend()+
    theme(plot.title = element_text(hjust = 0.5,size=15),
          axis.text.x = element_text(vjust = 0.5))
dev.off()
