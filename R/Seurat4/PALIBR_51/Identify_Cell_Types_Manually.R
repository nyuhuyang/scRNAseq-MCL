library(Seurat)
library(dplyr)
library(tidyr)
library(magrittr)
library(kableExtra)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
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


df_samples <- readxl::read_excel("doc/20220901_scRNAseq_info.xlsx") %>% as.data.frame
colnames(df_samples) %<>% tolower()
df_samples %<>% filter(sequence %in% "GEX") %>% filter(phase %in% "PALIBR_I")
nrow(df_samples)
df_samples$date %<>% gsub(" UTC","",.) %>% as.character()

meta.data <- readRDS(file = "output/MCL_51_20210724_meta.data_v2.rds")
colnames(meta.data) %<>% sub("Phase","cell cycle phase",.)
for(i in 1:length(df_samples$sample)){
    cells <- meta.data$orig.ident %in% df_samples$sample[i]
    print(table(cells))
    meta.data[cells,"response"] = df_samples$response[i]
}

saveRDS(meta.data,"output/MCL_51_20210724_meta.data_v2.rds")
