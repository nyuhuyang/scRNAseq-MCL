# ######################################################################
library(Seurat)
library(dplyr)
library(kableExtra)
library(magrittr)
library(ggplot2)
library(cowplot)
library(fgsea)
library(tibble)
library(ggsci)
library(fgsea)
library(openxlsx)
library(eulerr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

###############################
#### Fig 4
###############################
path = "Yang/PALIBR/Fig. 6/"
(load(file="data/MCL_41_harmony_20200225.Rda"))
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
object@meta.data[object$SCT_snn_res.0.8 %in% c(5,16), "cell.types"] = "Monocytes"
Idents(object) = "cell.types"
object %<>% subset(idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)
table(Idents(object))

samples = c("N01","PtU01","PtU02","PtU04",#"PtU03",
            "Pt17_2","Pt17_7","Pt17_31",#"Pt17_LN1",
            "Pt25_1","Pt25_24","Pt25_25Pd"#,"Pt25_SB1","Pt25_1_8","Pt25_AMB25Pd",
            #"Pt28_LN1","Pt28_1","Pt28_4","Pt28_28",
            #"Pt11_LN1", "Pt11_1","Pt11_14","Pt11_28",
            #"PtB13_Ibp","PtB13_Ib1","PtB13_IbR"
            )
samples %<>% factor(levels = as.character(samples))

Idents(object) = "cell.types"
gene_list <- list("T_cells_CD4+" = c("IL32", "CD52", "TXNIP", "EEF1G", "EML4","TSC22D3"),#"TNFAIP3",
                  "T_cells_CD8+" = c("IL32", "CD52", "LTB", "IFITM1", "JUN", "JUNB", "RPS26",
                             "EEF1G", "EML4", "TSC22D3"))#"PIK3R1", "TNFAIP3",
cell.types <- c("T_cells_CD8+","T_cells_CD4+")

#' adjust hex color scheme acorrding to the range, put white color at middle
#' @example skewed_corlor(sns.RdBu_r, Min = -1, Max = 2)
skewed_corlor <- function(hex_color, Min, Max){
        #if(Min > 0) return(hex_color[7:12])
        #if(Max < 0) return(hex_color[1:6])
        if(Min < 0 & Max > 0) {
                l = length(hex_color)
                Range = max(Max, -Min)*2
                remove <- (Min + Max) / Range * l
                if(remove > 0) return(hex_color[ceiling(remove):l])
                if(remove < 0) return(hex_color[1:(l+min(-1,remove)+1)])
        }
}
#sns.color_palette("RdBu_r", 15).as_hex()

# average scale before
for(i in seq_along(cell.types)){
        #sns.RdBu <- if(i ==1) sns.RdBu_r_199 else sns.RdBu_r_15
        features = gene_list[[cell.types[i]]]
        sub_object <- subset(object, idents = gsub("_CD",":CD",cell.types[i]))
        Idents(sub_object) = "orig.ident"
        sub_object %<>% subset(idents = samples[samples %in% Idents(sub_object)])
        sub_object %<>% ScaleData(features = features)
        exp = AverageExpression(sub_object[features,],
                                assays = "SCT") %>% .$SCT
        #if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        exp %<>% t %>% scale %>% t
        S = samples[samples %in% colnames(exp)]
        DoHeatmap.matrix(exp[,as.character(S)], features = features,
                         group.by = S,size = 4,angle = 45,
                         draw.lines =F, raster = FALSE,
                         colors = skewed_corlor(sns.RdBu_r_199,
                                                Min = min(exp[,as.character(S)]),
                                                Max = max(exp[,as.character(S)])),
                         width=ifelse(i == 1,8,6), height=ifelse(i == 1,6.5,3),res=600,
                         no.legend = T,
                         cex.row=10,
                         group.bar.height = 0,
                         do.print = T,
                         file.name = paste0("Heatmap_",cell.types[i],"_","average_scale.jpeg"),
                         save.path = paste0(path,"Heatmaps"))
        print(min(exp[,as.character(S)]))
        print(max(exp[,as.character(S)]))
}

# scale before average
for(i in seq_along(cell.types)){
        #sns.RdBu <- if(i ==1) sns.RdBu_r_199 else sns.RdBu_r_15
        features = gene_list[[cell.types[i]]]
        sub_object <- subset(object, idents = gsub("_CD",":CD",cell.types[i]))
        Idents(sub_object) = "orig.ident"
        sub_object %<>% subset(idents = samples[samples %in% Idents(sub_object)])
        sub_object %<>% ScaleData(features = features)
        exp = AverageExpression(sub_object[features,],
                                assays = "SCT",slot = "scale.data") %>% .$SCT
        #if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

        S = samples[samples %in% colnames(exp)]
        DoHeatmap.matrix(exp[,as.character(S)], features = features,
                         group.by = S,size = 4,angle = 45,
                         draw.lines =F, raster = FALSE,
                         colors = skewed_corlor(sns.RdBu_r_199,
                                                Min = min(exp[,as.character(S)]),
                                                Max = max(exp[,as.character(S)])),
                         width=ifelse(i == 1,8,6), height=ifelse(i == 1,6.5,3),res=600,
                         no.legend = T,
                         cex.row=10,
                         group.bar.height = 0,
                         do.print = T,
                         file.name = paste0("Heatmap_",cell.types[i],"_","scale_average.jpeg"),
                         save.path = paste0(path,"Heatmaps"))
        print(min(exp[,as.character(S)]))
        print(max(exp[,as.character(S)]))
}

for(i in seq_along(cell.types)){
        #sns.RdBu <- if(i ==1) sns.RdBu_r_199 else sns.RdBu_r_15
        features = gene_list[[cell.types[i]]]
        sub_object <- subset(object, idents = gsub("_CD",":CD",cell.types[i]))
        Idents(sub_object) = "orig.ident"
        sub_object %<>% subset(idents = samples[samples %in% Idents(sub_object)])
        sub_object %<>% ScaleData(features = features)
        exp = AverageExpression(sub_object[features,],
                                assays = "SCT") %>% .$SCT
        #if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

        S = samples[samples %in% colnames(exp)]
        DoHeatmap.matrix(exp[,as.character(S)], features = features,
                         group.by = S,size = 4,angle = 45,
                         draw.lines =F, raster = FALSE,
                         colors = sns.RdBu_r_199,
                         width=8, height=6.5,res=600,
                         no.legend = F,
                         cex.row=10,
                         group.bar.height = 0,
                         do.print = T,
                         save.path = paste0(path,"Heatmaps/PB_only/average/Doheatmap_",cell.types[i]))
        print(min(exp[,as.character(S)]))
        print(max(exp[,as.character(S)]))
        if(i == 2){
                DoHeatmap.matrix(exp[,as.character(S)], features = features,
                                 group.by = S,size = 4,angle = 45,
                                 draw.lines =F, raster = FALSE,
                                 colors = sns.RdBu_r_199,
                                 width=ifelse(i == 1,8,6), height=ifelse(i == 1,6.5,3),res=600,
                                 no.legend = F,
                                 cex.row=10,
                                 group.bar.height = 0,
                                 do.print = T,
                                 save.path = paste0(path,"Heatmaps/PB_only/average/Doheatmap_",cell.types[i],"_noLegend"))

        }
}

object@meta.data$cell.types %<>% gsub(":CD","_CD",.)
object@meta.data$cell.types %<>% gsub("MCL|B_cells","MCL_B_cells",.)
Idents(object) = "cell.types"
# markers with Log2FC > 0.1, adjust p-value < 0.05
cell.types <- c("T_cells_CD8+","T_cells_CD4+","NK_cells","Monocytes","MCL_B_cells")
for(i in seq_along(cell.types)){
        file_list <- list.files(path= paste0(path,"DE_results/", cell.types[i]),
                                pattern = "csv", full.names = TRUE, recursive = F)
        file_list = file_list[grep("-N01.csv",file_list)]
        length(file_list)
        DE_res_list <- lapply(file_list, read.csv)
        DE_res <- bind_rows(DE_res_list)
        DE_res <- DE_res[DE_res$avg_logFC > 0.1 & DE_res$p_val_adj < 0.05 & !(DE_res$cluster %in% paste0("N0",1:4)),]
        features <- unique(DE_res$gene)
        print(length(features))
        sub_object <- subset(object, idents = cell.types[i])
        Idents(sub_object) = "orig.ident"
        samples = unique(c(DE_res$cluster, paste0("N0",1:4)))
        print(length(samples))
        sub_object %<>% subset(idents = samples[samples %in% Idents(sub_object)])
        sub_object %<>% ScaleData(features = features)
        exp = AverageExpression(sub_object[features,],
                                assays = "SCT") %>% .$SCT
        exp %<>% t %>% scale %>% t

        samples = samples[samples %in% colnames(exp)]
        save.path = paste0(path,"Heatmaps/",cell.types[i],"/")
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

        jpeg(paste0(save.path,"heatmap2_",cell.types[i],".jpeg"), units="in", width=10, height=7,res=600)
        heatmap.2(as.matrix(exp),
                  breaks = seq(-2,2,length.out = 300),
                  cexRow = 0.5,
                  #col = sns.RdBu_r,
                  dendrogram = "both",
                  margins = c(7,5),
                  col = bluered(299),
                  key.xlab = "scale log nUMI",
                  Colv = T,
                  Rowv = T,
                  scale= "none",
                  trace = "none",
                  density.info="none",
                  main = cell.types[i])
        dev.off()
}

# markers with Log2FC > 0.1, adjust p-value < 0.05
cell.types <- c("Monocytes","NK_cells")
for(i in seq_along(cell.types)){
        file_list <- list.files(path= paste0(path, cell.types[i],"_DE"),
                                pattern = "csv", full.names = TRUE, recursive = TRUE)
        file_list = file_list[grep("-N01.csv",file_list)]
        DE_res_list <- lapply(file_list, read.csv)
        DE_res <- bind_rows(DE_res_list)
        DE_res <- DE_res[DE_res$avg_logFC > 0.1 & DE_res$p_val_adj < 0.05 & DE_res$cluster != "N01",]
        features <- unique(DE_res$gene)
        print(length(features))
        sub_object <- subset(object, idents = cell.types[i])
        Idents(sub_object) = "orig.ident"
        samples = c(unique(DE_res$cluster),"N01")
        sub_object %<>% subset(idents = samples[samples %in% Idents(sub_object)])
        sub_object %<>% ScaleData(features = features)
        exp = AverageExpression(sub_object[features,],
                                assays = "SCT") %>% .$SCT
        exp %<>% t %>% scale %>% t

        samples = samples[samples %in% colnames(exp)]
        save.path = paste0(path,"Heatmaps/",cell.types[i],"/")
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

        jpeg(paste0(save.path,"heatmap2_",cell.types[i],".jpeg"), units="in", width=10, height=7,res=600)
        heatmap.2(as.matrix(exp),
                  breaks = seq(-2,2,length.out = 300),
                  cexRow = 0.5,
                  #col = sns.RdBu_r,
                  dendrogram = "both",
                  margins = c(7,5),
                  col = bluered(299),
                  key.xlab = "scale log nUMI",
                  Colv = T,
                  Rowv = T,
                  scale= "none",
                  trace = "none",
                  density.info="none",
                  main = cell.types[i])
        dev.off()
}


# markers with Log2FC > 0.1, adjust p-value < 0.05
cell.types <- c("Monocytes","NK_cells")
for(i in seq_along(cell.types)){
        file_list <- list.files(path= paste0(path, cell.types[i],"_DE"),
                                pattern = "csv", full.names = TRUE, recursive = TRUE)
        file_list = file_list[grep("-N01.csv",file_list)]
        DE_res_list <- lapply(file_list, read.csv)
        DE_res <- bind_rows(DE_res_list)
        DE_res <- DE_res[DE_res$avg_logFC > 0.1 & DE_res$p_val_adj < 0.05 & DE_res$cluster != "N01",]
        features <- unique(DE_res$gene)
        print(length(features))
        sub_object <- subset(object, idents = cell.types[i])
        Idents(sub_object) = "orig.ident"
        samples = c(unique(DE_res$cluster),"N01")
        sub_object %<>% subset(idents = samples[samples %in% Idents(sub_object)])
        sub_object %<>% ScaleData(features = features)
        exp = AverageExpression(sub_object[features,],
                                assays = "SCT") %>% .$SCT
        exp %<>% t %>% scale %>% t

        samples = samples[samples %in% colnames(exp)]
        save.path = paste0(path,"Heatmaps/",cell.types[i],"/")
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

        jpeg(paste0(save.path,"heatmap2_",cell.types[i],"_long.jpeg"), units="in", width=10, height=70,res=600)
        heatmap.2(as.matrix(exp),
                  breaks = seq(-2,2,length.out = 300),
                  cexRow = 0.5,
                  #col = sns.RdBu_r,
                  dendrogram = "both",
                  margins = c(7,5),
                  col = bluered(299),
                  key.xlab = "scale log nUMI",
                  Colv = T,
                  Rowv = T,
                  scale= "none",
                  trace = "none",
                  density.info="none",
                  main = cell.types[i])
        dev.off()
}

# markers with Log2FC > 0.1, adjust p-value < 0.05
cell.types <- c("MCL_B_cells")
for(i in seq_along(cell.types)){
        file_list <- list.files(path= paste0(path, cell.types[i],"_DE"),
                                pattern = "csv", full.names = TRUE, recursive = F)
        file_list = file_list[grep("-N01.csv",file_list)]
        DE_res_list <- lapply(file_list, read.csv)
        DE_res <- bind_rows(DE_res_list)
        DE_res <- DE_res[DE_res$avg_logFC > 0.1 & DE_res$p_val_adj < 0.05 & DE_res$cluster != "N01",]
        features <- unique(DE_res$gene)
        print(length(features))
        sub_object <- subset(object, idents = c("MCL","B_cells"))
        Idents(sub_object) = "orig.ident"
        samples = c(unique(DE_res$cluster),"N01")
        sub_object %<>% subset(idents = samples[samples %in% Idents(sub_object)])
        sub_object %<>% ScaleData(features = features)
        exp = AverageExpression(sub_object[features,],
                                assays = "SCT") %>% .$SCT
        exp %<>% t %>% scale %>% t

        samples = samples[samples %in% colnames(exp)]
        save.path = paste0(path,"Heatmaps/",cell.types[i],"/")
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

        jpeg(paste0(save.path,"heatmap2_",cell.types[i],"_long.jpeg"), units="in", width=10, height=35,res=300)
        heatmap.2(as.matrix(exp),
                  breaks = seq(-2,2,length.out = 300),
                  cexRow = 0.5,
                  #col = sns.RdBu_r,
                  dendrogram = "both",
                  margins = c(7,5),
                  col = bluered(299),
                  key.xlab = "scale log nUMI",
                  Colv = T,
                  Rowv = T,
                  scale= "none",
                  trace = "none",
                  density.info="none",
                  main = cell.types[i])
        dev.off()
}
