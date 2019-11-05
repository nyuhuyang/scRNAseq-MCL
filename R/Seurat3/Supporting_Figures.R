########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
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
library(ggpubr)
source("../R/Seurat3_functions.R")
path <- "Yang/Figure 4S/Supplementary Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data
(load(file="data/MCL_V3_Harmony_43_20190627.Rda"))

df_samples <- readxl::read_excel("doc/191030_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:12)))
df_samples = df_samples[sample_n,]

object@meta.data$orig.ident %<>% plyr::mapvalues(from = unique(df_samples$sample),
                                                 to = unique(df_samples$publication.id))
table(object@meta.data$orig.ident)
NewNames = paste0(object@meta.data$orig.ident,"_",object@meta.data$Barcode)
object %<>% RenameCells(new.names = NewNames)
rownames(object@reductions$tsne@cell.embeddings) = colnames(object)

Idents(object) = "groups"
object %<>% subset(idents = c("AFT-03","AFT-04"),invert = T)
Idents(object) = "orig.ident"
(samples = df_samples$publication.id[df_samples$publication.id %in% object$orig.ident])
Idents(object) %<>% factor(levels = samples)
table(Idents(object))

# Extend Data
Idents(object) = "manual"
object %<>% sortIdent()

cell_Freq <- table(Idents(object)) %>% as.data.frame
cell_Freq$Percent <- prop.table(cell_Freq$Freq) %>% scales::percent()

cell_Freq = cell_Freq[order(cell_Freq$Var1),]
cell_Freq$col = ExtractMetaColor(object)
cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")
cell_Freq$Cell_Type %<>% gsub("_"," ",.)

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=6, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,x.text.angle = 45,
          ylab = "Cell Number",
          label = "Percent",
          sort.val = "desc",
          width = 1, size = 0.5,
          title = "Numbers of major cell types in total 43 samples")+NoLegend()+
        theme(plot.title = element_text(hjust = 0.5,size=15))
dev.off()

#=======
QC_list <- read.csv(paste0(path,"QC_list.csv"), stringsAsFactors = F)
jpeg(paste0(path,"mean.Reads.per.Cell.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(QC_list, x = "submitter", y= "mean.Reads.per.Cell",
         title = "Mean reads per cell in each scRNA-seq",
         xlab = "",ylab = "Mean Reads per Cell",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5,
         yscale = "log10")+
                 theme(plot.title = element_text(hjust = 0.5,size=15),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank())
dev.off()


jpeg(paste0(path,"UMI.per.Cell.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(QC_list, x = "submitter", y= "median.Transcripts.per.Cell",
         title = "Median transcripts per cell in each scRNA-seq",
         xlab = "",ylab = "Median UMI per Cell",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5,
         yscale = "log10")+
        theme(plot.title = element_text(hjust = 0.5,size=13),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()

jpeg(paste0(path,"cell.number.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(QC_list, x = "submitter", y= "cell.number",
         title = "Cell number in each scRNA-seq",
         xlab = "",ylab = "Cell Number",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5)+
        theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()

# ================
Idents(object) = "orig.ident"
(mito.features <- grep(pattern = "^MT-", x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = "^MT-")
g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=8),legend.position="none")
})

jpeg(paste0(path,"S2_nGene.jpeg"), units="in", width=7, height=5,res=600)
g2[[1]]+ggtitle("Distribution of Gene Number per Cells")+
        xlab("scRNA-seq samples")+
        ylab("Gene Number")+
        ylim(0,max(object$nGene)+100)+
        theme(plot.title = element_text(face = 'plain'))
dev.off()

jpeg(paste0(path,"S2_nUMI.jpeg"), units="in", width=7, height=5,res=600)
g2[[2]]+ggtitle("Distribution of mRNA Counts per Cells")+
        xlab("scRNA-seq samples")+
        ylab("mRNA Counts")+
        ylim(0,max(object$nUMI)+1000)+
        theme(plot.title = element_text(face = 'plain'))
dev.off()

jpeg(paste0(path,"S2_mito.jpeg"), units="in", width=7, height=5,res=600)
g2[[3]]+ggtitle("Distribution of mitochondrial gene percentage per Cells")+
        xlab("scRNA-seq samples")+
        ylab("Mitochondrial gene percentage %")+
        theme(plot.title = element_text(face = 'plain'))
dev.off()
