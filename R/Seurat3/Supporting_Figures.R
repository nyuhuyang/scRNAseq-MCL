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
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data
(load(file="data/MCL_V3_Harmony_43_20190627.Rda"))
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

jpeg("Yang/Figure 2/Supplementary/cell_type_numbers.jpeg", units="in", width=6, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,x.text.angle = 45,
          label = "Percent",
          sort.val = "desc",
          width = 1, size = 0.5,
          title = "Numbers of major cell types in total 43 samples")+NoLegend()+
        theme(plot.title = element_text(hjust = 0.5,size=15))
dev.off()

#=======
QC_list <- read.csv("Yang/Figure 2/Supplementary/QC_list.csv", stringsAsFactors = F)
jpeg("Yang/Figure 2/Supplementary/mean.Reads.per.Cell.jpeg", units="in", width=5, height=5,res=600)
ggviolin(QC_list, x = "submitter", y= "mean.Reads.per.Cell",
         title = "Mean reads per cell in each scRNA-seq",
         xlab = "",ylab = "mean Reads per Cell",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5,
         yscale = "log10")+
                 theme(plot.title = element_text(hjust = 0.5,size=15),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank())
dev.off()


jpeg("Yang/Figure 2/Supplementary/UMI.per.Cell.jpeg", units="in", width=5, height=5,res=600)
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

jpeg("Yang/Figure 2/Supplementary/cell.number.jpeg", units="in", width=5, height=5,res=600)
ggviolin(QC_list, x = "submitter", y= "cell.number",
         title = "Cell number in each scRNA-seq",
         xlab = "",ylab = "cell number",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5)+
        theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()
