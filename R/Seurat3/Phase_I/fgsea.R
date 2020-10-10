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
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- "Yang/Figure 2/Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)


choose <- c("X4clusters","X4cluster_vs_Normal")[1]
res = read.csv(file = paste0("Yang/Figure 2/Figure Sources/",
                             choose,"/",choose,"_41-FC0.csv"),
               row.names = 1, stringsAsFactors=F)
if(choose == "X4clusters") table(res$cluster)
if(choose == "X4cluster_vs_Normal"){
        table(res$cluster1.vs.cluster2)
        res$cluster1.vs.cluster2 = gsub("Normal vs. ","",res$cluster1.vs.cluster2)
}

head(res)
res = res[order(res["p_val_adj"]),]
head(res, 20)
(clusters <- unique(res$cluster))
hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v6.2.symbols.gmt")
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(hallmark) = gsub("\\_"," ",names(hallmark))

# Now, run the fgsea algorithm with 1000 permutations:
set.seed(100)
fgseaRes = FgseaDotPlot(stats=res, pathways=hallmark,
                        padj = 0.25,pval = 0.05,
                        order.yaxis.by = c("C4","NES"),decreasing = F,
                        order.xaxis = paste0("C",1:4),
                        Rowv = F,Colv = F,
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = "Hallmark",rotate.x.text = T,
                        title = "in B and MCL sub-clusters",
                        font.xtickslab=12, font.main=12, font.ytickslab = 10,
                        font.legend = list(size = 12),font.label = list(size = 12),
                        do.return = T,save.path = path, do.print = T,
                        width = 5,height = 6,hjust = 0.75)
write.csv(fgseaRes, file = paste0(path,choose,"_FDR0.25_pval0.05.csv"))
