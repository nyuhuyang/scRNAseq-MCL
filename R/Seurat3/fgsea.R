########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
library(ggplot2)
library(fgsea)
library(tibble)
library(ggpubr)
library(ggsci)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 3.1.1 load data
# Load some results from Seurat
res = read.csv(file="output/20190621/X5_clusters_FC0.1_markers.csv",
                        row.names = 1, stringsAsFactors=F)
table(res$cluster)
head(res)
res = res[order(res["p_val_adj"]),]
head(res, 20)
(clusters <- unique(res$cluster))
hallmark <- gmtPathways("../seurat_resources/msigdb/h.all.v6.2.symbols.gmt")
biocarta <- gmtPathways("../seurat_resources/msigdb/c2.cp.biocarta.v6.2.symbols.gmt")
kegg <- gmtPathways("../seurat_resources/msigdb/c2.cp.kegg.v6.2.symbols.gmt")
tft <- gmtPathways("../seurat_resources/msigdb/c3.tft.v6.2.symbols.gmt")
c6 <- gmtPathways("../seurat_resources/msigdb/c6.all.v6.2.symbols.gmt")
GO <- gmtPathways("../seurat_resources/msigdb/c5.all.v6.2.symbols.gmt")
allpathways <- c(hallmark,biocarta,kegg,c6)

hallmark %>% head() %>% lapply(head)
biocarta %>% head() %>% lapply(head)

# Now, run the fgsea algorithm with 1000 permutations:

for(i in 1:length(clusters)) FgseaBarplot(pathways=hallmark, stats=res, nperm=1000,
                                          cluster = i,
                                          sample="B_MCL cluster",pathway.name = "Hallmark", hjust=0.5)
for(i in 1:length(clusters)) FgseaBarplot(pathways=allpathways, stats=res, nperm=1000,show=50,
                                       cluster = clusters[i],sample = cell.line,hjust=0,
                                       pathway.name = "Hallmark, biocarta,and KEGG")

FgseaDotPlot(stats=res, pathways=hallmark, nperm=1000,padj = 0.25,pval = 0.05,
             order.by = c(4,"NES"),decreasing = F,
             size = "-log10(pval)", fill = "NES",sample = "each B_MCL clusters", 
             pathway.name = "Hallmark",rotate.x.text = F)

FgseaDotPlot(stats=res, pathways=allpathways, nperm=1000,padj = 0.1,pval = 0.02,
             order.by = c(4,"NES"),decreasing = F,
             size = "-log10(pval)", fill = "NES",sample = "each B_MCL clusters", 
             rotate.x.text = F, pathway.name = "Hallmark, biocarta,and KEGG")

df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
groups = c("Untreated","Pt-17","Pt-25")
for(i in 1:length(groups)){
        res_B = read.csv(file = paste0("output/20190622/B/B_MCL_DE/B_",groups[i],".csv"))
        res_T = read.csv(file = paste0("output/20190622/T/T_NK_DE/T_",groups[i],".csv"))
        
        (samples = df_samples$sample[df_samples$sample %in% unique(res_B$cluster)])
        #res_B$cluster %<>% factor(levels = samples)
        #res_T$cluster %<>% factor(levels = samples)
        
        FgseaDotPlot(stats=res_B, pathways=hallmark, nperm=1000,padj = 0.25,pval = 0.05,
                     order.by = c(4,"NES"),decreasing = F,
                     size = "-log10(pval)", fill = "NES",
                     sample = paste(groups[i],"B_MCL clusters"), 
                     pathway.name = "Hallmark",rotate.x.text = F)
        (samples = df_samples$sample[df_samples$sample %in% unique(res_T$cluster)])
        FgseaDotPlot(stats=res_T, pathways=hallmark, nperm=1000,padj = 0.25,pval = 0.05,
                     order.by = c(4,"NES"),decreasing = F,
                     size = "-log10(pval)", fill = "NES",
                     sample = paste(groups[i],"T_NK clusters"), 
                     pathway.name = "Hallmark",rotate.x.text = F)

}