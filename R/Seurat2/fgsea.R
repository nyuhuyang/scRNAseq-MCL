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
source("../R/Seurat_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 3.1.1 load data
# Load some results from Seurat
res = read.csv("output/20190603/X5_clusters_FC0.01_markers.csv")
table(res$cluster)
head(res)
res = res[order(res["p_val_adj"]),]
head(res, 20)

hallmark <- gmtPathways("../seurat_resources/msigdb/h.all.v6.2.symbols.gmt")
biocarta <- gmtPathways("../seurat_resources/msigdb/c2.cp.biocarta.v6.2.symbols.gmt")
kegg <- gmtPathways("../seurat_resources/msigdb/c2.cp.kegg.v6.2.symbols.gmt")
tft <- gmtPathways("../seurat_resources/msigdb/c3.tft.v6.2.symbols.gmt")
c6 <- gmtPathways("../seurat_resources/msigdb/c6.all.v6.2.symbols.gmt")
GO <- gmtPathways("../seurat_resources/msigdb/c5.all.v6.2.symbols.gmt")
pathways <- c(hallmark,biocarta,kegg,c6)

hallmark %>% head() %>% lapply(head)
biocarta %>% head() %>% lapply(head)
# Now, run the fgsea algorithm with 1000 permutations:

FgseaPlot <- function(pathways=hallmark, stats=res, nperm=1000,show=NULL,cluster = 1,
                      pathway.name = "Hallmark"){
        
        res = stats[order(stats["p_val_adj"]),]
        res1 = res[res$cluster == cluster,c("gene","avg_logFC")] %>% deframe

        fgseaRes <- fgsea(pathways=pathways, stats=res1, nperm=nperm)
        fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
        fgseaResTidy %<>% 
                dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
                arrange(padj)
        print(dim(fgseaResTidy))
        if(!is.null(show)){
                (N = nrow(fgseaResTidy)-1)
                fgseaResTidy =fgseaResTidy[c(1:(show/2),(N-show/2):N),]
        }
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,cluster,"-",pathway.name,".jpeg"), units="in", width=10, height=7,res=600)
        g <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
                geom_col(aes(fill=pval<0.05)) +
                guides(fill = guide_legend(reverse = TRUE))+
                coord_flip() +
                labs(x="Pathway", y="Normalized Enrichment Score",
                     title=paste(pathway.name,"pathways in B/MCL cluster",cluster)) + 
                theme_minimal()
        print(g)
        dev.off()
}

 

for(i in 1:5) FgseaPlot(pathways=hallmark, stats=res, nperm=1000,show=NULL,cluster = i,
                        pathway.name = "Hallmark")
for(i in 1:5) FgseaPlot(pathways=pathways, stats=res, nperm=1000,show=50,cluster = i,
                        pathway.name = "Hallmark, biocarta, KEGG, and oncogenic")
DotPlot <- function(){
                plot<- ggplot(data = data.plot, mapping = aes_string(x = "features.plot",
                                                                 y = "id")) + 
                geom_point(mapping = aes_string(size = "pct.exp",color = color.by)) + 
                scale.func(range = c(0, dot.scale),limits = c(scale.min, scale.max)) + 
                theme(axis.title.x = element_blank(),
                      labs(x = "Features", y = ifelse(test = is.null(x = split.by),
                                                      yes = "Identity", no = "Split Identity")) + 
                              theme_cowplot()
}