########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#devtools::install_github(repo = "ChristophH/sctransform", ref = "develop")
library(Seurat)
library(dplyr)
library(kableExtra)
library(magrittr)
library(ggplot2)
library(cowplot)
library(fgsea)
library(tibble)
library(ggsci)
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data
(load(file="data/MCL_V3_Harmony_43_20190627.Rda"))
# preprocess
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "res.0.6"
object %<>% subset(idents = 11,invert = T)
Idents(object) = "manual"
object <- subset(object,idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)
object %<>% sortIdent()
table(Idents(object))

#==== Figure 2-A ===========
TSNEPlot.1(object = object, label = T, label.repel = T, group.by = "manual",
           cols = ExtractMetaColor(object),no.legend = T,border = T,
           pt.size = 0.1,label.size = 5, do.print = T,do.return = F,
           title = "Cell type labeling by Blueprint + Encode + MCL",title.size = 15,
           units= "cm",width=14, height=14)

#==== Figure 2-B ===========
features <- FilterGenes(object,c("CD19","CCND1","SOX11",
                                 "CD3D","CD4","CD8A",
                                 "GNLY","KLRC1","NCAM1",
                                 "MS4A7","CD14","FCGR1A"))
FeaturePlot.1(object,features = features, pt.size = 0.5, cols = c("gray90", "red"), alpha = 1,
              threshold = 1, strip.text.size = 30, border = T,do.print = T, do.return = F,ncol = 3, 
              units = "in",width=9, height=12)

#==== Figure 2-C ===========
table(Idents(object))
object@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",object@meta.data$orig.ident)
Idents(object) = "orig.ident"

df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
tests <- paste0("test",c(2))
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df <- as.data.frame(df_samples[sample_n,])
        samples <- unique(df$sample)
        rownames(df) = samples
        
        samples <- c(ifelse(length(samples)>5,NA,"Normal"),df$sample[order(df$tsne)])
        print(samples <- samples[!is.na(samples)])
        
        subset_object <- subset(object, idents = samples)
        subset_object$orig.ident %<>% factor(levels = samples)
        
        Idents(subset_object) = "manual"
        
        subset_object %<>% sortIdent()
        TSNEPlot.1(subset_object, pt.size =0.3, 
                   strip.text.size = min(240/max(stringr::str_length(samples)),30),
                   group.by = "manual",split.by = "orig.ident",
                   cols = ExtractMetaColor(subset_object), ncol = length(samples),
                   unique.name = T, do.print = T,do.return = F,border = T,
                   width=length(samples)*2+2, height=3)
}

#==== Figure 2-D ===========
(load(file = "data/B_cells_MCL_43_20190713.Rda"))
TSNEPlot.1(B_cells_MCL, pt.size = 0.2,label = T, label.repel = T,
           do.print = T,no.legend = T,border = T,alpha = 1,
           label.size = 5, repel = T, title = "5 clusters in B/MCL cells",
           width=7, height=7)

#==== Figure 2-E ===========
X5_clusters_markers = read.csv("output/20190717/X5_clusters_FC0.2_markers.csv",row.names = 1)

markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
(MT_gene <- grep("^MT-",X5_clusters_markers$gene))
X5_clusters_markers = X5_clusters_markers[-MT_gene,]
Top_n = 40
top = X5_clusters_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
B_cells_MCL %<>% ScaleData(features=unique(c(as.character(top$gene),markers)))
DoHeatmap.1(B_cells_MCL, marker_df = X5_clusters_markers, features = markers, Top_n = Top_n,
            do.print=T, angle = 0, group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0, label=F, cex.row= 2, legend.size = 0,width=10, height=6.5,
            pal_gsea = FALSE,
            title = paste("Top",Top_n,"differentially expressed genes in B and MCL clusters"))

#==== Figure 2-F ===========
res = read.csv(file="output/20190621/X5_clusters_FC0.1_markers.csv",
               row.names = 1, stringsAsFactors=F)
table(res$cluster)
head(res)
res = res[order(res["p_val_adj"]),]
head(res, 20)
(clusters <- unique(res$cluster))
hallmark <- gmtPathways("../seurat_resources/msigdb/h.all.v6.2.symbols.gmt")
FgseaDotPlot(stats=res, pathways=hallmark, nperm=1000,padj = 0.25,pval = 0.05,
             order.by = c(4,"NES"),decreasing = F,
             size = "-log10(pval)", fill = "NES",sample = "each B_MCL clusters", 
             pathway.name = "Hallmark",rotate.x.text = F)
