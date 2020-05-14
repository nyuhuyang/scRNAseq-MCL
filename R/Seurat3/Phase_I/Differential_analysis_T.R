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
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# load data
(load(file="data/MCL_41_harmony_20200225.Rda"))
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "cell.types"

# select pair
opts = data.frame(cell.types = rep(c("T_cells:CD8+","T_cells:CD4+"),  each = 15),
                  ident.1 = rep(c("Pt17_7","Pt17_12","Pt17_31","Pt28_4","Pt28_28","Pt25_1",
                                  "Pt25_1_8","Pt25_24","Pt25_25Pd","Pt25_25Pd","Pt11_28",
                                  "PtU01","PtU02","PtU03","PtU04"),2),
                  ident.2 = rep(c("Pt17_2","Pt17_2","Pt17_2","Pt28_1","Pt28_1","N01",
                                  "Pt25_1","Pt25_1","Pt25_1","Pt25_24","Pt11_14",
                                  "N01", "N01", "N01","N01"),2),
                  stringsAsFactors = F)
for(i in c(12:15,27,28,30)){
        (opt = opts[i,])
        sub_object <- subset(object, idents = opt$cell.types)
        Idents(sub_object) = "orig.ident"
        sub_object %<>% subset(idents = c(opt$ident.1,opt$ident.2))
        save.path = paste0(path,sub(":","_",opt$cell.types), opt$ident.1, "-",opt$ident.2,"/")
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        file.copy(from=paste0(path,sub(":","_",opt$cell.types), opt$ident.1, "-",opt$ident.2,".csv"),
                  to=paste0(save.path,basename(save.path),".csv"))
        
        T_markers = read.csv(file= paste0(save.path,basename(save.path),".csv"),
                                      row.names = 1, stringsAsFactors=F)
        table(T_markers$cluster)
        markers <- FilterGenes(sub_object,c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B","GZMK"))
        (MT_gene <- grep("^MT-",T_markers$gene))
        if(length(MT_gene) >0 ) T_markers = T_markers[-MT_gene,]
        Top_n = 40

        top = T_markers %>% group_by(cluster) %>%
                top_n(40, avg_logFC)
        unique(top$cluster)
        #top = top[order(top$cluster),]
        write.csv(top,paste0(save.path,"Top40_",sub(":","_",opt$cell.types),
                             opt$ident.1, "-",opt$ident.2,".csv"))
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = sub_object), 2),
                     markers)
        #DoHeatmap.1======
        # raw heatmap
        featuresNum <- make.unique(features, sep = ".")
        exp = AverageExpression(sub_object[features,], 
                                assays = "SCT") %>% .$SCT
        exp %<>% MakeUniqueGenes(features = features)
        exp[tail(VariableFeatures(object = sub_object), 2),] =0

        (group.by = c(opt$ident.2, opt$ident.1))
        DoHeatmap.matrix(exp, features = featuresNum,
                         group.by = group.by,
                         size = 6,angle = 0,label =F,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         width=2.5, height=10,res=600,no.legend = F,
                         cex.row=5,
                         group.colors = gg_color_hue(2),
                         do.print = T,
                         unique.name = "cell.types",
                         save.path = paste0(save.path,"Heatmap_top40_",sub(":","_",opt$cell.types),
                                            opt$ident.1, "-",opt$ident.2,"_raw"))
        # scale heatmap
        sub_object %<>% ScaleData(features = features)
        scale_exp = AverageExpression(sub_object[features,], 
                                      assays = "SCT", slot = "scale.data") %>% .$SCT
        scale_exp %<>% MakeUniqueGenes(features = features)
        scale_exp[tail(VariableFeatures(object = sub_object), 2),] =0
        DoHeatmap.matrix(scale_exp, features = featuresNum,
                         group.by = group.by,
                         size = 6,angle = 0,label =F,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         width=2.5, height=10,res=600,no.legend = F,
                         cex.row=5,
                         group.colors = gg_color_hue(2),
                         do.print = T,
                         unique.name = "cell.types",
                         save.path = paste0(save.path,"Heatmap_top40_",sub(":","_",opt$cell.types),
                                            opt$ident.1, "-",opt$ident.2,"_scale"))
        #T_markers = T_markers[T_markers$cluster %in% opt$ident.1,]
        avg_logFC = T_markers[T_markers$cluster %in% opt$ident.2,"avg_logFC"]
        T_markers[T_markers$cluster %in% opt$ident.2,"avg_logFC"] = avg_logFC * -1
        p <- VolcanoPlots(data = T_markers, cut_off_value = 0.05, cut_off = "p_val", cut_off_logFC = 0.1,
                          top = 20, cols = c("#2a52be","#d2dae2","#d9321f"),alpha=1, size=2,
                          legend.size = 12)+
                ggtitle(paste0(opt$ident.1, " \\ ",opt$ident.2, " in ", sub(":","_",opt$cell.types)))+
                theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"),
                      legend.position="bottom")
        jpeg(paste0(save.path,"VolcanoPlots_",sub(":","_",opt$cell.types),
                    opt$ident.1, "-",opt$ident.2,".jpeg"), units="in", width=10, height=7,res=600)
        print(p)
        dev.off()
        Progress(i, nrow(opts))
 }
