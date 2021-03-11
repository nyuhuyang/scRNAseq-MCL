# ######################################################################
# conda activate r3.6.2
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
# load data
(load(file="data/MCL_41_harmony_20200225.Rda"))
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "cell.types"

# select pair
opts = data.frame(cell.types = rep(c("T_cells:CD8+","T_cells:CD4+","NK_cells","Monocytes"),  each = 24),
                  ident.1 = rep(c("Pt17_7","Pt17_12","Pt17_31","Pt17_31","Pt28_4","Pt28_28",
                                  "Pt25_1_8","Pt25_24","Pt25_25Pd","Pt25_25Pd","Pt11_28",
                                  "Pt25_1","PtU01","PtU02","PtU03","PtU04",
                                  "Pt17_LN1","Pt17_2","Pt17_7","Pt17_12","Pt17_31",
                                  "PtB13_Ibp","PtB13_Ib1","PtB13_IbR"),4),
                  ident.2 = rep(c("Pt17_2","Pt17_2","Pt17_2","Pt17_7","Pt28_1","Pt28_1",
                                  "Pt25_1","Pt25_1","Pt25_1","Pt25_24","Pt11_14",
                                  rep("N01",13)),4),
                  stringsAsFactors = F)
read.path = "Yang/Figure Sources/"
#for(i in c(11:23,34:46)){ #N01
#index = grep("_",opts$ident.2)
index = c(3:4,18,8:10,12,c(3:4,18,8:10,12)+24)
for(i in index){
        print(opt <- opts[i,])
        DE_res = paste0(read.path,sub(":","_",opt$cell.types),"/DE_files/",
                        sub(":","_",opt$cell.types),opt$ident.1, "-",opt$ident.2,".csv")

        sub_object <- subset(object, idents = opt$cell.types)
        Idents(sub_object) = "orig.ident"
        sub_object %<>% subset(idents = c(opt$ident.1, opt$ident.2))
        save.path = paste0(read.path,sub(":","_",opt$cell.types),"/",
                           sub(":","_",opt$cell.types), opt$ident.1, "-",opt$ident.2,"/")
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

        if(!file.exists(DE_res)) {
                print(paste(i, "file don't exists!"))
                markers <- FindAllMarkers.UMI(sub_object,
                                              logfc.threshold = 0,
                                              only.pos = F,
                                              return.thresh = 1,
                                              test.use = "MAST",
                                              latent.vars = "nFeature_SCT")
                write.csv(x = markers,file = DE_res)
        }

        deg = read.csv(file= DE_res, row.names = 1, stringsAsFactors=F)
        table(deg$cluster)
        markers <- FilterGenes(sub_object,c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B","GZMK"))
        (MT_gene <- grep("^MT-",deg$gene))
        if(length(MT_gene) >0 ) deg = deg[-MT_gene,]
        Top_n = 40

        top = deg %>% group_by(cluster) %>%
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
                         width=2.5, height=10,res=600,no.legend = F,
                         cex.row=5,
                         group.colors = c("#E6AB02","#2055da"),
                         colors = sns.RdBu_r_199,
                         do.print = T,
                         file.name = paste0("Heatmap_top40_",sub(":","_",opt$cell.types),
                                            opt$ident.1, "-",opt$ident.2,"_raw.jpeg"),
                         save.path = save.path)
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
                         width=2.5, height=10,res=600,no.legend = F,
                         cex.row=5,
                         group.colors = c("#E6AB02","#2055da"),
                         colors = sns.RdBu_r_199,
                         do.print = T,
                         file.name = paste0("Heatmap_top40_",sub(":","_",opt$cell.types),
                                            opt$ident.1, "-",opt$ident.2,"_scale.jpeg"),
                         save.path = save.path)
        deg = deg[deg$cluster %in% opt$ident.1,]
        #avg_logFC = deg[deg$cluster %in% opt$ident.2,"avg_logFC"]
        #deg[deg$cluster %in% opt$ident.2,"avg_logFC"] = avg_logFC * -1
        p <- VolcanoPlots(data = deg, cut_off_value = 0.05, cut_off = "p_val",
                          sort.by = "avg_logFC", cut_off_logFC = 0.1,
                          top = 20, cols = c("#2a71b2","#d2dae2","#ba2832"),alpha=1, size=2,
                          legend.size = 12)+
                ggtitle(paste0(opt$ident.1, " \\ ",opt$ident.2, " in ", sub(":","_",opt$cell.types)))+
                theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"),
                      legend.position="bottom",
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))
        jpeg(paste0(save.path,"VolcanoPlots_",sub(":","_",opt$cell.types),
                    opt$ident.1, "-",opt$ident.2,".jpeg"), units="in", width=10, height=7,res=600)
        print(p)
        dev.off()

        p1 <- VolcanoPlots(data = deg, cut_off_value = 0.05, cut_off = "p_val",
                          sort.by = "avg_logFC", cut_off_logFC = 0.1,
                          top = 0, cols = c("#2a71b2","#d2dae2","#ba2832"),alpha=1, size=2,
                          legend.size = 12)+
                ggtitle(paste0(opt$ident.1, " \\ ",opt$ident.2, " in ", sub(":","_",opt$cell.types)))+
                theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"),
                      legend.position="bottom",
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
        jpeg(paste0(save.path,"VolcanoPlots_",sub(":","_",opt$cell.types),
                    opt$ident.1, "-",opt$ident.2,"_noGene.jpeg"), units="in", width=10, height=7,res=600)
        print(p1)
        dev.off()

        Progress(i, length(index))
 }

# relocate DE_files
read.path = "Yang/Figure Sources/"
cell.types <- c("T_cells_CD8+","T_cells_CD4+")
for(cell.type in cell.types){

        csv_list <- list.files(paste0(read.path,cell.type),pattern = cell.type,
                               all.files = T,full.names = T,recursive = T)
        csv_list = grep("Top40_T_cells_|Heatmap_|VolcanoPlots", csv_list, value=T , invert = T)
        save.path = paste0("Yang/Figure Sources/",cell.type,"/DE_files/")
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

        for(csv in csv_list){
                file.rename(from = csv,
                            to = paste0(save.path, basename(csv)))
        }
}


