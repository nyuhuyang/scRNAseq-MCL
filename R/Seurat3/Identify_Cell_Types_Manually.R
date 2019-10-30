library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
marker_path <- paste0(path,"markers/")
if(!dir.exists(marker_path))dir.create(marker_path, recursive = T)

# ======== 2.1 =========== test with known markers==================
(load(file="data/MCL_V3_Harmony_43_20190627.Rda"))
DefaultAssay(object) <- "RNA"
#df_markers <- readxl::read_excel("doc/Lynch mouse model scRNAseq genes of interest 073119.xlsx")

df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx",sheet = "Human.sub")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
markers = df_markers[,grep("Alias",colnames(df_markers),invert = T)]
marker.list <- df2list(markers)

marker.list %<>% lapply(function(x) x[1:18]) %>% 
     lapply(function(x) FilterGenes(object,x)) %>% 
     lapply(function(x) x[!is.na(x)]) %>% 
    lapply(function(x) x[1:min(length(x),12)])
marker.list <- marker.list[!is.na(marker.list)]
marker.list %>% list2df %>% t %>% kable() %>% kable_styling()

object@meta.data = object@meta.data[,grep("ident",colnames(object@meta.data),invert = T)]
Idents(object) = "integrated_snn_res.0.6"
for(i in 1:length(marker.list)){
    if(length(marker.list[[i]]) == 0) next
    p <- lapply(marker.list[[i]], function(marker) {
        FeaturePlot(object = object, feature = marker,pt.size = 0.5,
                    reduction="tsne", label = T)+
            NoLegend()+
            ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+
            theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"))
    })
    jpeg(paste0(path,"markers/",names(marker.list)[i],".jpeg"),units="in", width=10, height=7,res=600)
    print(do.call(cowplot::plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
              theme(plot.title = element_text(hjust = 0.5,size = 20)))
    dev.off()
    print(paste0(i,":",length(marker.list)))
}

#=============== multiple color in the single Featureplot===================================
(load(file = "data/B_cells_MCL_43_20190917.Rda"))
Idents(B_cells_MCL) = "orig.ident"
samples = c("Pt-10-LN-C2","Pt-11-LN-C1", "Pt-17-LN-C1")
features.list = lapply(list(c("EZH2","E2F1"),
                            c("EZH2","PCNA"),
                            c("EZH2","CDK1"),
                            c("EZH2","EZH1")),
                       function(x) FilterGenes(B_cells_MCL,x,unique = F))
cols.use.list = list(c("#b88801","#2c568c", "#E31A1C"),
                     c("#b88801","#2c568c", "#E31A1C"),
                     c("#b88801","#2c568c", "#E31A1C"),
                     c("#b88801","#2c568c", "#E31A1C"))
for(s in samples){
    s_path <- paste0(path,s,"/")
    if(!dir.exists(s_path)) dir.create(s_path, recursive = T)
    subset_object = subset(B_cells_MCL, idents = s)
    for(i in 1:length(features.list)){
        # FeaturePlot.2
        g <- FeaturePlot.2(object = subset_object, features = features.list[[i]],do.return = T,
                           overlay = T,cols = c("#d8d8d8",cols.use.list[[i]]),
                           pt.size = c(1,3), alpha = c(0.5, 1))
        jpeg(paste0(s_path,s,"_",paste(features.list[[i]],collapse = "_"),".jpeg"), 
             units="in", width=7, height=7,res=600)
        print(g+theme(plot.title = element_text(hjust = 0.5,size = 20),
                      legend.position="bottom",
                      legend.text=element_text(size=15))+
                  guides(colour = guide_legend(override.aes = list(size=5)), 
                         shape = guide_legend(override.aes = list(size=5))))
        dev.off()
        # percentage
        Idents(subset_object) = "X5_clusters"
        subset_object_4 = subset(subset_object, idents = 4)
        
        features_var <- FetchData(subset_object_4,features.list[[i]])
        df <- table(features_var[,1]>0,features_var[,2]>0) %>% prop.table(margin = NULL) %>%
            as.data.frame() %>% spread(Var2,Freq)
        df$Var1 = plyr::mapvalues(df$Var1, c(FALSE, TRUE), 
                                  paste(features.list[[i]][1],c("== 0","> 0")))
        colnames(df) = c(features.list[[i]][1],
                         paste(features.list[[i]][2],c("== 0","> 0")))
        write.csv(df,file = paste0(s_path,s,"_",paste(features.list[[i]],
                                                      collapse = "_"),".csv"))
        
    }
}
Idents(subset_object) = "orig.ident"
RidgePlot(subset_object, features = c("EZH1", "EZH2"), ncol = 1)
FeatureScatter(B_cells_MCL, feature1 = "EZH2",feature2 = "CDK1",slot = "data")+NoLegend()
