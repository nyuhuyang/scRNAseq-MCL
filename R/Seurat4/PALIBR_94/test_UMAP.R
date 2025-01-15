invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

#======1.2 load  Seurat =========================

object <- readRDS(file = "data/MCL_SCT_94_20240615.rds")
meta.data <- readRDS(file = "output/MCL_SCT_94_20240615_meta.data_v2.rds")
if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}

test_df = data.frame(min_dist = rep(c(0.2,0.4,0.6),each = 3),
                     spread = rep(c(0.6,1.0,1.4),times = 3),
                     npcs = rep(c(50,60,70,80,90,100),each = 9))
test_df = bind_rows(list(test_df,test_df))
test_df$nfeatures = rep(c(2000,3000),each =54)


object$X6cluster_AIM74 %<>% factor(levels = sort(unique(object$X6cluster_AIM74)))
cols_AIM74= ExtractMetaColor(object, group.by = "X6cluster_AIM74")
object$X6cluster_MCL61 %<>% factor(levels = sort(unique(object$X6cluster_MCL61)))
cols_MCL61= ExtractMetaColor(object, group.by = "X6cluster_MCL61")


for(args in 1:108){
    nfeatures <- test_df[args,"nfeatures"]
    spread <- test_df[args,"spread"]
    min.dist <- test_df[args,"min_dist"]
    npcs <- test_df[args,"npcs"]
    file.name = paste0("npcs",npcs,"_dist.",min.dist,"_spread.",spread)

    read.path <- paste0("output/20240917/",nfeatures)
    if(! file.exists(paste0(read.path, "/reductions_",file.name,".rds"))) {
        print(paste(args,"don't exist"))
        next
    } else print(paste(args,"exist"))
}


for(args in 1:65){
    print(args)
    nfeatures <- test_df[args,"nfeatures"]
    spread <- test_df[args,"spread"]
    min.dist <- test_df[args,"min_dist"]
    npcs <- test_df[args,"npcs"]
    file.name = paste0("npcs",npcs,"_dist.",min.dist,"_spread.",spread)

    read.path <- paste0("output/20240917/",nfeatures)
    if(! file.exists(paste0(read.path, "/reductions_",file.name,".rds"))) {
        print(paste(args,"don't exist"))
        next
    } else print(paste(args,"exist"))
    object@reductions <- readRDS(file = paste0(read.path, "/reductions_",file.name,".rds"))
    sub_object_AIM74 <- subset(object, X6cluster_AIM74 != "MCL")

    save.path <- paste0("output/20241018/",nfeatures)
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)


    plot <- UMAPPlot(object, group.by = "X6cluster_AIM74",raster=FALSE,label = T,
                     cols= cols_AIM74)+ ggtitle(paste(file.name,"copy AIM74"))+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))+
        NoLegend()
    plot1 <- UMAPPlot(sub_object_AIM74, group.by = "X6cluster_AIM74",raster=FALSE,label = T,
                      cols= cols_AIM74)+ ggtitle(paste(file.name,"copy AIM74, rm MCL"))+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))+
        NoLegend()
    jpeg(paste0(save.path,"/", file.name,"_AIM74.jpeg"),units="in", width=10, height=7,res=600)
    print(plot)
    dev.off()
    jpeg(paste0(save.path,"/", file.name,"_AIM74_noMCL.jpeg"),units="in", width=10, height=7,res=600)
    print(plot1)
    dev.off()
    sub_object_MCL61 <- subset(object, X6cluster_MCL61 != "MCL")

    plot2 <- UMAPPlot(object, group.by = "X6cluster_MCL61",raster=FALSE,label = T,
                      cols= cols_AIM74)+ ggtitle(paste(file.name,"copy MCL61"))+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))+
        NoLegend()

    plot3 <- UMAPPlot(sub_object_MCL61, group.by = "X6cluster_MCL61",raster=FALSE,label = T,
                      cols= cols_MCL61)+ ggtitle(paste(file.name,"copy MCL61, rm MCL"))+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))+
        NoLegend()
    jpeg(paste0(save.path,"/", file.name,"_MCL61.jpeg"),units="in", width=10, height=7,res=600)
    print(plot2)
    dev.off()
    jpeg(paste0(save.path,"/", file.name,"_MCL61_noMCL.jpeg"),units="in", width=10, height=7,res=600)
    print(plot3)
    dev.off()

    plot4 <- FeaturePlot(object, features = c("PCNA","E2F1","BCL2A1","BCL2","IRF4","MKI67"),raster=FALSE)
    jpeg(paste0(save.path,"/", file.name,"_features.jpeg"),units="in", width=7, height=9,res=600)
    print(plot4)
    dev.off()
    Progress(args, 65)
}
