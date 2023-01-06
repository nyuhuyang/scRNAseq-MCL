invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

#======1.2 load  Seurat =========================

object <- readRDS(file = "data/MCL_SCT_87_20220901.rds")
meta.data <- readRDS(file = "output/MCL_SCT_87_20220901_meta.data_v3.rds")
if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}
cols= ExtractMetaColor(object, group.by = "X6cluster_AIM74")


test_df = data.frame(min_dist = rep(c(0.2,0.4,0.6),each = 3),
                     spread = rep(c(0.6,1.0,1.4),times = 3),
                     npcs = rep(c(50,60,70,80,90,100),each = 9))
test_df = bind_rows(list(test_df,test_df))
test_df$nfeatures = rep(c(2000,3000),each =54)

for(args in 1:108){
    nfeatures <- test_df[args,"nfeatures"]
    spread <- test_df[args,"spread"]
    min.dist <- test_df[args,"min_dist"]
    npcs <- test_df[args,"npcs"]
    file.name = paste0("npcs",npcs,"_dist.",min.dist,"_spread.",spread)

    read.path <- paste0("output/20220915/",nfeatures)
    object@reductions <- readRDS(file = paste0(read.path, "/reductions_",file.name,".rds"))

    save.path <- paste0("output/20220916/",nfeatures)
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

    sub_object <- subset(object, X6cluster_AIM74 != "MCL")

    plot <- UMAPPlot(object, group.by = "X6cluster_AIM74",raster=FALSE,label = T,
                     cols= cols)+ ggtitle(file.name)+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))+
        NoLegend()
    plot1 <- UMAPPlot(sub_object, group.by = "X6cluster_AIM74",raster=FALSE,label = T,
                      cols= cols)+ ggtitle(file.name)+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))+
        NoLegend()
    jpeg(paste0(save.path,"/", file.name,".jpeg"),units="in", width=10, height=7,res=600)
    print(plot)
    dev.off()
    jpeg(paste0(save.path,"/", file.name,"_noMCL.jpeg"),units="in", width=10, height=7,res=600)
    print(plot1)
    dev.off()
    Progress(args, 108)
}


cols= ExtractMetaColor(object, group.by = "X6cluster_MCL61")

for(args in 1:108){
    nfeatures <- test_df[args,"nfeatures"]
    spread <- test_df[args,"spread"]
    min.dist <- test_df[args,"min_dist"]
    npcs <- test_df[args,"npcs"]
    file.name = paste0("npcs",npcs,"_dist.",min.dist,"_spread.",spread)

    read.path <- paste0("output/20220915/",nfeatures)
    object@reductions <- readRDS(file = paste0(read.path, "/reductions_",file.name,".rds"))

    save.path <- paste0("output/20220916/",nfeatures)
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

    sub_object <- subset(object, X6cluster_MCL61 != "MCL")

    plot2 <- UMAPPlot(sub_object, group.by = "X6cluster_MCL61",raster=FALSE,label = T,
                      cols= cols)+ ggtitle(file.name)+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))+
        NoLegend()
    jpeg(paste0(save.path,"/", file.name,"_noMCL_MCL61.jpeg"),units="in", width=10, height=7,res=600)
    print(plot2)
    dev.off()
    Progress(args, 108)
}

#=============  nfeautre =3000,  ==================
object <- readRDS(file = "data/MCL_SCT_87_20220901.rds")
meta.data = readRDS(file = "output/MCL_B_only_SCT_87_20220901_meta.data_v2.rds")

for(args in 1:108){
    nfeatures <- test_df[args,"nfeatures"]
    spread <- test_df[args,"spread"]
    min.dist <- test_df[args,"min_dist"]
    npcs <- test_df[args,"npcs"]
    file.name = paste0("npcs",npcs,"_dist.",min.dist,"_spread.",spread)

    read.path <- paste0("output/20220915/",nfeatures)
    object@reductions <- readRDS(file = paste0(read.path, "/reductions_",file.name,".rds"))

    save.path <- paste0("output/20221007/",nfeatures)
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

    sub_object <- subset(object, cells = rownames(meta.data))
    sub_object@meta.data = meta.data
    sub_object <- subset(sub_object, subset = SCT_snn_res.0.09 %in% 0:4)

    plot1 <- UMAPPlot(sub_object, group.by = "SCT_snn_res.0.09",raster=FALSE,label = T)+
        ggtitle(file.name)+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))+
        NoLegend()
    jpeg(paste0(save.path,"/", file.name,"_res0.09_0_1_2_3_4.jpeg"),units="in", width=10, height=7,res=600)
    print(plot1)
    dev.off()
    Progress(args, 108)
}


#=================
cols= ExtractMetaColor(object, group.by = "celltype.l1")

for(args in 1:108){
    nfeatures <- test_df[args,"nfeatures"]
    spread <- test_df[args,"spread"]
    min.dist <- test_df[args,"min_dist"]
    npcs <- test_df[args,"npcs"]
    file.name = paste0("npcs",npcs,"_dist.",min.dist,"_spread.",spread)

    read.path <- paste0("output/20220915/",nfeatures)
    save.path <- paste0("output/20221007/",nfeatures)
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

    object@reductions <- readRDS(file = paste0(read.path, "/reductions_",file.name,".rds"))

    plot1 <- UMAPPlot(object, group.by = "celltype.l1",raster=FALSE,label = T,cols= cols)+
        ggtitle(file.name)+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))+
        NoLegend()
    jpeg(paste0(save.path,"/", file.name,"_celltype.l1.jpeg"),units="in", width=10, height=7,res=600)
    print(plot1)
    dev.off()
    Progress(args, 108)
}

#=================
cols= ExtractMetaColor(object, group.by = "celltype.l1")

for(args in 1:108){
    nfeatures <- test_df[args,"nfeatures"]
    spread <- test_df[args,"spread"]
    min.dist <- test_df[args,"min_dist"]
    npcs <- test_df[args,"npcs"]
    file.name = paste0("npcs",npcs,"_dist.",min.dist,"_spread.",spread)

    read.path <- paste0("output/20220915/",nfeatures)
    save.path <- paste0("output/20221007/",nfeatures)
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

    object@reductions <- readRDS(file = paste0(read.path, "/reductions_",file.name,".rds"))

    plot1 <- UMAPPlot(object, group.by = "celltype.l1",raster=FALSE,label = T,cols= cols)+
        ggtitle(file.name)+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))+
        NoLegend()
    jpeg(paste0(save.path,"/", file.name,"_celltype.l1.jpeg"),units="in", width=10, height=7,res=600)
    print(plot1)
    dev.off()
    Progress(args, 108)
}
