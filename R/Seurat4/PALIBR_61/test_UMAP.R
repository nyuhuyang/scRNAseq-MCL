# test
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)


object = readRDS(file = "data/MCL_61_20220331.rds")
object@meta.data = readRDS(file = "MCL_61_20220331_metadata.rds")

test_df = data.frame(min_dist = rep(c(0.2,0.4,0.6),each = 3),
                     spread = rep(c(0.6,1.0,1.4),times = 3),
                     npcs = rep(c(50,60,70,80,90,100),each = 9))
test_df = bind_rows(list(test_df,test_df,test_df))
test_df$sample = rep(c("all", "AFT12.Pt22_64", "Pt22_64"),each =54)

file.names = paste0("cs",test_df$npcs,"_dist.",test_df$min_dist,"_spread.",test_df$spread)

#=========== all + AFT12.Pt22 ========================
save.path <- "output/20220401/all"
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)


for (i in 1:54){
    print(test_df[i,])
    file.name = file.names[i]
    if(!file.exists(paste0("output/20220331/umap_",file.name,".rds"))) next
    umap = readRDS(paste0("output/20220331/umap_",file.name,".rds"))
    object@reductions[[paste0("umap_",file.name)]] = umap
    object@reductions[["umap"]] =  umap[[1]]
    colnames(object@reductions[["umap"]]@cell.embeddings) %<>% gsub(".*UMAP","UMAP",.)
    object@reductions[["umap"]]@key %<>% gsub(".*UMAP","UMAP",.)

    UMAPPlot.1(object, group.by = "cell.types",do.print = T,raster=FALSE,no.legend = T,
           title = file.name,file.name = paste0(file.name,".jpeg"),
           save.path = save.path)
}

#=========== AFT12.Pt22_64 ========================
save.path <- "output/20220401/AFT12.Pt22_64"
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

object %<>% subset(subset = orig.ident %in% "Pt22_67", invert = T)
object$orig.ident  %<>% droplevels()

for (i in 55:108){
    print(test_df[i,])
    file.name = file.names[i-54]
    if(!file.exists(paste0("output/20220401/AFT12.Pt22_64/umap_",file.name,".rds"))) next
    umap = readRDS(paste0("output/20220401/AFT12.Pt22_64/umap_",file.name,".rds"))
    object@reductions[[paste0("umap_",file.name)]] = umap
    object@reductions[["umap"]] =  umap[[1]]
    colnames(object@reductions[["umap"]]@cell.embeddings) %<>% gsub(".*UMAP","UMAP",.)
    object@reductions[["umap"]]@key %<>% gsub(".*UMAP","UMAP",.)

    UMAPPlot.1(object, group.by = "cell.types",do.print = T,raster=FALSE,
               title = file.name,file.name = paste0(file.name,".jpeg"),
               save.path = save.path)
}


#=========== Pt22_64 ========================
save.path <- "output/20220401/Pt22_64"
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

object %<>% subset(subset = orig.ident %in% c(grep("AFT12",unique(object$orig.ident), value = T),
                                              "Pt22_67"), invert = T)
object$orig.ident  %<>% droplevels()

for (i in 109:162){
    print(test_df[i,])
    file.name = file.names[i-108]
    if(!file.exists(paste0("output/20220401/Pt22_64/umap_",file.name,".rds"))) next
    umap = readRDS(paste0("output/20220401/Pt22_64/umap_",file.name,".rds"))
    object@reductions[[paste0("umap_",file.name)]] = umap
    object@reductions[["umap"]] =  umap[[1]]
    colnames(object@reductions[["umap"]]@cell.embeddings) %<>% gsub(".*UMAP","UMAP",.)
    object@reductions[["umap"]]@key %<>% gsub(".*UMAP","UMAP",.)

    UMAPPlot.1(object, group.by = "cell.types",do.print = T,raster=FALSE,
               title = file.name,file.name = paste0(file.name,".jpeg"),
               save.path = save.path)
}
