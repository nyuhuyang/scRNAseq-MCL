invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

# read sample summary list
df_samples <- readxl::read_excel("output/20220318/20220318_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples %<>% filter(sequence %in% "GEX") %>% filter(phase %in% "PALIBR_I")
nrow(df_samples)
df_samples$date %<>% gsub(" UTC","",.) %>% as.character()
#======================================
object = readRDS(file = "data/MCL_61_20220331.rds")
meta.data = readRDS(file = "data/MCL_61_20220331_metadata.rds")
meta.data$patient = plyr::mapvalues(meta.data$orig.ident,
                                    from = df_samples$sample,
                                    to = df_samples$patient)
meta.data$response =  plyr::mapvalues(meta.data$orig.ident,
                                      from = df_samples$sample,
                                      to = df_samples$response)
meta.data$treatment =  plyr::mapvalues(meta.data$orig.ident,
                                       from = df_samples$sample,
                                       to = df_samples$treatment)
meta.data$Mean.Reads.per.Cell %<>% gsub(",","",.) %>% as.integer()
meta.data$Number.of.Reads %<>% gsub(",","",.) %>% as.integer()
meta.data$Sequencing.Saturation %<>% gsub("%","",.) %>% as.numeric()
meta.data$barcode  = rownames(meta.data)

# find X4cluster
meta_data = read.csv("shinyApp/PALIBR_I_51/cnv_meta_data.csv",row.names = 1)
meta_data$barcode  = rownames(meta_data)

meta.data %<>% left_join(meta_data[,c("X4cluster","barcode")], by = "barcode")
meta.data %<>% tibble::column_to_rownames(var = "barcode")

table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data
object$response %<>% factor(levels = c("Normal","Untreated","CR","PR","PD"))
object$treatment %<>% factor(levels = c("Normal","Untreated","PALIBR+Ibrutinib"))

DefaultAssay(object) = "SCT"
s.genes <- cc.genes$s.genes %>% gsub("MLF1IP","CENPU",.)
g2m.genes <- cc.genes$g2m.genes %>% plyr::mapvalues(from = c("FAM64A", "HN1"),
                                                    to = c("PIMREG","JPT1"))
object %<>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
colnames(object@meta.data) %<>% sub("Phase","cell cycle phase",.)



test_df = data.frame(min_dist = rep(c(0.2,0.4,0.6),each = 3),
                     spread = rep(c(0.6,1.0,1.4),times = 3),
                     npcs = rep(c(50,60,70,80,90,100),each = 9))
test_df = bind_rows(list(test_df,test_df))
test_df$nfeatures = rep(c(2000,3000),each =54)

for(i in c(71,72)){
    print(nfeatures <- test_df[i,"nfeatures"])
    print(spread <- test_df[i,"spread"])
    print(min.dist <- test_df[i,"min_dist"])
    print(npcs <- test_df[i,"npcs"])
    umap.name = paste0("nfeatures",nfeatures,"npcs",npcs,"dist.",min.dist,"spread.",spread) %>% gsub("\\.","",.)
    file.name = paste0("output/20220405/",nfeatures,"/umap_cs",npcs,"_dist.",min.dist,"_spread.",spread,".rds")
    umap  = readRDS(file.name)[[1]]
    umap@key = paste0(umap.name,"_")
    colnames(umap@cell.embeddings) = paste0(umap.name,"_",1:2)
    object[[paste0("umap_",umap.name)]] <- umap
    Progress(i,72)
}
object[["umap"]] = object[["umap_nfeatures3000npcs60dist06spread1"]]
object[["umap"]]@key = "UMAP_"
colnames(object[["umap"]]@cell.embeddings) = c("UMAP_1","UMAP_2")
#system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))

object %<>% FindNeighbors(reduction = "umap",dims = 1:2)
resolutions = c( 0.01, 0.1, 0.2, 0.5,0.8)
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    Progress(i,length(resolutions))
}

UMAPPlot.1(object, do.print = T,raster=FALSE,cols = Singler.colors,label = T,label.repel = T,group.by = "SCT_snn_res.0.1",no.legend=T)
UMAPPlot.1(object, do.print = T,raster=FALSE,cols = Singler.colors,label = T,label.repel = T,group.by = "SCT_snn_res.0.2",no.legend=T)
UMAPPlot.1(object, do.print = T,raster=FALSE,cols = Singler.colors,label = T,label.repel = T,group.by = "SCT_snn_res.0.5",no.legend=T)


meta.data = object@meta.data
meta.data$X6cluster = meta.data$cell.types
meta.data[meta.data$SCT_snn_res.0.1 %in% c("0","1","11","12","19"),"X6cluster"] = "1"
meta.data[meta.data$SCT_snn_res.0.1 %in% "10","X6cluster"] = "2"
meta.data[meta.data$SCT_snn_res.0.1 %in% "13","X6cluster"] = "3"
meta.data[meta.data$SCT_snn_res.0.1 %in% "14","X6cluster"] = "4"
meta.data[meta.data$SCT_snn_res.0.1 %in% "9","X6cluster"] = "5"
meta.data[meta.data$SCT_snn_res.0.1 %in% "18","X6cluster"] = "6"
meta.data[meta.data$SCT_snn_res.0.1 %in% "16","X6cluster"] = "Pt02"
meta.data$X6cluster %<>% gsub("Monocytes","Monocytes:CD14+",.)
meta.data[meta.data$SCT_snn_res.0.2 %in% "6" &
              meta.data$X6cluster %in% "Monocytes:CD14+","X6cluster"] = "Monocytes:CD16+"

saveRDS(meta.data, file = "data/MCL_61_20220331_metadata.rds")


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
