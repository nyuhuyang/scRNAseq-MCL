path <- paste0("output/",gsub("-","",Sys.Date()),"/UMAP")
if(!dir.exists(path)) dir.create(path, recursive = T)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

object = readRDS(file = "data/MCL_61_20220318.rds")

test_df = data.frame(min_dist = rep(c(0.2,0.4,0.6),each = 3),
                     spread = rep(c(0.6,1.0,1.4),times = 3),
                     npcs = rep(c(50,60,70,80,90,100),each = 9))
file.names = paste0("cs",test_df$npcs,"_dist.",test_df$min_dist,"_spread.",test_df$spread)

for (i in 1:nrow(test_df)){
    print(test_df[i,])
    file.name = file.names[i]
    if(!file.exists(paste0("output/20220331/Pt22_67/umap_",file.name,".rds"))) next
    umap = readRDS(paste0("output/20220331/Pt22_67/umap_",file.name,".rds"))
    object@reductions[[paste0("umap_",file.name)]] = umap
    object@reductions[["umap"]] =  umap[[1]]
    colnames(object@reductions[["umap"]]@cell.embeddings) %<>% gsub(".*UMAP","UMAP",.)
    object@reductions[["umap"]]@key %<>% gsub(".*UMAP","UMAP",.)

        UMAPPlot.1(object, group.by = "cell.types",do.print = T,
               title = file.name,file.name = paste0(file.name,".jpeg"),
               save.path = path)
}
