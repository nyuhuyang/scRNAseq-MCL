invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","magrittr","monocle3",
                   "SeuratData","SeuratWrappers","patchwork"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


#=========================
object <- readRDS(file = "data/MCL_51_20210724.rds")
meta.data <- readRDS(file = "output/MCL_51_20210724_meta.data_v2.rds")

if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}

object %<>% subset(subset =  discard == FALSE)
# CNV
CNV_meta = read.csv("shinyApp/PALIBR_I_51/cnv_meta_data.csv",row.names = 1)
if(all(colnames(object) == rownames(CNV_meta))){
    print("all cellID match!")
    object[["CNV_burden"]] = round(CNV_meta$cnv_score,digits = 6)
    lvl = sort(unique(object$`CNV_burden`))
    object$`CNV_burden` %<>% factor(levels = lvl)
}
object[["harmony.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                    key = "harmonyUMAP_", assay = DefaultAssay(object))
object %<>% subset(subset =  UMAP_1 > 0
              & discard == FALSE
              & tSNE_2 < tSNE_1 +1
              & UMAP_2 <10
              #& orig.ident != "Pt2_30"
              & Doublets == "Singlet"
              #& X4cluster %in% c("1","2","3","4")
              & cell.types  %in% c("B_cells","MCL")
)
npcs <- 80
object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs,k.param = 20)
object %<>% FindClusters(resolution = 0.05, algorithm = 2)
colnames(object@meta.data) %<>% sub("SCT_snn_res.0.05","harmony_Louvain.refine_SCT_snn_res.0.05",.)

cds <- as.cell_data_set(object)

object %<>% FindNeighbors(reduction = "umap",dims = 1:2,k.param = 20)
object %<>% FindClusters(resolution = 0.05, algorithm = 4)
cds <- as.cell_data_set(object)
cds %<>% cluster_cells(reduction_method = "UMAP",cluster_method = "leiden", k=60, random_seed = 101)
plot_cells(cds, show_trajectory_graph = FALSE)
if(all(colnames(object) == names(cds@clusters$UMAP$clusters))){
    object$UMAP_leiden_SCT_snn_res.0.05 <- cds@clusters$UMAP$clusters
}
object$X4cluster_v2 <- plyr::mapvalues(object$harmony_Louvain.refine_SCT_snn_res.0.05,
                                       from = 0:8,
                                       to = c(1,2,3,4,5,5,5,5,5))
object$X4cluster_v3 <- plyr::mapvalues(object$UMAP_leiden_SCT_snn_res.0.05,
                                       from = 1:9,
                                       to = c(1,1,2,3,4,5,5,5,5))
saveRDS(object@meta.data, file = "output/MCL_B_51_20230113_meta.data.rds")

#============= after running DE2 with X4cluster, v1, v2, v3 on PALIBR_I_51_MCL shinyApp ===========
object = readRDS(file = "data/MCL_51_20210724.rds")
meta.data <- readRDS(file = "output/MCL_B_51_20230113_meta.data.rds")

if(all(rownames(meta.data) %in% colnames(object))){
    print("all cellID within!")
    object %<>% subset(cells = rownames(meta.data))
    object@meta.data = meta.data
}
object %<>% subset(subset = X4cluster  %in% c("1","2","3","4"))

excel_files <- list.files(path = "output/20230126",pattern = ".xlsx",full.names = TRUE)
degs_list <- pbapply::pblapply(excel_files,function(FILE){
    readxl::read_excel(path = FILE)
})
Top_n = 150
Cluster_DEGs <- pbapply::pblapply(c("1","2","3","4"),function(cl){
        deg1 <- degs_list[[1]] %>% filter(group == cl & avg_log2FC > 0.1 & p_val_adj < 0.01) %>%
            .[["genes"]]
        deg2 <- degs_list[[2]] %>% filter(group == cl & avg_log2FC > 0.1 & p_val_adj < 0.01) %>%
            .[["genes"]]
        deg3 <- degs_list[[3]] %>% filter(group == cl & avg_log2FC > 0.1 & p_val_adj < 0.01) %>%
            .[["genes"]]
        head(intersect(intersect(deg1,deg2),deg3),Top_n)
}) %>% unlist %>% grep("^MT-",.,invert = TRUE,value = TRUE)
write.table(Cluster_DEGs,file = "output/20230126/top150_degs.csv",row.names = FALSE,quote = FALSE, col.names = FALSE)


Top_n = 40
markers <- c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11")
markers = markers[markers %in% rownames(object)]
degs <- bind_rows(degs_list)
top <- bind_rows(degs_list) %>% group_by(group) %>%
    distinct(genes,avg_log2FC) %>%
    top_n(Top_n, avg_log2FC)
top <- top[grep("^MT-",top$genes,invert = TRUE),]
table(top$group)
features = c(as.character(top$genes),
             tail(VariableFeatures(object = object), 2),
             markers)
sub_obj <- object[features,object$response != "Normal"]


ident <- "X4cluster"
Idents(sub_obj) <- ident
X4 <- sub_obj@meta.data[,ident]
C2 <- sub_obj$orig.ident %in% c("Pt3_BMA72_6","Pt3_72_6") & sub_obj$X4cluster %in% c("1","2","3")
C4 <- sub_obj$orig.ident %in% c("Pt25_25","Pt25_AMB25") & sub_obj$X4cluster %in% "2"

sub_obj@meta.data[C2,ident] <- "2"
sub_obj@meta.data[C4,ident] <- "4"

X4_1 <- sub_obj@meta.data[,ident]
table(X4,X4_1)
Idents(sub_obj) <- ident
sub_obj %<>% ScaleData(features=features)
featuresNum <- make.unique(features, sep = ".")
sub_obj %<>% MakeUniqueGenes(features = features)


DoHeatmap.2(object =sub_obj, group.by = c(ident,"orig.ident"),
            features = featuresNum,
            do.print=T, angle = 45, group.bar = TRUE, title.size = 20, no.legend = FALSE,size=20,hjust = 0.5,
            group.bar.height = 0.02, label=TRUE, cex.row= 8, legend.size = 12,width=14, height=12,
            group1.colors = c('#40A635','#FE8205','#8861AC','#E83C2D'),
            save.path = path,file.name = paste0("Heatmap_top",Top_n,"_",ident,"_orig.ident_legend.jpeg"),
            title = paste("Top",Top_n,"DE genes in 4 B/MCL cells clusters"))

DoHeatmap.2(object =sub_obj, group.by = c(ident,"patient"),
            features = featuresNum,
            do.print=T, angle = 45, group.bar = TRUE, title.size = 20, no.legend = FALSE,size=20,hjust = 0.5,
            group.bar.height = 0.02, label=TRUE, cex.row= 4, legend.size = 12,width=14, height=13,
            group1.colors = c('#40A635','#FE8205','#8861AC','#E83C2D'),
            save.path = path,file.name = paste("Heatmap_top",Top_n,"_",ident,"_patient_legend~.jpeg"),
            title = paste("Top",Top_n,"DE genes in 4 B/MCL cells clusters"))



#Top_n = 150
DoHeatmap.2(object =sub_obj, group.by = c(ident,"orig.ident"),
            features = featuresNum,
            do.print=T, angle = 45, group.bar = TRUE, title.size = 20, no.legend = FALSE,size=20,hjust = 0.5,
            group.bar.height = 0.02, label=TRUE, cex.row= 4.5, legend.size = 12,width=28, height=25,
            group1.colors = c('#40A635','#FE8205','#8861AC','#E83C2D'),
            save.path = path,file.name = paste0("Heatmap_top",Top_n,"_",ident,"_orig.ident_legend.jpeg"),
            title = paste("Top",Top_n,"DE genes in 4 B/MCL cells clusters"))

DoHeatmap.2(object =sub_obj, group.by = c(ident,"patient"),
            features = featuresNum,
            do.print=T, angle = 45, group.bar = TRUE, title.size = 20, no.legend = FALSE,size=20,hjust = 0.5,
            group.bar.height = 0.02, label=TRUE, cex.row= 2, legend.size = 12,width=14, height=13,
            group1.colors = c('#40A635','#FE8205','#8861AC','#E83C2D'),
            save.path = path,file.name = paste("Heatmap_top",Top_n,"_",ident,"_patient_legend~.jpeg"),
            title = paste("Top",Top_n,"DE genes in 4 B/MCL cells clusters"))


UMAPPlot.1(object,group.by = ident,cols = c('#40A635','#FE8205','#8861AC','#E83C2D'),do.print  = T)
TSNEPlot.1(object,group.by = ident,cols = c('#40A635','#FE8205','#8861AC','#E83C2D'),do.print  = T,
           width=8, height=10)
DoHeatmap.1(object =object, features = featuresNum, Top_n = Top_n,group.by = "X4cluster",
            do.print=T, angle = 0, group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0, label=TRUE, cex.row= 3, legend.size = 0,width=6.5, height=10,
            pal_gsea = FALSE,
            save.path = path,file.name = "Heatmap_top40_4cluster.jpeg",
            title = paste("Top",Top_n,"DE genes in 4, B/MCL cells cluster"))

meta.data <- readRDS(file = "output/MCL_B_51_20230113_meta.data.rds")
C2 <- meta.data$orig.ident %in% c("Pt3_BMA72_6","Pt3_72_6") & meta.data$X4cluster %in% c("1","2","3")
C4 <- meta.data$orig.ident %in% c("Pt25_25","Pt25_AMB25") & meta.data$X4cluster %in% "2"

X4 <- meta.data[,"X4cluster"]
meta.data[C2,"X4cluster"] <- "2"
meta.data[C4,"X4cluster"] <- "4"
X4_1 <- meta.data[,"X4cluster"]
table(X4,X4_1)
saveRDS(meta.data, file = "output/MCL_B_51_20230113_meta.data_v2.rds")

#============= after running DE2 with orig.ident_SCT_snn_res.0.05 on PALIBR_I_51_MCL shinyApp ===========
DEGs <- readxl::read_excel("output/20230118/2023-01-18-SCT_snn_res.0.05_0_1_2_3 PtU01_PtU02_PtU03_PtU04_Pt01_Pt02_Pt09_Pt10_Pt11_Pt13_Pt15_Pt16_Pt17_Pt18_Pt19_Pt20_Pt25_Pt27_Pt28_PtB13 .xlsx")
degs <- DEGs %>% group_by(group) %>% arrange(desc(scores), .by_group = TRUE) %>% top_n(125,scores)


length(unique(degs$genes))
with_MT <- unique(degs$genes)
without_MT <- grep("^MT-|^RPL|^RPS",with_MT,value = TRUE,invert = TRUE)
gene_df <- list("with_MT" = with_MT,"without_MT" = without_MT) %>% list2df()
write.csv(gene_df,file = "output/20230118/degs.csv")

#============ after generating heatmap using above unique(degs$genes) for orig.ident_SCT_snn_res.0.05 =======
heatmap <- readxl::read_excel("output/20230119/2023-01-19-Heatmap_orig.ident_SCT_snn_res.0.05.xlsx")
Sample_cl <- colnames(heatmap)
C1 <- Sample_cl[grep("Pt27_1_0",Sample_cl):grep("Pt1_87_3",Sample_cl)]
C2 <- Sample_cl[grep("Pt3_72_6_1",Sample_cl):grep("Pt13_1b_1",Sample_cl)]
C3 <- Sample_cl[grep("PtU04_1",Sample_cl):grep("Pt9_60_2",Sample_cl)]
C4 <- Sample_cl[grep("Pt3_72_6_3",Sample_cl):grep("Pt18_5_8_3",Sample_cl)]

meta.data <- readRDS(file = "output/MCL_B_51_20230113_meta.data.rds")
meta.data$orig.ident_SCT_snn_res.0.05 <- paste0(meta.data$orig.ident,"_",meta.data$SCT_snn_res.0.05)
meta.data$X4cluster_v2 <- plyr::mapvalues(meta.data$orig.ident_SCT_snn_res.0.05,
                                          from = c(C1,C2,C3,C4),
                                          to = c(rep("C1",length(C1)),rep("C2",length(C2)),rep("C3",length(C3)),rep("C4",length(C4))))

meta.data$X4cluster_v2[grep("C1|C2|C3|C4",meta.data$X4cluster_v2,invert = TRUE)] <- "others"
saveRDS(meta.data,"output/MCL_B_51_20230113_meta.data_v2.rds")
