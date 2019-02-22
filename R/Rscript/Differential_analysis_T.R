
library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
source("../R/Seurat_functions.R")
source("R/util.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file="data/MCL_Harmony_24_20190128.Rda"))
# select 1/4 of cell from control
# in Identify_Cell_Types_Manually.R 2.2
object <- ScaleDown(object = object)

# T cells only ================
object <- SetAllIdent(object, id="res.0.6")
table(object@ident)
T_cells_MCL <- SubsetData(object, ident.use = c(2,3,8))
T_cells_MCL <- SetAllIdent(T_cells_MCL, id="singler1sub")
table(T_cells_MCL@ident)
T_cells_MCL <- SubsetData(T_cells_MCL,
                          ident.remove = c("NK_cells","CLP","MCL","MPP",
                                           "MEP","GMP","B_cells:Memory",
                                           "B_cells:Naive_B_cells"))
table(T_cells_MCL@meta.data$singler1main)
T_cells_MCL <- SetAllIdent(T_cells_MCL, id="singler1main")
T_cells_MCL <- SubsetData(T_cells_MCL, ident.remove = c("Monocytes","MCL"))
T_cells_MCL <- SetAllIdent(T_cells_MCL, id="singler1sub")
T_cells_MCL %<>% FindClusters(reduction.type = "harmony", resolution = 0.3, 
                              dims.use = 1:50,
                              save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                              force.recalc = TRUE, print.output = FALSE)
T_cells_MCL@ident <- plyr::mapvalues(x = T_cells_MCL@ident,
                                     from = c(1,2,3,0),
                                     to = c(1,2,3,4))
T_cells_MCL@ident %<>% factor(levels = 1:4)
T_cells_MCL <- StashIdent(object = T_cells_MCL, save.name = "X4_clusters")
T_cells_MCL@meta.data$X4_orig.ident = paste(T_cells_MCL@meta.data$orig.ident,
                                            T_cells_MCL@meta.data$X4_clusters, sep = "_")

###############################
# Doheatmap for Normal / MCL
###############################
df_samples <- readxl::read_excel("doc/190126_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% paste0("test",2:7))
(samples <- df_samples$sample[sample_n])
# remove samples with low T cells======
table_df <- table(T_cells_MCL@meta.data$orig.ident) %>% as.data.frame
keep <- table_df[table_df$Freq > 100,"Var1"] %>% as.character()
(samples <- samples[samples %in% keep])
T_cells_MCL %<>% SetAllIdent(id = "orig.ident")
for(sample in samples[1]){
    subset.MCL <- SubsetData(T_cells_MCL, ident.use = c(sample,"Normal"))
    #---FindAllMarkers.UMI---- "Keep the shared X4 cluster only"
    subset.MCL %<>% SetAllIdent(id = "X4_orig.ident")
    x4_cluster <- subset.MCL@ident %>% unique
    x4_cluster = x4_cluster[-grep("^Normal",x4_cluster)] %>% as.character %>% 
        gsub('.*\\_',"",.) %>% as.numeric %>% sort
    print(ident.1 <- paste(sample,x4_cluster,sep="_"))
    print(ident.2 <- paste("Normal",x4_cluster,sep="_"))
    
    subset.MCL <- SubsetData(subset.MCL, ident.use = c(ident.1,ident.2))
    subfolder <- paste0(path,sample,"_vs_Normal")
    gde.markers <- FindPairMarkers(subset.MCL, ident.1 = ident.1, 
                                   ident.2 = ident.2,only.pos = FALSE,
                                   logfc.threshold = 0.005,min.cells.group =3,
                                   min.pct = 0.01,
                                   return.thresh = 0.5)
}

###############################
# Doheatmap for MCL.1 / MCL.2
###############################
# remove samples with low T cells======
T_cells_MCL %<>% SetAllIdent(id = "orig.ident")
samples1 <- c("Pt-11-C28","Pt-17-C7","Pt-17-C31","AFT-03-C1D8")
samples2 <- c("Pt-11-C14","Pt-17-C2","Pt-17-C7","AFT-03-C1D1")
for(i in 1:length(samples1)){
    subset.MCL <- SubsetData(T_cells_MCL, ident.use = c(samples1[i],samples2[i]))

    #---FindAllMarkers.UMI---- "Keep the shared X4 cluster only"
    subset.MCL %<>% SetAllIdent(id = "X4_orig.ident")
    print(ident.1 <- paste(samples1[i],1:4,sep="_"))
    print(ident.2 <- paste(samples2[i],1:4,sep="_"))
    subset.MCL <- SubsetData(subset.MCL, ident.use = c(ident.1,ident.2))
    gde.markers <- FindPairMarkers(subset.MCL, ident.1 = ident.1, 
                                   ident.2 = ident.2,only.pos = FALSE,
                                   logfc.threshold = 0.005,min.cells.group =3,
                                   min.pct = 0.01,
                                   return.thresh = 0.5)
}