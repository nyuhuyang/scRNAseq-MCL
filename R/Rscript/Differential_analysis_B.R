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
B_cells_MCL <- SubsetData(object, ident.use = c(0,1,4,5,6,9))
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1main")
table(B_cells_MCL@ident)
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B_cells","HSC","MCL"))
table(B_cells_MCL@meta.data$singler1sub)
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1sub")
(ident.remove <- grep("T_cells.",B_cells_MCL@meta.data$singler1sub, value = T) %>% unique)
B_cells_MCL <- SubsetData(B_cells_MCL, ident.remove =  ident.remove)

B_cells_MCL %<>% FindClusters(reduction.type = "harmony", resolution = 0.2, 
                              dims.use = 1:50,
                              save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                              force.recalc = TRUE, print.output = FALSE)
B_cells_MCL@ident <- plyr::mapvalues(x = B_cells_MCL@ident,
                                     from = c(0,1,2,3,4,5),
                                     to = c(1,3,4,5,2,6))
B_cells_MCL@ident %<>% factor(levels = 1:6)
B_cells_MCL <- StashIdent(object = B_cells_MCL, save.name = "X6_clusters")
B_cells_MCL@meta.data$X6_orig.ident = paste(B_cells_MCL@meta.data$orig.ident,
                                            B_cells_MCL@meta.data$X6_clusters, sep = "_")
B_cells_MCL@meta.data$X6_orig.ident = gsub('^Normal_.*', 'Normal', B_cells_MCL@meta.data$X6_orig.ident)
remove(object);GC();GC();GC();GC();GC();GC();GC();GC();GC();
###############################
# Doheatmap for Normal / MCL
###############################
df_samples <- readxl::read_excel("doc/190126_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% paste0("test",2:7))
(samples <- df_samples$sample[sample_n])
# remove samples with low T cells======
table_df <- table(B_cells_MCL@meta.data$orig.ident) %>% as.data.frame
keep <- table_df[table_df$Freq > 100,"Var1"] %>% as.character()
(samples <- samples[samples %in% keep])
B_cells_MCL %<>% SetAllIdent(id = "orig.ident")
for(sample in samples[1]){
    subset.MCL <- SubsetData(B_cells_MCL, ident.use = c(sample,"Normal"))
    #---FindAllMarkers.UMI---- "Keep the shared X4 cluster only"
    subset.MCL %<>% SetAllIdent(id = "X6_orig.ident")
    x6_cluster <- subset.MCL@ident %>% unique
    x6_cluster = x6_cluster[-grep("^Normal",x6_cluster)] %>% as.character %>% 
        gsub('.*\\_',"",.) %>% as.numeric %>% sort
    print(ident.1 <- paste(sample,x6_cluster,sep="_"))
    print(ident.2 <- rep("Normal",length(ident.1)))
    
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
B_cells_MCL %<>% SetAllIdent(id = "orig.ident")
samples1 <- c("Pt-11-C28","Pt-17-C7","Pt-17-C31","AFT-03-C1D8")
samples2 <- c("Pt-11-C14","Pt-17-C2","Pt-17-C7","AFT-03-C1D1")
for(i in 1:length(samples1)){
    subset.MCL <- SubsetData(B_cells_MCL, ident.use = c(samples1[i],samples2[i]))

    #---FindAllMarkers.UMI---- "Keep the shared X4 cluster only"
    subset.MCL %<>% SetAllIdent(id = "X6_orig.ident")
    print(ident.1 <- paste(samples1[i],1:4,sep="_"))
    print(ident.2 <- paste(samples2[i],1:4,sep="_"))
    subset.MCL <- SubsetData(subset.MCL, ident.use = c(ident.1,ident.2))
    gde.markers <- FindPairMarkers(subset.MCL, ident.1 = ident.1, 
                                   ident.2 = ident.2,only.pos = FALSE,
                                   logfc.threshold = 0.005,min.cells.group =3,
                                   min.pct = 0.01,
                                   return.thresh = 0.5)
}
