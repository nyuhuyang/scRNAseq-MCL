
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
NK <- SubsetData(object, ident.use = c(2))
NK <- SetAllIdent(NK, id="singler1sub")
table(NK@ident)
NK <- SubsetData(NK,ident.use = c("NK_cells"))
table(NK@meta.data$singler1main)
NK <- SetAllIdent(NK, id="singler1sub")
NK %<>% FindClusters(reduction.type = "harmony", resolution = 0.3, 
                              dims.use = 1:50,
                              save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                              force.recalc = TRUE, print.output = FALSE)
NK@ident <- plyr::mapvalues(x = NK@ident,
                            from = c(0,1),
                            to = c(1,2))
NK <- StashIdent(object = NK, save.name = "X_clusters")
NK@meta.data$X_orig.ident = paste(NK@meta.data$orig.ident,
                                            NK@meta.data$X_clusters, sep = "_")

###############################
# DE for Normal / MCL
###############################
df_samples <- readxl::read_excel("doc/190126_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% paste0("test",2:7))
(samples <- df_samples$sample[sample_n])
# remove samples with low T cells======
table_df <- table(NK@meta.data$orig.ident) %>% as.data.frame
keep <- table_df[table_df$Freq > 10,"Var1"] %>% as.character()
(samples <- samples[samples %in% keep])
NK %<>% SetAllIdent(id = "orig.ident")
for(sample in samples){
    subset.MCL <- SubsetData(NK, ident.use = c(sample,"Normal"))
    #---FindAllMarkers.UMI---- "Keep the shared X4 cluster only"
    subset.MCL %<>% SetAllIdent(id = "X_orig.ident")
    x4_cluster <- subset.MCL@ident %>% unique
    x4_cluster = x4_cluster[-grep("^Normal",x4_cluster)] %>% as.character %>%
        gsub('.*\\_',"",.) %>% as.numeric %>% sort
    print(ident.1 <- paste(sample,x4_cluster,sep="_"))
    print(ident.2 <- paste("Normal",x4_cluster,sep="_"))
    
    #subset.MCL <- SubsetData(subset.MCL, ident.use = c(ident.1,ident.2))
    subfolder <- paste0(path,"NK/",sample,"_vs_Normal/")
    gde.markers <- FindPairMarkers(subset.MCL, ident.1 = ident.1,
                                   ident.2 = ident.2,only.pos = FALSE,
                                   logfc.threshold = -Inf,min.cells.group =3,
                                   min.pct = -Inf,
                                   return.thresh = Inf,
                                   save.path = subfolder)
}

###############################
# DE for MCL.1 / MCL.2
###############################
# remove samples with low T cells======
NK %<>% SetAllIdent(id = "orig.ident")
samples1 <- c("Pt-11-C28","Pt-17-C7","Pt-17-C31","AFT-04-C1D8")
samples2 <- c("Pt-11-C14","Pt-17-C2","Pt-17-C7","AFT-04-C1D1")
for(i in 1:length(samples1)){
    subset.MCL <- SubsetData(NK, ident.use = c(samples1[i],samples2[i]))

    #---FindAllMarkers.UMI---- "Keep the shared X4 cluster only"
    subset.MCL %<>% SetAllIdent(id = "X_orig.ident")
    print(ident.1 <- paste(samples1[i],1:2,sep="_"))
    print(ident.2 <- paste(samples2[i],1:2,sep="_"))
    #subset.MCL <- SubsetData(subset.MCL, ident.use = c(ident.1,ident.2))
    subfolder <- paste0(path,"NK/",samples1[i],"_vs_",samples2[i],"_T/")
    gde.markers <- FindPairMarkers(subset.MCL, ident.1 = ident.1, 
                                   ident.2 = ident.2,only.pos = FALSE,
                                   logfc.threshold = -Inf,min.cells.group =3,
                                   min.pct = -Inf,
                                   return.thresh = Inf,
                                   save.path = subfolder)
}
