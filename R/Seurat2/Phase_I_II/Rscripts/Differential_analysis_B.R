 invisible(sapply(c("Seurat","dplyr","tidyr","magrittr","dplyr","gplots"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
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
(load(file="data/MCL_Harmony_36_20190420.Rda"))
# select 1/4 of cell from control
# B cells only ================
object <- SetAllIdent(object, id="res.0.6")
table(object@ident)
#TSNEPlot(object,do.label = T)
B_cells_MCL <- SubsetData(object, ident.use = c(0,1,6,7,10,12,13,14,17))
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1main")
table(B_cells_MCL@ident)
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B_cells","MCL"))
table(B_cells_MCL@meta.data$singler1sub) %>% as.data.frame %>%
.[.[,"Freq"] >0,]
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1sub")
(keep <- grep("B_cells.|MCL",B_cells_MCL@meta.data$singler1sub, value = T) %>% unique)
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use =  keep)

g_Harmony <- TSNEPlot.1(object = B_cells_MCL, do.label = T, group.by = "ident",
                        do.return = TRUE, no.legend = T,
                        colors.use = ExtractMetaColor(B_cells_MCL),
                        pt.size = 1,label.size = 5 )+
                        ylim(-25,25)+
                ggtitle("Tsne plot of all B and MCL cells")+
                theme(plot.title = element_text(hjust = 0.5,size = 18))

jpeg(paste0(path,"TSNEplot-B_cells_MCL_sub.jpeg"), units="in", width=10, height=7,res=600)
print(g_Harmony)
dev.off()

#remove(object);GC()
system.time(
            B_cells_MCL <- RunTSNE(B_cells_MCL, reduction.use = "harmony", dims.use = 1:75,
            perplexity = 30, do.fast = TRUE))
system.time(
            B_cells_MCL %<>% FindClusters(reduction.type = "harmony", resolution = 0.2, dims.use = 1:75,
            save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
            force.recalc = TRUE, print.output = FALSE))
TSNEPlot(B_cells_MCL,do.label = T,label.size=5,do.print=T)
B_cells_MCL@ident <- plyr::mapvalues(x = B_cells_MCL@ident,
                                     from = c(0,1,2,3,4,5,6,7,8,9,10),
                                     to = c(1,3,4,2,5,5,5,5,5,5,5))
B_cells_MCL@ident %<>% factor(levels = 1:5)
B_cells_MCL <- StashIdent(object = B_cells_MCL, save.name = "X5_clusters")
B_cells_MCL@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",B_cells_MCL@meta.data$orig.ident)
B_cells_MCL@meta.data$X5_orig.ident = paste(B_cells_MCL@meta.data$orig.ident,
                                            B_cells_MCL@meta.data$X5_clusters, sep = "_")
B_cells_MCL@meta.data$X5_orig.ident = gsub('^Normal_.*', 'Normal', B_cells_MCL@meta.data$X5_orig.ident)
###############################
# Doheatmap for Normal / MCL
###############################
df_samples <- readxl::read_excel("doc/190406_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% paste0("test",8))
(samples <- df_samples$sample[sample_n])
# remove samples with low B cells======
table_df <- table(B_cells_MCL@meta.data$orig.ident) %>% as.data.frame
keep <- table_df[table_df$Freq > 100,"Var1"] %>% as.character()
(samples <- samples[samples %in% keep])
B_cells_MCL %<>% SetAllIdent(id = "orig.ident")
for(sample in samples[1:length(samples)]){
    subset.MCL <- SubsetData(B_cells_MCL, ident.use = c(sample,"Normal"))
    subset.MCL %<>% SetAllIdent(id = "X5_orig.ident")
    # remove cluster with less than 3 cells======
    table_subset.MCL <- table(subset.MCL@meta.data$X5_orig.ident) %>% as.data.frame
    keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character()
    subset.MCL <- SubsetData(subset.MCL, ident.use = keep.MCL)

    x6_cluster <- subset.MCL@ident %>% unique
    x6_cluster = x6_cluster[-grep("^Normal",x6_cluster)] %>% as.character %>% 
        gsub('.*\\_',"",.) %>% as.numeric %>% sort
    print(ident.1 <- paste(sample,x6_cluster,sep="_"))
    print(ident.2 <- rep("Normal",length(ident.1)))
    
    subfolder <- paste0(path,sample,"_vs_Normal/")

    gde.markers <- FindPairMarkers(subset.MCL, ident.1 = ident.1, 
                                   ident.2 = ident.2,only.pos = T,
                                   logfc.threshold = 1.005,min.cells.group =3,
                                   min.pct = 0.01,
                                   return.thresh = 0.05,
                                   save.path = subfolder)
}

###############################
# Doheatmap for MCL.1 / MCL.2
###############################
# remove samples with low B cells======
B_cells_MCL %<>% SetAllIdent(id = "orig.ident")
samples1 <- c("Pt-11-C1","Pt-11-C14","Pt-11-C28","Pt-11-C28","Pt-17-C2","Pt-17-C7","Pt-17-C31","Pt-17-C31","AFT-03-C1D8","Pt-AA13-Ib-1")
samples2 <- c("Pt-11-LN-C1","Pt-11-C1","Pt-11-C1","Pt-11-C14","Pt-17-LN-C1","Pt-17-C2","Pt-17-C2","Pt-17-C7","AFT-03-C1D1","Pt-AA13-Ib-p")

for(i in 1:length(samples1)){
    subset.MCL <- SubsetData(B_cells_MCL, ident.use = c(samples1[i],samples2[i]))
    subset.MCL %<>% SetAllIdent(id = "X5_orig.ident")
    # remove cluster with less than 3 cells======
    table_subset.MCL <- table(subset.MCL@meta.data$X5_orig.ident) %>% as.data.frame
    keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character()
    subset.MCL <- SubsetData(subset.MCL, ident.use = keep.MCL)
    
    x6_cluster <- subset.MCL@ident %>% unique %>% 
        gsub('.*\\_',"",.) %>% as.numeric %>% sort %>% .[duplicated(.)]
    
    print(ident.1 <- paste(samples1[i],x6_cluster,sep="_"))
    print(ident.2 <- paste(samples2[i],x6_cluster,sep="_"))

    subfolder <- paste0(path,"20190222_B/",samples1[i],"_vs_",samples2[i],"/")
    gde.markers <- FindPairMarkers(subset.MCL, ident.1 = ident.1, 
                                   ident.2 = ident.2,only.pos = FALSE,
                                   logfc.threshold = 0.005,min.cells.group =3,
                                   min.pct = 0.01,
                                   return.thresh = 0.5,
                                   save.path = subfolder)
}
