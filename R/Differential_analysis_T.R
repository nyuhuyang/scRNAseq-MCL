########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file="data/MCL_Harmony_12_20181121.Rda"))

# T cells only ================
object <- SetAllIdent(object, id="res.0.6")
table(object@ident)
TSNEPlot.1(object,do.label = T)
T_cells_object <- SubsetData(object, ident.use = c(3,4,5,8))
T_cells_object <- SetAllIdent(T_cells_object, id="singler1sub")
table(T_cells_object@ident)

T_cells_object <- SubsetData(T_cells_object, 
                          ident.remove = c("NK_cells","MCL:2PB1","MCL:8PB1",
                                           "MCL:15PB1","MCL:27PB1",
                                           "MEP","GMP","B_cells:Memory",
                                           "B_cells:Naive_B_cells"))
#table(T_cells_MCL@meta.data$singler1main)
#T_cells_MCL <- SetAllIdent(T_cells_MCL, id="singler1main")
#T_cells_MCL <- SubsetData(T_cells_MCL, ident.use = c("T_cells"))
#T_cells_MCL <- SetAllIdent(T_cells_MCL, id="singler2sub")
p1 <- TSNEPlot.1(T_cells_MCL, do.label = T, do.return = T, pt.size = 0.5, 
                 colors.use = ExtractMetaColor(T_cells_MCL), no.legend =T)
T_cells_MCL %<>% FindClusters(reduction.type = "harmony", resolution = 0.3, 
                              dims.use = 1:50,
                              save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                              force.recalc = TRUE, print.output = FALSE)
p2 <- TSNEPlot.1(T_cells_MCL, do.label = T, do.return = T, pt.size = 0.5, 
                 no.legend =T)
T_cells_MCL <- StashIdent(object = T_cells_MCL, save.name = "4_clusters")

jpeg(paste0(path,"/T_cells_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1, p2)+
        ggtitle("T cells")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))
dev.off()

T_cells_MCL <- SetAllIdent(T_cells_MCL, id="singler2main")
T_cells_MCL <- SubsetData(T_cells_MCL, ident.remove = c("HSC","MCL"))
p3 <- TSNEPlot.1(T_cells_MCL, do.label = T, do.return = T, pt.size = 0.5, 
                 colors.use = ExtractMetaColor(T_cells_MCL), no.legend =T)
jpeg(paste0(path,"/T_cells_TSNEPlot1.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1+ggtitle("Sub cell types")+
                  theme(text = element_text(size=15),							
                        plot.title = element_text(hjust = 0.5,size = 18,
                                                  face = "bold")),
          p3+ggtitle("main cell types")+
                  theme(text = element_text(size=15),							
                        plot.title = element_text(hjust = 0.5,size = 18,
                                                  face = "bold")))
dev.off()
##############################
# SplitTSNEPlot
###############################
table(T_cells_MCL@meta.data$orig.ident)
table(T_cells_MCL@ident)

df_samples <- readxl::read_excel("doc/181227_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
T_cells_MCL <- SetAllIdent(T_cells_MCL, id="singler2sub")
tests <- paste0("test",c(7))
control = c("DJ","MD","BH","NZ")
for(test in tests){
        sample_n = which(df_samples$tests %in% c(control,test))
        samples <- unique(df_samples$sample[sample_n])
        print(paste(c(control,samples), collapse = " "))
        
        cell.use <- rownames(T_cells_MCL@meta.data)[T_cells_MCL@meta.data$orig.ident %in%
                                                            c(control,samples)]
        subset.T_cells_MCL <- SubsetData(T_cells_MCL, cells.use = cell.use)
        
        g <- SplitTSNEPlot(subset.T_cells_MCL,group.by = "ident",split.by = "orig.ident",
                           select.plots = c(1:4,6,5,7),
                           no.legend = T,do.label =F,label.size=3,
                           return.plots =T, label.repel = T,force=2)
        jpeg(paste0(path,test,"_TSNEPlot~.jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, g))
        dev.off()
}

# SingleFeaturePlot.1 for Normal / MCL ================
(markers <-  HumanGenes(MCL,c("CD3D","CD3G","CD8A","CD4", "CD28",
                             "SELL", "CD69", "HLA-DRB1","HLA-DRA")))
(markers <-  HumanGenes(MCL,c("PDCD1","CD274","PDCD1LG2")))
df_samples <- readxl::read_excel("doc/181128_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
tests <- paste0("test",c(3:4))
control = "MD"
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- unique(df_samples$sample[sample_n])
        print(paste(c(control,samples), collapse = " "))
        
        cell.use <- rownames(T_cells_MCL@meta.data)[T_cells_MCL@meta.data$orig.ident %in%
                                                            c(control,samples)]
        subset.T_cells_MCL <- SubsetData(T_cells_MCL, cells.use = cell.use)
        SplitSingleFeaturePlot(subset.T_cells_MCL, 
                               #select.plots = c(1,5,4,3,2),
                               group.by = "ident",split.by = "orig.ident",
                               no.legend = T,label.size=3,do.print =T,
                               markers = markers, threshold = 0.5)
}
# Table with quantification of T cell subsets (absolute # and % of total)
T_cells_MCL <- SetAllIdent(T_cells_MCL, id = "singler2sub")
table(T_cells_MCL@ident,
      T_cells_MCL@meta.data$orig.ident) %>% #prop.table(margin = 2) %>%
        t %>% kable %>% kable_styling
T_cells_MCL.copy <- T_cells_MCL
sub('\\_.*', '', x)
T_cells_MCL.copy@meta.data$singler2sub = gsub('_effector_memory|_central_memory|_Central_memory','',
                        T_cells_MCL.copy@meta.data$singler2sub)
table(T_cells_MCL.copy@meta.data$singler2sub,
      T_cells_MCL@meta.data$orig.ident) %>% prop.table(margin = 2) %>%
        t %>% kable %>% kable_styling

table(T_cells_MCL@meta.data$orig.ident) %>% kable %>% kable_styling
# Doheatmap for Normal / MCL ================
markers <-  HumanGenes(T_cells_MCL,c("CCND1","CCND2","CCND3","CD19","MS4A1","CD79A","CD5","CD40",
                                     "CDK4","CDK6","PCNA","CDK1","SOX11",
                                     "RB1","TP53","ATM","MYC","MTAP",
                                     "FOXO1","FOXO3","CD28"))
(samples <- df_samples$sample[5:7])
control <- "MD"
for(sample in samples){
        cell.use <- rownames(T_cells_MCL@meta.data)[T_cells_MCL@meta.data$orig.ident %in% c(control,sample)] #
        subset.MCL <- SubsetData(T_cells_MCL, cells.use = cell.use)
        #---plot_grid(p1, p2)----
        p1 <- TSNEPlot.1(subset.MCL, do.return = T, pt.size = 0.5, do.label = T, 
                         group.by = "ident",no.legend =T )
        subset.MCL <- SetAllIdent(subset.MCL, id= "singler2sub")
        p2 <- TSNEPlot.1(subset.MCL, do.label = F, do.return = T, pt.size = 0.5, 
                         colors.use = ExtractMetaColor(subset.MCL), no.legend =T)
        jpeg(paste0(path,control,"_",sample,"_","TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
        print(plot_grid(p1, p2)+
                      ggtitle(paste(control, "vs.",sample, "in T cells"))+
                      theme(text = element_text(size=15),							
                            plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) )
        dev.off()
        #---FindAllMarkers.UMI----
        subset.MCL@meta.data$X4_clusters <- paste(subset.MCL@meta.data$orig.ident, 
                                                  subset.MCL@meta.data$X4_clusters, sep = "_")
        subset.MCL <- SetAllIdent(subset.MCL,id = "X4_clusters")
        print(table(subset.MCL@ident))
        #---DoHeatmap.1----
        test_markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.4, only.pos = T,
                                           test.use = "MAST")
        
        g <- DoHeatmap.1(subset.MCL, test_markers, add.genes = markers, Top_n = 25,
                         ident.use = paste(control, "vs.",sample, "in T cells"),
                         group.label.rot = T,cex.row = 2,remove.key =F,title.size = 12)
        jpeg(paste0(path,"/DoHeatmap_",control,"_",sample,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(g)
        dev.off()
        #---DoHeatmap vertical bar----
        g1 <- MakeCorlorBar(subset.MCL, test_markers,Top_n = 25, add.genes = markers, do.print = F,
                            do.return = T)
        jpeg(paste0(path,"/DoHeatmap_",control,"_",sample,"_legend.jpeg"),
             units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
}


# Doheatmap for Normal / MCL.1 / MCL.2 ================
control <- "MD"
tests <- c("test3","test4")
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$sample[sample_n]
        cell.use <- rownames(T_cells_MCL@meta.data)[T_cells_MCL@meta.data$orig.ident %in%
                                                            c(control,samples)]
        subset.MCL <- SubsetData(T_cells_MCL, cells.use = cell.use)
        #---plot_grid(p1, p2)----
        p1 <- TSNEPlot.1(subset.MCL, do.return = T, pt.size = 0.5, do.label = T, 
                         group.by = "ident",no.legend =T )
        subset.MCL <- SetAllIdent(subset.MCL, id= "singler2sub")
        p2 <- TSNEPlot.1(subset.MCL, do.label = F, do.return = T, pt.size = 0.5, 
                         colors.use = ExtractMetaColor(subset.MCL), no.legend =T)
        jpeg(paste0(path,paste0(samples,collapse = "_"),"TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
        print(plot_grid(p1, p2)+
                      ggtitle(paste(control, "vs.",paste0(samples,collapse = " vs. "), "in T cells"))+
                      theme(text = element_text(size=15),							
                            plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) )
        dev.off()
        #---FindAllMarkers.UMI----
        subset.MCL@meta.data$X4_clusters <- paste(subset.MCL@meta.data$orig.ident, 
                                                  subset.MCL@meta.data$X4_clusters, sep = ".")
        subset.MCL <- SetAllIdent(subset.MCL,id = "X4_clusters")
        table(subset.MCL@ident)
        test_markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.4, only.pos = T,
                                           test.use = "MAST")
        subset.MCL <- SetAllIdent(subset.MCL, id="4_clusters")
        #---DoHeatmap.1----
        g <- DoHeatmap.1(subset.MCL, test_markers, add.genes = markers, Top_n = 15,
                         ident.use = paste(control, "vs.",paste0(samples,collapse = " & "), "in MCL"),
                         group.label.rot = T,cex.row = 3,remove.key =F,
                         title.size = 12)
        jpeg(paste0(path,"/DoHeatmap_",control,"_",paste0(samples,collapse = "_"),".jpeg"),
             units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
        #---DoHeatmap vertical bar----
        g1 <- MakeCorlorBar(subset.MCL, test_markers,Top_n = 15, add.genes = markers, do.print = F,
                            do.return = T)
        jpeg(paste0(path,"/DoHeatmap_",control,"_",paste0(samples,collapse = "_"),"legend.jpeg"),
             units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
}
