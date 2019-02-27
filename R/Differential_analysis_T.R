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
TSNEPlot.1(object,do.label = T)
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
p1 <- TSNEPlot.1(T_cells_MCL, do.label = T, do.return = T, pt.size = 0.5, 
                 colors.use = ExtractMetaColor(T_cells_MCL), no.legend =T)
T_cells_MCL %<>% FindClusters(reduction.type = "harmony", resolution = 0.3, 
                              dims.use = 1:50,
                              save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                              force.recalc = TRUE, print.output = FALSE)
T_cells_MCL@ident <- plyr::mapvalues(x = T_cells_MCL@ident,
                                from = c(1,2,3,0),
                                to = c(1,2,3,4))
T_cells_MCL@ident %<>% factor(levels = 1:4)

p2 <- TSNEPlot.1(T_cells_MCL, do.label = T, do.return = T, pt.size = 0.5, 
                 no.legend =T)
T_cells_MCL <- StashIdent(object = T_cells_MCL, save.name = "X4_clusters")

jpeg(paste0(path,"/T_cells_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1, p2)+
        ggtitle("T cells")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))
dev.off()

T_cells_MCL <- SetAllIdent(T_cells_MCL, id="singler1main")

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
#T_cells_MCL@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",T_cells_MCL@meta.data$orig.ident)
T_cells_MCL %<>% SetAllIdent(id="orig.ident")

df_samples <- readxl::read_excel("doc/190126_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower

tests <- paste0("test",2:7)
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df <- as.data.frame(df_samples[sample_n,])
        samples <- unique(df$sample)
        rownames(df) = samples
        
        print(samples <- c("Normal",df$sample[order(df$tsne)]))
        
        g <- lapply(samples,function(sample) {
                SubsetData(T_cells_MCL, ident.use = sample) %>%
                        SetAllIdent(id = "singler1sub") %>%
                        TSNEPlot.1(no.legend = T,do.label =F,label.size=3,size=20,
                                   colors.use = ExtractMetaColor(.),
                                   return.plots =T, label.repel = T,force=2)+
                        ggtitle(sample)+theme(text = element_text(size=20),
                                              plot.title = element_text(hjust = 0.5))
        })
        jpeg(paste0(path,test,"_Tcell_TSNEPlot.jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, c(g, nrow = 2)))
        dev.off()
}

##############################
# Split UMI
###############################
T_cells_MCL@meta.data$X4_orig.ident = paste(T_cells_MCL@meta.data$orig.ident,
                                            T_cells_MCL@meta.data$X4_clusters, sep = "_")
T_cells_MCL %<>% SetAllIdent(id = "X4_orig.ident")
All_AveExp <- AverageExpression(T_cells_MCL)
write.csv(All_AveExp,file = paste0(path,"All_AveExp.csv"))

VarGenes_AveExp <- AverageExpression(T_cells_MCL,genes.use = T_cells_MCL@var.genes)
write.csv(VarGenes_AveExp,file = paste0(path,"VarGenes_AveExp.csv"))

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
for(sample in samples){
    subset.MCL <- SubsetData(T_cells_MCL, ident.use = c(sample,"Normal"))
    #---FindAllMarkers.UMI---- "Keep the shared X4 cluster only"
    subset.MCL %<>% SetAllIdent(id = "X4_orig.ident")
    x4_cluster <- subset.MCL@ident %>% unique
    x4_cluster = x4_cluster[-grep("^Normal",x4_cluster)] %>% as.character %>% 
        gsub('.*\\_',"",.) %>% as.numeric %>% sort
    print(ident.1 <- paste(sample,x4_cluster,sep="_"))
    print(ident.2 <- paste("Normal",x4_cluster,sep="_"))

    subset.MCL <- SubsetData(subset.MCL, ident.use = c(ident.1,ident.2))
    gde.markers <- FindPairMarkers(subset.MCL, ident.1 = ident.1, 
                                   ident.2 = ident.2,
                                   logfc.threshold = 0.005,min.cells.group =3,
                                   min.pct = 0.01,
                                   return.thresh = 0.5)
}
    (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
    if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
    GC();GC();GC();GC();GC();GC();GC();GC();
    subset.MCL <- ScaleData(subset.MCL)

    #---SplitTSNEPlot---- "must keep the order"
    #subset.MCL <- SetAllIdent(subset.MCL, id= "singler1sub")
    #SplitTSNEPlot(subset.MCL,do.return = FALSE,do.print = TRUE, do.label = F)
    #subset.MCL <- SetAllIdent(subset.MCL,id = "X4_clusters")
    #SplitTSNEPlot(subset.MCL,do.return = FALSE,do.print = TRUE, do.label = F)
    for(meta in c("singler1sub","X4_orig.ident")){
        g <- lapply(c("Normal",sample),function(s) {
            SubsetData(subset.MCL, ident.use = s) %>%
                SetAllIdent(id = meta) %>%
                TSNEPlot.1(no.legend = T, do.label =F, label.size=3,size=20,
                           colors.use = ExtractMetaColor(.),
                           return.plots =T, label.repel = T,force=2)+
                ggtitle(s)+theme(text = element_text(size=20),
                                 plot.title = element_text(hjust = 0.5))
        })
        jpeg(paste0(path,"N_",sample,"_",meta,"_Tcell_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
        print(do.call(plot_grid, c(g, nrow = 1)))
        dev.off()
    }}
    #---DoHeatmap.1----
    subset.MCL@ident = factor(subset.MCL@ident, levels = c(dent.1,dent.2))
    markers <-  HumanGenes(subset.MCL,c("CD3D","CD3G","CD4","CD8A","FCGR3A","FCGR3B","NCAM1","KLRC1"))
    g <- DoHeatmap.1(subset.MCL, gde.markers, add.genes = markers, Top_n = 25,
                     use.scaled = T,
                     ident.use = paste("Normal vs.",sample, "in T cells"),
                     group.label.rot = T,cex.row = 3,remove.key =F,title.size = 12)
    jpeg(paste0(path,"DoHeatmap_Normal_",sample,".jpeg"), units="in", width=10, height=7,
         res=600)
    print(g)
    dev.off()
    #---DoHeatmap vertical bar----
    g1 <- MakeCorlorBar(subset.MCL, gde.markers,Top_n = 25, add.genes = markers, do.print = F,
                        do.return = T)
    jpeg(paste0(path,"DoHeatmap_Normal_",sample,"_legend.jpeg"),
         units="in", width=10, height=7,res=600)
    print(g1)
    dev.off()
}

###############################
# Doheatmap for MCL.1 / MCL.2
###############################
# remove samples with low T cells======
T_cells_MCL %<>% SetAllIdent(id = "orig.ident")
samples1 <- c("Pt-11-C14","Pt-17-C2","Pt-17-C7","AFT-03-C1D1")
samples2 <- c("Pt-11-C28","Pt-17-C7","Pt-17-C31","AFT-03-C1D8")

for(i in 1:length(samples1)){
    subset.MCL <- SubsetData(T_cells_MCL, ident.use = c(samples1[i],samples2[i]))
    #---SplitTSNEPlot---- "must keep the order"
    #subset.MCL <- SetAllIdent(subset.MCL, id= "singler1sub")
    #SplitTSNEPlot(subset.MCL,do.return = FALSE,do.print = TRUE, do.label = F)
    #subset.MCL <- SetAllIdent(subset.MCL,id = "X4_clusters")
    #SplitTSNEPlot(subset.MCL,do.return = FALSE,do.print = TRUE, do.label = F)
    for(meta in c("singler1sub","X4_orig.ident")){
        g <- lapply(c(samples1[i],samples2[i]),function(s) {
            SubsetData(subset.MCL, ident.use = s) %>%
                SetAllIdent(id = meta) %>%
                TSNEPlot.1(no.legend = T, do.label =F, label.size=3,size=20,
                           colors.use = ExtractMetaColor(.),
                           return.plots =T, label.repel = T,force=2)+
                ggtitle(s)+theme(text = element_text(size=20),
                                 plot.title = element_text(hjust = 0.5))
        })
        jpeg(paste0(path,paste(samples1[i],samples2[i],collapse = "_"),
                    "_",meta,"_Tcell_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
        print(do.call(plot_grid, c(g, nrow = 1)))
        dev.off()
    }
    #---FindAllMarkers.UMI---- "Keep the shared X4 cluster only"
    subset.MCL %<>% SetAllIdent(id = "X4_orig.ident")
    print(dent.1 <- paste(samples1[i],1:4,sep="_"))
    print(dent.2 <- paste(samples2[i],1:4,sep="_"))
    subset.MCL <- SubsetData(subset.MCL, ident.use = c(dent.1,dent.2))
    gde.markers <- FindPairMarkers(subset.MCL, ident.1 = c(dent.1,dent.2), 
                                   ident.2 = c(dent.2,dent.1),
                                   logfc.threshold = 0.2,min.cells.group =3,
                                   return.thresh = 0.05)
    write.csv(gde.markers, paste0(path,samples1[i],"_",samples2[i],".csv"))
    (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
    if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
    GC();GC();GC();GC();GC();GC();GC();GC();
    subset.MCL <- ScaleData(subset.MCL)
    #---DoHeatmap.1----
    subset.MCL@ident = factor(subset.MCL@ident, levels = c(dent.1,dent.2))
    markers <-  HumanGenes(subset.MCL,c("CD3D","CD3G","CD4","CD8A","FCGR3A","FCGR3B","NCAM1","KLRC1"))
    g <- DoHeatmap.1(subset.MCL, gde.markers, add.genes = markers, Top_n = 25,
                     use.scaled = T,
                     ident.use = paste(samples1[i],"vs.",samples2[i], "in T cells"),
                     group.label.rot = T,cex.row = 3,remove.key =F,title.size = 12)
    jpeg(paste0(path,"DoHeatmap_Normal_",samples1[i],"_",samples2[i],".jpeg"), units="in", width=10, height=7,
         res=600)
    print(g)
    dev.off()
    #---DoHeatmap vertical bar----
    g1 <- MakeCorlorBar(subset.MCL, gde.markers,Top_n = 25, add.genes = markers, do.print = F,
                        do.return = T)
    jpeg(paste0(path,"DoHeatmap_Normal_",samples1[i],"_",samples2[i],"_legend.jpeg"),
         units="in", width=10, height=7,res=600)
    print(g1)
    dev.off()
}

##############################
# SingleFeaturePlot.1 for Normal / MCL
###############################
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

object@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",object@meta.data$orig.ident)
object %<>% SetAllIdent(id="orig.ident")
Normal <- SubsetData(object, ident.use= "Normal")
Normal_Exp <- AverageExpression(Normal)



#========
# Add FoldChange column ===========
output <- "/Users/yah2014/Downloads/20190124_MCL/"
files <- list.files(output)
(files =files[grep(".csv",files)])
for(file in files){
    gde.markers = read.csv(paste0(output,file))
    gde.markers$avg_log2FC = log2(exp(1)) * gde.markers$avg_logFC
    gde.markers = gde.markers[,c( "gene","p_val","avg_log2FC","FC","pct.1","pct.2",
                        "p_val_adj","UMI.1","UMI.2","cluster1.vs.cluster2")]
    write.csv(gde.markers, paste0(output,file))
}

