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
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file="data/MCL_V3_Harmony_43_20190610.Rda"))
# T cells only ================
object <- SetAllIdent(object, id="singler1sub")
TSNEPlot.1(object,do.label = T,no.legend=T,col.use = ExtractMetaColor(object),do.print = T)

object <- SetAllIdent(object, id="res.0.6")
table(object@ident)
B_cells_MCL <- SubsetData(object, ident.use = c(0,1,5,6,9,13,14,18,19,20))
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1main")
table(B_cells_MCL@ident)
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B_cells","MCL","HSC"))
#table(B_cells_MCL@meta.data$singler1sub) %>% as.data.frame %>%
#        .[.[,"Freq"] >0,]
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1sub") 
(keep <- grep("B_cells.|MCL|MEP|CLP|HSC|CMP|GMP",B_cells_MCL@meta.data$singler1sub, value = T) %>% unique)
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use =  keep)

B_cells_MCL@meta.data$cell.type = gsub("B_cells.*","B_cells",B_cells_MCL@meta.data$singler1sub)
B_cells_MCL@meta.data$cell.type = gsub("MEP|CLP|HSC|CMP|GMP","MCL/HSC",B_cells_MCL@meta.data$cell.type)
table(B_cells_MCL@meta.data$cell.type)

B_cells_MCL@meta.data$cell.type.colors = gsub("B_cells","#33A02C",B_cells_MCL@meta.data$cell.type)
B_cells_MCL@meta.data$cell.type.colors = gsub("MCL/HSC","#4038b0",B_cells_MCL@meta.data$cell.type.colors)
B_cells_MCL@meta.data$cell.type.colors = gsub("MCL","#ff0000",B_cells_MCL@meta.data$cell.type.colors)

table(B_cells_MCL@meta.data$cell.type.colors)


B_cells_MCL <- SetAllIdent(B_cells_MCL, id="cell.type")
g <- TSNEPlot.1(object = B_cells_MCL, do.label = F, group.by = "ident",
                do.return = TRUE, no.legend = F, 
                colors.use = ExtractMetaColor(B_cells_MCL),
                pt.size = 1,label.size = 3 )+
        ggtitle("Tsne plot of all B and MCL cells")+
        theme(plot.title = element_text(hjust = 0.5,size = 18)) 

jpeg(paste0(path,"TSNEplot-B_cells.jpeg"), units="in", width=10, height=7,res=600)
print(g)
dev.off()

TSNEPlot.1(object,do.label = F,no.legend=F,do.print = T)


##############################
# re-scale, PCA, harmony and tsne
##############################
#B_cells_MCL <- NormalizeData(B_cells_MCL)
#B_cells_MCL <- FindVariableGenes(object = B_cells_MCL, mean.function = ExpMean, 
#                            dispersion.function = LogVMR, 
#                            x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.5)
#length(B_cells_MCL@var.genes)

#B_cells_MCL %<>% ScaleData
#B_cells_MCL %<>% RunPCA(pc.genes = B_cells_MCL@var.genes, pcs.compute = 100, do.print = F)

#pcs = 1:75
#jpeg(paste0(path,"S1_RunHarmony_B.jpeg"), units="in", width=10, height=7,res=600)
#system.time(B_cells_MCL %<>% RunHarmony("orig.ident", dims.use = pcs,
#                                   theta = 2, plot_convergence = TRUE,
#                                   nclust = 50, max.iter.cluster = 100))
#dev.off()

#system.time(
#        B_cells_MCL <- RunTSNE(B_cells_MCL, reduction.use = "harmony", dims.use = 1:75,
#                                  perplexity = 30, do.fast = TRUE))
system.time(
        B_cells_MCL %<>% FindClusters(reduction.type = "harmony", resolution = 0.8, dims.use = 1:75,
                                 save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                                 force.recalc = TRUE, print.output = FALSE))

B_cells_MCL <- SetAllIdent(B_cells_MCL, id="res.0.8")
TSNEPlot.1(B_cells_MCL,do.label = F,#colors.use = ExtractMetaColor(B_cells_MCL),
           no.legend = F,label.size=5,do.print=T,do.return = F)
#save(B_cells_MCL, file = "data/B_cells_MCL_43_20190519.Rda")

#(load("data/B_cells_MCL_43_20190519.Rda"))
meta.data = cbind.data.frame(B_cells_MCL@meta.data, B_cells_MCL@dr$tsne@cell.embeddings)
meta.data[(meta.data$tSNE_1 > -10) & (meta.data$res.0.6 == 9),"res.0.6"]=21
meta.data[(meta.data$tSNE_1 < -17.5) & (meta.data$res.0.6 == 9),"res.0.6"]=21
meta.data[(meta.data$res.0.8 %in% 5),"res.0.6"]=22
meta.data[(meta.data$tSNE_2 > -15) & (meta.data$res.0.8 %in% 5),"res.0.6"]=6
B_cells_MCL@meta.data = meta.data

B_cells_MCL <- SetAllIdent(B_cells_MCL, id="res.0.6")
B_cells_MCL <- SubsetData(B_cells_MCL, ident.remove = 21)

g1 <- TSNEPlot.1(B_cells_MCL,do.label = T, do.return=T)
B_cells_MCL@meta.data$X6_clusters <- plyr::mapvalues(x = B_cells_MCL@meta.data$res.0.6,
                                                     from = c(0,1,5,6,9,13,14,18,19,20,22),
                                                     to =   c(1,2,3,4,5,2, 1, 1, 2, 1,6))


B_cells_MCL <- SetAllIdent(B_cells_MCL, id="X6_clusters")
B_cells_MCL@ident %<>% factor(levels = 1:6)

# tsne plot
g2 <- TSNEPlot.1(B_cells_MCL, do.label = T, do.return=T)
jpeg(paste0(path,"rename_tsne_B.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g1,g2) +  ggtitle("Rename the clusters")+
        theme(text = element_text(size=15),
              plot.title = element_text(hjust = 0.5))
dev.off()
SingleFeaturePlot.1(B_cells_MCL, feature = "CD274", threshold = NULL,do.print = T)

save(B_cells_MCL, file = "data/B_cells_MCL_43_20190608.Rda")
(load(file = "data/B_cells_MCL_43_20190603.Rda"))

######
cluster_6 <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$X6_clusters %in% 6]
write.csv(cluster_6,paste0(path,"cluster_6_cell_names.csv"))

# rename cluster6
cluster_6_cell_names = read.csv(file="output/20190608/cluster_6_cell_names.csv",row.names=1)
cluster_6_cell_names <- as.character(cluster_6_cell_names$x)
meta.data = B_cells_MCL@meta.data
meta.data$X6_clusters = meta.data$X5_clusters     
meta.data[cluster_6_cell_names,"X6_clusters"]=6            
table(meta.data$X6_clusters)
B_cells_MCL@meta.data = meta.data
##############
# DE genes between Clusters 5 cluster top 50
###############
B_cells_MCL %<>% SetAllIdent(id="X6_clusters")
table(B_cells_MCL@ident)
X6_clusters_markers <- FindAllMarkers.UMI(B_cells_MCL,logfc.threshold = -Inf,only.pos = FALSE, 
                                          min.pct = 0.01,return.thresh = 1)
X6_clusters_markers <- FindAllMarkers.UMI(B_cells_MCL,logfc.threshold = 0.2,only.pos = T, 
                                          min.pct = 0.5,return.thresh = 0.05)
write.csv(X6_clusters_markers,paste0(path,"X6_clusters_FC0.2_markers.csv"))

X6_clusters_markers = read.csv("output/20190531/X6_clusters_FC0.1_markers.csv",row.names = 1)
#B_cells_MCL %<>% ScaleData()

markers <- c("BTK","CCND1","CCND2","CCND3","CD3D","CD5","CD8A","CDK4","IL2RA",
             "MS4A1","PDCD1","RB1","SOX11")
(MT_gene <- grep("^MT-",X6_clusters_markers$gene))
X6_clusters_markers = X6_clusters_markers[-MT_gene,]
        
g <- DoHeatmap.1(B_cells_MCL, X6_clusters_markers,add.genes = markers, Top_n = 50,
                 use.scaled = T,cex.col = 8,group.cex = 10,
                 title = "Top 50 differentially expressed genes in each B/MCL clusters",
                 group.label.rot = F,cex.row = 2,remove.key =F,title.size = 12)
jpeg(paste0(path,"DE_X6clusters_top50.jpeg"), units="in", width=10, height=7,
     res=600)
print(g)+ theme(strip.text.x = element_text(margin=margin(t = 30, r = 0, b = 0, l = 0)))
dev.off()

g3 <- MakeCorlorBar(B_cells_MCL, X6_clusters_markers,,add.genes = markers, Top_n = 50, do.print = F,
                    do.return = T)
jpeg(paste0(path,"DE_clusters_top50_legend.jpeg"),units="in", width=10, height=7,res=600)
print(g3)
dev.off()

MakeHCorlorBar(B_cells_MCL, group_by = "X6_clusters", remove.legend = F,
              file_name = "DE_clusters_celltype_top50",split.by = "cell.type",do.print = TRUE,do.return=FALSE)


# remove cluster with less than 3 cells======
# no scale down, keep the normal cells.
B_cells_MCL@meta.data$old.ident = B_cells_MCL@meta.data$orig.ident
B_cells_MCL@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",B_cells_MCL@meta.data$orig.ident)
B_cells_MCL@meta.data$X6_orig.ident = paste(B_cells_MCL@meta.data$orig.ident,
                                            B_cells_MCL@meta.data$X6_clusters, sep = "_")
B_cells_MCL@meta.data$X6_orig.ident = gsub('^Normal_.*', 'Normal', B_cells_MCL@meta.data$X6_orig.ident)

table(B_cells_MCL@meta.data$X6_orig.ident)
B_cells_MCL %<>% SetAllIdent(id="X6_orig.ident")
table_B_cells_MCL <- table(B_cells_MCL@meta.data$X6_orig.ident) %>% as.data.frame
(keep.MCL <- table_B_cells_MCL[table_B_cells_MCL$Freq > 2,"Var1"] %>% as.character())
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = keep.MCL)
B_cells_MCL_exp <- AverageExpression(B_cells_MCL)
write.csv(B_cells_MCL_exp,paste0(path,"B_MCL_exp.csv"))

# add color to cluster
gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
gg_color_hue(8)
B_cells_MCL <- AddMetaColor(object = B_cells_MCL, label= "X6_clusters", colors = gg_color_hue(8))

B_cells_MCL %<>% SetAllIdent(id = "orig.ident")
df_samples <- readxl::read_excel("doc/190429_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))

#meta.data = B_cells_MCL@meta.data
#(remove <- colnames(meta.data) %in% "singler1main" %>% which)
#meta.data = meta.data[,-remove]
#B_cells_MCL@meta.data = meta.data
tests <- paste0("test",c(2:10,12,11))
tsnes <- c("X6_clusters")#"cell.type")#,"X6_clusters")
for(tsne in tsnes){
        folder_path <- paste0(path,tsne,"/")
        if(!dir.exists(folder_path)) dir.create(folder_path, recursive = T)
        for(test in tests){
                sample_n = which(df_samples$tests %in% test)
                df <- as.data.frame(df_samples[sample_n,])
                samples <- unique(df$sample)
                rownames(df) = samples
                
                samples <- c(ifelse(length(samples)>5,NA,"Normal"),df$sample[order(df$tsne)])
                print(samples <- samples[!is.na(samples)])
                
                g <- lapply(samples[1:2],function(sample) {
                        SubsetData(B_cells_MCL, ident.use = sample) %>%
                                SetAllIdent(id = tsne) %>%
                                TSNEPlot.1(no.legend = T,do.label =F,label.size=3,size=20,
                                           #colors.use = ExtractMetaColor(.),
                                           do.return = T, label.repel = T,force=2,do.print = F)+
                                ggtitle(sample)+theme(text = element_text(size=15),
                                                      plot.title = element_text(hjust = 0.5))
                })
                jpeg(paste0(folder_path,test,"_Plots_L.jpeg"), units="in", width=10, height=7,
                     res=600)
                print(do.call(cowplot::plot_grid, c(g, nrow = ifelse(length(samples[1:2])>2,2,1))))
                dev.off()
        }
}

###############################
# Doheatmap for Normal / MCL
###############################
df_samples <- readxl::read_excel("doc/190429_scRNAseq_info.xlsx",sheet = "heatmap")
list_samples <- df2list(df_samples)
print(list_samples %>% unlist %>% as.vector %>% unique %in% 
              B_cells_MCL@meta.data$orig.ident)

B_cells_MCL %<>% SetAllIdent(id = "orig.ident")
for(sample in list_samples$MCL){
        subset.MCL <- SubsetData(B_cells_MCL, ident.use = c("Normal",sample))

        # SplitTSNEPlot======
        subset.MCL <- SetAllIdent(subset.MCL,id = "X6_orig.ident")
        SplitTSNEPlot(subset.MCL, select.plots = order(c("Normal",sample)),
                      do.return = FALSE,do.print = TRUE)
        
        # remove cluster with less than 3 cells======
        subset.MCL %<>% SetAllIdent(id = "X6_orig.ident")
        table_subset.MCL <- table(subset.MCL@meta.data$X5_orig.ident) %>% as.data.frame
        keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character()
        subset.MCL <- SubsetData(subset.MCL, ident.use = keep.MCL)
        
        X5_cluster <- subset.MCL@ident %>% unique
        X5_cluster = X5_cluster[-grep("^Normal",X5_cluster)] %>% as.character %>% 
                gsub('.*\\_',"",.) %>% as.numeric %>% sort
        print(ident.1 <- rep("Normal",length(X5_cluster)))
        print(ident.2 <- paste(sample,X5_cluster,sep="_"))
        # FindAllMarkers.UMI======
        subfolder <- paste0(path,"Normal_vs_",sample,"/")
        
        gde.markers <- FindPairMarkers(subset.MCL, ident.1 = c(ident.1,ident.2),
                                       ident.2 = c(ident.2,ident.1), only.pos = T,
                                       logfc.threshold = 1.005,min.cells.group =3,
                                       min.pct = 0.01,return.thresh = 0.05,
                                       save.path = subfolder)
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        subset.MCL <- ScaleData(subset.MCL)
        markers <- c("BTK","CCND1","CCND2","CCND3","CD3D","CD5","CD8A","CDK4","IL2RA",
                     "MS4A1","PDCD1","RB1","SOX11")
        Top_n =25
        g1 <- DoHeatmap.1(subset.MCL, gde.markers, add.genes = markers, Top_n = 25,
                         use.scaled = T,group.order = c("Normal",ident.2),
                         title=paste("Top",Top_n,"DE genes","in Normal B/MCL cells vs.",sample, "B/MCL cells"),
                         group.label.rot = T,cex.row = 3,remove.key =F,title.size = 12)
        jpeg(paste0(subfolder,"DoHeatmap_Normal_",sample,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(g1)
        dev.off()
        #---DoHeatmap vertical bar----
        g2 <- MakeCorlorBar(subset.MCL, gde.markers,Top_n = 25, add.genes = markers, do.print = F,
                            do.return = T)
        jpeg(paste0(subfolder,"DoHeatmap_Normal_",sample,"_legend.jpeg"),units="in", width=10, height=7,res=600)
        print(g2)
        dev.off()
}


# Doheatmap for MCL.1 / MCL.2 ================
df_samples <- readxl::read_excel("doc/190429_scRNAseq_info.xlsx",sheet = "heatmap")
list_samples <- df2list(df_samples)
print(list_samples %>% unlist %>% as.vector %>% unique %in% 
              B_cells_MCL@meta.data$orig.ident)

B_cells_MCL %<>% SetAllIdent(id = "orig.ident")
for(i in 1:length(list_samples$MCL.1)){
        
        samples1 = list_samples$MCL.1[i]
        samples2 = list_samples$MCL.2[i]
        
        subset.MCL <- SubsetData(B_cells_MCL, ident.use = c(samples1,samples2))

        #---SplitTSNEPlot----
        subset.MCL %<>% SetAllIdent(id = "X5_orig.ident")
        select.plots =1:2
        if(all(c(samples1,samples1) != sort(c(samples1,samples2)))){
                select.plots =2:1
        }
        SplitTSNEPlot(subset.MCL, do.return = FALSE,do.print = TRUE,select.plots = select.plots)
        # remove cluster with less than 3 cells======

        table_subset.MCL <- table(subset.MCL@meta.data$X5_orig.ident) %>% as.data.frame
        keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character()
        
        X5_cluster <- subset.MCL@ident %>% unique %>% 
                gsub('.*\\_',"",.) %>% as.numeric %>% sort %>% .[duplicated(.)]
        
        print(ident.1 <- paste(samples1,X5_cluster,sep="_"))
        print(ident.2 <- paste(samples2,X5_cluster,sep="_"))
        
        subset.MCL <- SubsetData(subset.MCL, ident.use = c(ident.1,ident.2))
        
        subfolder <- paste0(path,samples1,"_vs_",samples2,"/")
        gde.markers <- FindPairMarkers(subset.MCL, ident.1 = c(ident.1,ident.2), 
                                       ident.2 = c(ident.2,ident.1), only.pos = T,
                                       logfc.threshold = 0.125,min.cells.group =3,
                                       min.pct = 0.01,
                                       return.thresh = 0.5,
                                       save.path = subfolder)
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()

        #---DoHeatmap.1----
        subset.MCL <- ScaleData(subset.MCL)
        markers <- c("BTK","CCND1","CCND2","CCND3","CD3D","CD5","CD8A","CDK4","IL2RA",
                     "MS4A1","PDCD1","RB1","SOX11")
        Top_n = 25
        g <- DoHeatmap.1(subset.MCL, gde.markers, add.genes = markers, Top_n = 25,
                         group.order = c(ident.1,ident.2),
                         title=paste("Top",Top_n,"DE genes","in",samples1,"vs.",samples2, "in B/MCL cells"),
                         group.label.rot = T,cex.row = 3,remove.key =F,title.size = 12)
        jpeg(paste0(subfolder,"DoHeatmap_",paste(c(samples1,samples2),collapse = "_"),
                    ".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
        #---DoHeatmap vertical bar----
        g1 <- MakeCorlorBar(subset.MCL, gde.markers,Top_n = 25, add.genes = markers, do.print = F,
                            do.return = T)
        jpeg(paste0(subfolder,"DoHeatmap_",paste(c(samples1,samples2),collapse = "_"),
                    "_legend.jpeg"),units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
}

# heatmap.2 for Normal / MCL ================
markers <-  HumanGenes(B_cells_MCL,c("CD19","MS4A1","CD79A","CD5","CD40","CDK4"))
control <- "MD"
tests <- c("test3","test4")
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$samples[sample_n]
        cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% c(control,sample)] #
        subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
        test_markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.1, only.pos = T,
                                           test.use = "MAST")
        top <- test_markers %>% group_by(cluster) %>% top_n(50, avg_logFC)
        y = subset.MCL@scale.data[unique(c(markers,top$gene)),]
        ## Column clustering (adjust here distance/linkage methods to what you need!)
        hc <- hclust(as.dist(1-cor(as.matrix(y), method="spearman")), method="complete")
        cc = gsub("_.*","",hc$labels)
        cc = gsub(control,"#B3DE69",cc)
        cc = gsub(sample,"#195016",cc)
        
        jpeg(paste0(path,"/Heatmap2_",control,"_",sample,".jpeg"), units="in", width=10, height=7,res=600)
        heatmap.2(as.matrix(y),
                  Colv = as.dendrogram(hc), Rowv= FALSE,
                  ColSideColors = cc, trace ="none",labCol = FALSE,dendrogram = "column",#scale="row",
                  key.xlab = "scale log nUMI",
                  cexRow = 0.5,
                  margins = c(2,5),
                  #scale = "row",
                  breaks = seq(-3,3,length.out = 101),
                  col = bluered,
                  main = paste(control, "vs.",sample, "in B cells and MCL cells"))
        par(lend = 1)           # square line ends for the color legend
        legend(0, 0.8,       # location of the legend on the heatmap plot
               legend = c(control, sample), # category labels
               col = c("#B3DE69", "#195016"),  # color key
               lty= 1,             # line style
               lwd = 10            # line width
        )
        dev.off()
}


# heatmap.2 ================
markers <-  HumanGenes(B_cells_MCL,c("CD19","MS4A1","CD79A","CD5","CD40","CDK4"))
tests <- c("test3","test4")
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df_samples[sample_n,] %>% kable() %>% kable_styling()
        samples <- df_samples$samples[sample_n]
        cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
        subset.MCL <- SetAllIdent(subset.MCL,id = "orig.ident")
        test_markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.1, only.pos = T,
                                           test.use = "MAST")
        top <- test_markers %>% group_by(cluster) %>% top_n(30, avg_logFC)
        y = subset.MCL@scale.data[unique(c(markers,top$gene)),]
        ## Column clustering (adjust here distance/linkage methods to what you need!)
        hc <- hclust(as.dist(1-cor(as.matrix(y), method="spearman")), method="complete")
        cc = gsub("_.*","",hc$labels)
        cc = gsub(samples[1],"#B3DE69",cc)
        cc = gsub(samples[3],"#E31A1C",cc)
        cc = gsub(samples[2],"#195016",cc)
        
        jpeg(paste0(path,"/Heatmap2_MD_",samples[2],"_",samples[3],".jpeg"), units="in", width=10, height=7,res=600)
        heatmap.2(as.matrix(y),
                  Colv = as.dendrogram(hc), Rowv= FALSE,
                  ColSideColors = cc, trace ="none",labCol = FALSE,dendrogram = "column",
                  key.xlab = "scale log nUMI",
                  cexRow = 0.5,
                  margins = c(2,5),
                  #scale = "row",
                  breaks = seq(-3,3,length.out = 101),
                  col = bluered,
                  main = paste(samples[1],"vs.",samples[2],"vs.",samples[3], "in B cells"))
        par(lend = 1)           # square line ends for the color legend
        legend(0, 0.8,       # location of the legend on the heatmap plot
               legend = c(samples[1], samples[2], samples[3]), # category labels
               col = c("#B3DE69", "#195016","#E31A1C"),  # color key
               lty= 1,             # line style
               lwd = 10            # line width
        )
        dev.off()
}


# heatmap.2 for MCL bulk RNA ================
markers <-  c("CD19","MS4A1","CD79A","CD5","CD40","CDK4","SOX11","ITGA4","CCND1","CCND2")
markers <-  unique(MCL_markers$gene)[1:1000]
markers <- markers[markers %in% rownames(X181120_MCL_WTS)]

y = X181120_MCL_WTS[markers,]
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(dist(1-cor(as.matrix(y), method="spearman")), method="complete")

jpeg(paste0(path,"/Heatmap2_MCL_bulk_10.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(as.matrix(y),
          Colv = as.dendrogram(hc),
          trace ="none",dendrogram = "both",
          cexRow = 1,
          margins = c(5,5),
          scale = "row",
          col = bluered,
          main = paste("MCL bulk RNA-seq with marker genes"))
dev.off()

# heatmap.2 for MCL scRNA================
MCL_markers_list <- list()
cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% samples]
subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
subset.MCL <- SetAllIdent(subset.MCL,id = "orig.ident")
TSNEPlot(subset.MCL)
MCL_exp <- AverageExpression(subset.MCL)

(samples <- c("DJ","MD","Pt-1294","Pt-RM","Pt-MS","Pt-LM","Pt-1475"))
markers <-  unique(MCL_markers$gene)[1:100]
markers <- markers[markers %in% rownames(MCL_exp)]
y = MCL_exp[markers,]

## Column clustering (adjust here distance/linkage methods to what you need!)
#hc <- hclust(dist(1-cor(as.matrix(y), method="spearman")), method="complete")

jpeg(paste0(path,"/Heatmap2_MCL_sc_100.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(as.matrix(y),
          #Colv = as.dendrogram(hc),
          Colv=FALSE,
          trace ="none",dendrogram = "row",
          cexRow = 0.3,
          margins = c(5,5),
          scale = "row",
          col = bluered,
          main = paste("MCL scRNA RNA top 100 genes"))
dev.off()


# heatmap.2 for MCL bulk RNA + scRNA================
All_MCL <- merge(MCL_exp,log1p(X181120_MCL_WTS), by="row.names")
rownames(All_MCL) <-All_MCL$Row.names
All_MCL <- All_MCL[,-1]
column.sum <- colSums(All_MCL)

jpeg(paste0(path,"/All_MCL.jpeg"), units="in", width=10, height=7,res=600)
par(mfrow=c(2,1))
boxplot(All_MCL, ylab= "log UMI or FPKM")
title(main = "merge MCL bulk and RNA")
plot(x = column.sum)
dev.off()

#markers <-  c("CD19","MS4A1","CD79A","CD5","CD40","CDK4","SOX11","ITGA4","CCND1","CCND2")
markers <-  unique(MCL_markers$gene)[1:1000]
markers <- markers[markers %in% rownames(All_MCL)]

y = All_MCL[markers,]
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(dist(1-cor(as.matrix(y), method="spearman")), method="complete")

jpeg(paste0(path,"/Heatmap2_MCL_sc_bulk_1000.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(as.matrix(y),
          Colv = as.dendrogram(hc),
          trace ="none",dendrogram = "both",
          cexRow = 0.03,
          margins = c(5,5),
          scale = "row",
          col = bluered,
          main = paste("MCL bulk and scRNA top 1000 genes"))
dev.off()

# check HSC===========
B_cells_MCL <- SetAllIdent(B_cells_MCL,id='cell.type')
cell.type_markers <- FindAllMarkers(B_cells_MCL,logfc.threshold=0.5)
PROM1_markers <- FindMarkers.UMI(B_cells_MCL,logfc.threshold = 0, genes.use = "PROM1",
                                 min.pct = 0,ident.1 = "MCL/HSC",return.thresh = 1)

write.table(cell.type_markers,paste0(path,"cell.type_markers.txt"))


object@meta.data$cell.type = gsub("B_cells.*","B_cells",object@meta.data$singler1sub)
object@meta.data$cell.type = gsub("T_cells.*","T_cells",object@meta.data$cell.type)
object@meta.data$cell.type = gsub("MPP|MEP|CLP|HSC|CMP|GMP","MCL.HSC",object@meta.data$cell.type)

table_subset.MCL <- table(object@meta.data$cell.type) %>% as.data.frame 
(keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character())
object <- SetAllIdent(object, id="cell.type")
table(object@ident)
major_cell <- SubsetData(object,ident.use=keep.MCL)
major_cell_exp <- AverageExpression(major_cell)
write.csv(major_cell_raw_exp, paste0(path,"major_cell_raw_exp.csv"))



PROM1_markers <- FindMarkers.UMI(object,logfc.threshold = 0, genes.use = "PROM1",
                                 min.pct = 0,ident.1 = "MCL.HSC",return.thresh = 1)
