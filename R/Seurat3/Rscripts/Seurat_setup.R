########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#devtools::install_github(repo = "ChristophH/sctransform", ref = "develop")
invisible(lapply(c("Seurat","dplyr","kableExtra","ggplot2","harmony","cowplot","sctransform"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("doc")) dir.create("doc")

########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/190406_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:12)))
df_samples = df_samples[sample_n,]
attach(df_samples)
df_samples
samples = sample

#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_36_20190410.Rda"))
names(sce_list)
object_list <- lapply(sce_list, as.Seurat)

for(i in 1:length(samples)){
        object_list[[i]]@meta.data$tests <- df_samples$tests[i]
        object_list[[i]]@meta.data$conditions <- df_samples$conditions[i]
        object_list[[i]]@meta.data$projects <- df_samples$project[i]
        object_list[[i]]@meta.data$groups <- df_samples$group[i]
        object_list[[i]]@meta.data$tissues <- df_samples$tissue[i]
        object_list[[i]]@meta.data$tsne <- df_samples$tsne[i]
        
}
#========1.3 merge ===================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), object_list)
remove(sce_list,object_list);GC()
save(object, file = "data/MCL_Harmony_36_20190411.Rda")

#======1.4 mito, QC, filteration =========================
# store mitochondrial percentage in object meta data
object <- PercentageFeatureSet(object = object, pattern = "^MT-", col.name = "percent.mt")
Idents(object) = factor(Idents(object),levels = samples)
(load(file = "output/20190410/g1_36_20190410.Rda"))

(remove <- which(colnames(object@meta.data) %in%c("is_cell_control",
                                                  "pct_counts_in_top_500_features_Mito")))
meta.data = object@meta.data[,-seq(remove[1], remove[2], by=1)]
object@meta.data = meta.data 

object %<>% subset(subset = nFeature_RNA > 500 & nCount_RNA > 800 & percent.mt < 50)
# FilterCellsgenerate Vlnplot before and after filteration
g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
    VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=8),legend.position="none")
})

save(g2,file= paste0(path,"g2_36_20190410.Rda"))
jpeg(paste0(path,"S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nFeature_RNA before filteration")+
                    scale_y_log10(limits = c(100,10000)),
                g2[[1]]+ggtitle("nFeature_RNA after filteration")+
                    scale_y_log10(limits = c(100,10000))))
dev.off()
jpeg(paste0(path,"S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nCount_RNA before filteration")+
                    scale_y_log10(limits = c(500,100000)),
                g2[[2]]+ggtitle("nCount_RNA after filteration")+ 
                    scale_y_log10(limits = c(500,100000))))
dev.off()
jpeg(paste0(path,"S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % before filteration")+
                    ylim(c(0,50)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,50))))
dev.off()

#======1.5 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../R/seurat_resources/regev_lab_cell_cycle_genes.txt")
s.genes <- FilterGenes(object=object ,cc.genes[1:43])
g2m.genes <- FilterGenes(object,cc.genes[44:97])
object <- CellCycleScoring(object = object, s.features = s.genes, g2m.features = g2m.genes)

#======1.6 NormalizeData and ScaleData =========================
# run sctransform============
object <- SCTransform(object = object, vars.to.regress = "orig.ident", verbose = T)


#======1.6 PCA Determine statistically significant principal components=======================
# Run the standard workflow for visualization and clustering
object <- RunPCA(object, features = VariableFeatures(object), npcs = 100, verbose = F)

jpeg(paste0(path,"S1_PCElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 100)
dev.off()
jpeg(paste0(path,"S1_PCHeatmap.jpeg"), units="in", width=10, height=7,res=600)
DimHeatmap(object, dims = c(1:3,38:40,48:50), cells = 500, balanced = TRUE)
dev.off()

#object <- JackStraw(object, num.replicate = 100,dims = 50)
#object <- ScoreJackStraw(object, dims = 1:50)
#p4 <- JackStrawPlot(object = object, dims = 30:40)
#jpeg(paste0(path,"/S1_JackStrawPlot.jpeg"), units="in", width=10, height=7,res=600)
#p4
#dev.off()

dim = 50
object <- FindNeighbors(object, reduction = "pca", dims = 1:dim)
object <- FindClusters(object, reduction = "pca", resolution = 0.8)

#======1.7 RunHarmony=======================
object <- RunTSNE(object, reduction = "pca",dims = 1:dim)
p0 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = T,
                 label.size = 4, repel = T)+ NoLegend()+
        ggtitle("Clustering without harmonization")+
        theme(plot.title = element_text(hjust = 0.5))
system.time(object <- RunHarmony.1(object, group.by.vars= "orig.ident", dims.use = 1:dim,
                                   theta = 0, plot_convergence = TRUE,
                                   nclust = 100, max.iter.cluster = 100))

#========1.6 Seurat tSNE Functions for Integrated Analysis Using Harmony Results=======
object <- FindNeighbors(object, reduction = "harmony", dims = 1:dim)
object <- FindClusters(object, reduction = "harmony", resolution = 0.8)
object@reductions$tsne = NULL
object <- RunTSNE(object = object,reduction = "harmony", dims = 1:dim, verbose = FALSE)

p1 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = T,
                 label.size = 4, repel = T)+ NoLegend()

p2 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T)+ NoLegend()
p3 <- TSNEPlot.1(object, group.by="RNA_snn_res.0.8",pt.size = 1,label = T,
                 label.size = 4, repel = T)+ NoLegend()

jpeg(paste0(path,"S1_pca_vs_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p0+ggtitle("Clustering without harmonization")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18)),
          p1+ggtitle("Clustering with harmonization")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()

jpeg(paste0(path,"S1_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p2+ggtitle("group by samples")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18)),
          p3+ggtitle("group by clusters")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()

g_Harmony <- TSNEPlot.1(object = object, label = T, group.by = "ident",
                        label.size = 6 )+ NoLegend()+
        ggtitle("TSNE plot of all clusters")+
        theme(plot.title = element_text(hjust = 0.5,size = 18))

jpeg(paste0(path,"TSNEplot-Harmony.jpeg"), units="in", width=10, height=7,res=600)
print(g_Harmony)
dev.off()
save(object, file = "data/MCL_full_Harmony_36_20190412.Rda")
saveRDS(object@assays$RNA@scale.data, file = "data/MCL.scale.data_36_20190412.Rds")
object@assays$RNA@scale.data =matrix()
save(object, file = "data/MCL_Harmony_36_20190412.Rda")
