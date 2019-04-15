########################################################################
#
#  0 setup environment, install libraries if nLynchessary, load libraries
# 
# ######################################################################

invisible(lapply(c("R.utils","Seurat","magrittr","dplyr","harmony"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data/")) dir.create("data")

########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
df_samples <- readxl::read_excel("doc/190406_scRNAseq_info.xlsx")
df_samples
colnames(df_samples) = tolower(colnames(df_samples))
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:10)))
df_samples = df_samples[sample_n,]
attach(df_samples)
samples = df_samples$sample
(new_samples <- samples[as.character(df_samples$date) %in% "2019-04-06"])

args <- commandArgs(trailingOnly = TRUE)
args[1] = as.integer(args[1])
new_sample <- new_samples[args[1]]
print(paste("add new sample:",new_sample))
sample_pth <- paste0(path,"/",new_sample)
if(!dir.exists(sample_pth)) dir.create(sample_pth, recursive = T)


(load( file = "data/MCL_raw_36_20190412.Rda"))
object %<>% SetAllIdent(id="orig.ident")
object <- SubsetData(object, ident.remove = new_sampels[-args[1]])
#======1.5 FindVariableGenes=======================
object <- NormalizeData(object = object)
jpeg(paste0(sample_pth,new_sample,"S1_dispersion.jpeg"), units="in", width=10, height=7,res=600)
object <- FindVariableGenes(object = object, mean.function = ExpMean, 
                            dispersion.function = LogVMR, do.plot = T, 
                            x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.5)
dev.off()
length(object@var.genes)

#======1.6 PCA =========================
object %<>% ScaleData
object %<>% RunPCA(pc.genes = object@var.genes, pcs.compute = 100, do.print = F)

jpeg(paste0(sample_pth,"S1_PCElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
PCElbowPlot(object, num.pc = 100)
dev.off()
GC()
pcs =1:75
#======1.6 RunHarmony=======================
jpeg(paste0(sample_pth,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony("orig.ident", dims.use = pcs,
                                theta = 2, plot_convergence = TRUE,
                                nclust = 50, max.iter.cluster = 100))
dev.off()

object@ident %<>% factor(levels = samples)
p1 <- DimPlot(object = object, reduction.use = "harmony", pt.size = 0.3, no.legend = T,
              group.by = "orig.ident", do.return = T)
p2 <- VlnPlot(object = object, features.plot = "Harmony1", 
              group.by = "orig.ident", do.return = TRUE, x.lab.rot = T)
jpeg(paste0(sample_pth,"S1_Harmony_vplot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1,p2)
dev.off()

jpeg(paste0(sample_pth,"S1_Harmony_DimHeatmap.jpeg"), units="in", width=10, height=7,res=600)
DimHeatmap(object = object, reduction.type = "harmony", cells.use = 500, 
           dim.use =c(1:3, 48:50,73:75), do.balanced = TRUE)
dev.off()

#========1.6 Seurat tSNE Functions for Integrated Analysis Using Harmony Results=======
system.time(
        object %<>% RunTSNE(reduction.use = "harmony", dims.use = pcs, do.fast = TRUE))
system.time(
        object %<>% FindClusters(reduction.type = "harmony", resolution = 0.6, dims.use = pcs,
                              save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                              force.recalc = TRUE, print.output = FALSE))


p3 <- TSNEPlot(object, do.return = T, pt.size = 0.3, group.by = "orig.ident", no.legend = T)
p4 <- TSNEPlot(object, do.label = T, do.return = T, pt.size = 0.3, no.legend = T)

jpeg(paste0(sample_pth,"S1_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3+ggtitle("group by samples")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p4+ggtitle("group by clusters")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")))
dev.off()

g_Harmony <- TSNEPlot.1(object = object, do.label = T, group.by = "ident",
                        do.return = TRUE, no.legend = T, 
                        #colors.use = ExtractMetaColor(object),
                        pt.size = 1,label.size = 6 )+
        ggtitle("Tsne plot of all clusters")+
        theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(sample_pth,"TSNEplot-Harmony.jpeg"), units="in", width=10, height=7,res=600)
print(g_Harmony)
dev.off()

#saveRDS(object@scale.data, file = "data/MCL.scale.data_Harmony_36_20190412.rds")
#object@scale.data = NULL; GC()
#save(object, file = "data/MCL_Harmony_36_20190413.Rda")
