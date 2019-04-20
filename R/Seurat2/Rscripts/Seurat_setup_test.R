
########################################################################
#
#  0 setup environment, install libraries if nLynchessary, load libraries
# 
# ######################################################################
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("R.utils","Seurat","magrittr","dplyr","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data/")) dir.create("data")
set.seed=101
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
args <- slurm_arrayid
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/190406_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:10)))
df_samples = df_samples[sample_n,]
samples = df_samples$sample
(new_samples <- samples[as.character(df_samples$date) %in% "2019-04-06"])

args[1] = as.integer(args[1])
new_sample1 <- new_samples[c(args[1],4:6)]
print(paste("add new sample:",paste(new_sample1[1],",Pt-27")))
sample_pth <- paste0(path,new_sample1[1],"_Pt-27")
if(!dir.exists(sample_pth)) dir.create(sample_pth, recursive = T)

rm_samples <- new_samples[!(new_samples %in% new_sample1)]
print(rm_samples)
print(paste("rm new sample:",paste(rm_samples,collapse = ",")))

#======1.2 load  data =========================
(load(file = "data/MCL_Harmony_36_20190420.Rda"))
object %<>% SetAllIdent(id="orig.ident")
object <- SubsetData(object, ident.remove = (rm_samples))

#======1.5 FindVariableGenes=======================
object <- NormalizeData(object = object)
jpeg(paste0(sample_pth,"S1_dispersion.jpeg"), units="in", width=10, height=7,res=600)
object <- FindVariableGenes(object = object, mean.function = ExpMean, 
                            dispersion.function = LogVMR, do.plot = F, 
                            x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.5)
dev.off()
length(object@var.genes)

#======1.6 PCA =========================
object %<>% ScaleData
object %<>% RunPCA(pc.genes = object@var.genes, pcs.compute = 75, do.print = F)

GC()
pcs =1:75
#======1.6 RunHarmony=======================
jpeg(paste0(sample_pth,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object <- RunHarmony(object, "orig.ident", dims.use = pcs,
                                   theta = 2, plot_convergence = TRUE,
                                   nclust = 50, max.iter.cluster = 100))
dev.off()

#========1.6 Seurat tSNE Functions for Integrated Analysis Using Harmony Results=======
system.time(
        object %<>% RunTSNE(reduction.use = "harmony", dims.use = pcs, do.fast = TRUE))
system.time(
        object %<>% FindClusters(reduction.type = "harmony", resolution = 0.6, dims.use = pcs,
                                 save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                                 force.recalc = TRUE, print.output = FALSE))


p3 <- TSNEPlot(object, do.return = T, pt.size = 0.3, group.by = "orig.ident")
p4 <- TSNEPlot(object, do.label = T, do.return = T, pt.size = 0.3)

jpeg(paste0(sample_pth,"S1_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3+ggtitle("group by samples")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p4+ggtitle("group by clusters")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")))
dev.off()

g_Harmony <- TSNEPlot.1(object = object, do.label = T, group.by = "ident",
                        do.return = TRUE, no.legend = F, 
                        #colors.use = ExtractMetaColor(object),
                        pt.size = 1,label.size = 6 )+
        ggtitle("Tsne plot of all clusters")+
        theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(sample_pth,"TSNEplot-Harmony.jpeg"), units="in", width=10, height=7,res=600)
print(g_Harmony)
dev.off()

#object@scale.data = NULL; GC()
#save(object, file = "data/MCL_Harmony_30_20190320~.Rda")
