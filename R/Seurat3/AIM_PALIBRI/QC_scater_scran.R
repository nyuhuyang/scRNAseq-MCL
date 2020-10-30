########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0
invisible(lapply(c("R.utils","Seurat","dplyr","kableExtra","ggplot2","scater",
                   "scran","BiocSingular","Matrix","cowplot"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("doc")) dir.create("doc")
save.path = "Yang/PALIBR_phasI_AIM/QC/"
########################################################################
#
#  1 Data preprocessing
#
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20200713_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower

#keep = c("Normal","Untreated","Pt-113","AIM13","AIM17","AIM24")
#df_samples = df_samples[df_samples$patient %in% keep,]
(attach(df_samples))
samples = df_samples$sample
df_samples$species = df_samples$organism
df_samples %>% kable %>% kable_styling()

# check missing data
current <- list.files("data/scRNAseq")
(current <- current[!grepl(".Rda|RData",current)])
(missing_data <- df_samples$sample.id[!(df_samples$sample.id %in% current)])

message("Loading the datasets")
## Load the dataset
Seurat_raw <- list()
Seurat_list <- list()
for(i in seq_along(df_samples$sample)){
        Seurat_raw[[i]] <- Read10X(data.dir = paste0("data/scRNAseq/",
                                                     df_samples$sample.id[i],"/",
                                                     df_samples$read.path[i]))
        colnames(Seurat_raw[[i]]) = paste0(df_samples$sample[i],"_",colnames(Seurat_raw[[i]]))
        Seurat_list[[i]] <- CreateSeuratObject(Seurat_raw[[i]],
                                               project = df_samples$`sample name`[i],
                                                 min.cells = 0,
                                               min.features = 0)
        Seurat_list[[i]]@meta.data$tests <- df_samples$tests[i]
        Progress(i, length(df_samples$sample))
}
remove(Seurat_raw);GC()

#========1.1.3 g1 QC plots before filteration=================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), Seurat_list)
remove(Seurat_list);GC()

# read and select mitochondial genes
mito = "^MT-"
message("mito.genes:")

(mito.features <- grep(pattern = mito, x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = mito)
Idents(object) = "orig.ident"
Idents(object) %<>% factor(levels = df_samples$`sample name`)
g1 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 1, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=10,angle = 90),legend.position="none")
})
save(g1,file= paste0(path,"g1","_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))
#load(paste0(save.path,"g1_58_20201009.Rda"))
#============1.2 scatter ======================
Seurat_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()

for(i in 1:length(df_samples$sample)){
        high.mito <- isOutlier(Seurat_list[[i]]$percent.mt, nmads=3, type="higher")
        low.lib <- isOutlier(log10(Seurat_list[[i]]$nCount_RNA), type="lower", nmad=3)
        low.genes <- isOutlier(log10(Seurat_list[[i]]$nFeature_RNA), type="lower", nmad=3)
        discard <- high.mito | low.lib | low.genes
        print(data.frame(HighMito= sum(high.mito),LowLib=sum(low.lib),
                         LowNgenes=sum(low.genes),Discard=sum(discard)))
        Seurat_list[[i]] <- Seurat_list[[i]][,!discard]
        #print(summary(!discard))
        print(i)
}

object <- Reduce(function(x, y) merge(x, y, do.normalize = F), Seurat_list)
remove(Seurat_list);GC()

object %<>% subset(subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 50)
# FilterCellsgenerate Vlnplot before and after filteration
Idents(object) = "orig.ident"
Idents(object) %<>% factor(levels = df_samples$`sample name`)

g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 1, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=10,angle = 90),legend.position="none")
})
save(g2,file= paste0(path,"g2","_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))
#load(paste0(save.path,"g2_58_20201009.Rda"))

jpeg(paste0(path,"S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nFeature_RNA before filteration")+
                        scale_y_log10(limits = c(100,10000))+
                        theme(axis.text.x = element_text(size=5),
                              plot.title = element_text(hjust = 0.5)),
                g2[[1]]+ggtitle("nFeature_RNA after filteration")+
                        scale_y_log10(limits = c(100,10000))+
                        theme(axis.text.x = element_text(size=5),
                              plot.title = element_text(hjust = 0.5))))
dev.off()
jpeg(paste0(path,"S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nCount_RNA before filteration")+
                        scale_y_log10(limits = c(500,100000))+
                        theme(axis.text.x = element_text(size=5),
                              plot.title = element_text(hjust = 0.5)),
                g2[[2]]+ggtitle("nCount_RNA after filteration")+
                        scale_y_log10(limits = c(500,100000))+
                        theme(axis.text.x = element_text(size=5),
                              plot.title = element_text(hjust = 0.5))))
dev.off()
jpeg(paste0(path,"S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % before filteration")+
                        ylim(c(0,50))+
                        theme(axis.text.x = element_text(size=5),
                              plot.title = element_text(hjust = 0.5)),
                g2[[3]]+ggtitle("mito % after filteration")+
                        ylim(c(0,50))+
                        theme(axis.text.x = element_text(size=5),
                              plot.title = element_text(hjust = 0.5))))
dev.off()

#====
format(object.size(object),unit = "GB")
save(object, file = "data/MCL_AIM_58_20201009.Rda")
