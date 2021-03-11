########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#SBATCH --mem=128G  # memory requested, units available: K,M,G,T

library(Seurat)
library(dplyr)
library(MAST)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

# load data
(load(file = "data/B_cells_MCL_43_20190917.Rda"))

df_samples <- readxl::read_excel("doc/191030_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:12)))
df_samples = df_samples[sample_n,]

B_cells_MCL@meta.data$orig.ident %<>% plyr::mapvalues(from = unique(df_samples$sample),
                                                 to = unique(df_samples$publication.id))
table(B_cells_MCL@meta.data$orig.ident)
NewNames = paste0(B_cells_MCL@meta.data$orig.ident,"_",B_cells_MCL@meta.data$Barcode)
B_cells_MCL %<>% RenameCells(new.names = NewNames)
rownames(B_cells_MCL@reductions$tsne@cell.embeddings) = colnames(B_cells_MCL)

Idents(B_cells_MCL) = "groups"
B_cells_MCL %<>% subset(idents = "Pt-25")
##############
# DE genes between Clusters 5 cluster top 40
###############
Idents(B_cells_MCL) <- "X5_clusters"
B_cells_MCL %<>% subset(idents = args)
(category = rev(c("Pt25_SB1", "Pt25_1","Pt25_1_8","Pt25_24","Pt25_25","Pt25_AMB25")))
B_cells_MCL@meta.data$orig.ident %<>% factor(levels = category)
Idents(B_cells_MCL) = "orig.ident"
table(Idents(B_cells_MCL))
X5_clusters_markers <- FindAllMarkers.UMI(B_cells_MCL,logfc.threshold = 0.01,only.pos = FALSE, 
                                          min.pct = 0.1,return.thresh = 1)
write.csv(X5_clusters_markers,paste0(path,"Pt-25_Clusters_",args,"_FC0.01_markers.csv"))

X5_clusters_markers = read.csv(paste0("Yang/B_longitudinal_heatmap/DE_analysis_files/Pt-25_Clusters_",
                                      args,"_FC0.01_markers.csv"),
                               row.names = 1, stringsAsFactors = F)

DT = data.frame(category = category,
                int = 1:length(category),
                stringsAsFactors = F)
X5_clusters_markers = left_join(X5_clusters_markers, DT, by = c("cluster" = "category"))
X5_clusters_markers = X5_clusters_markers[order(X5_clusters_markers$int),]
markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
(MT_gene <- grep("^MT-",X5_clusters_markers$gene))
X5_clusters_markers = X5_clusters_markers[-MT_gene,]
Top_n = 40
top = X5_clusters_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
features = c(as.character(top$gene),
             tail(VariableFeatures(object = B_cells_MCL), 2),
             markers)
B_cells_MCL %<>% ScaleData(features=features)
featuresNum <- make.unique(features, sep = ".")
B_cells_MCL = MakeUniqueGenes(object = B_cells_MCL, features = features)
DoHeatmap.1(B_cells_MCL, features = featuresNum, Top_n = Top_n,
            do.print=T, angle = 0, group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0, label=F, cex.row= 2, legend.size = 0,width=10, height=6.5,
            pal_gsea = FALSE,
            title = paste("Top",Top_n,"DE genes in Longtitudinal Pt25 in cluster",args))
