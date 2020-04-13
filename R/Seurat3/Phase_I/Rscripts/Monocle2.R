invisible(lapply(c("Seurat","monocle","dplyr",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- "output/202000413_monocle2/"
if(!dir.exists(path)) dir.create(path, recursive = T)
#SBATCH --mem=32G
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))

# load data
object = readRDS(file = "data/MCL_41_B_20200225.rds")
Idents(object) = "orig.ident"
object %<>% AddMetaColor(label= "X4clusters", colors = gg_color_hue(4))
object$X4_orig.ident = paste(object$orig.ident,
                             object$X4clusters, sep = "_")
sample_pairs = list(c("N01","N02","N03","N04"),
                    c("PtU01","PtU02","PtU03","PtU04"),
                    c("Pt11_LN1","Pt11_1","Pt11_14","Pt11_28"),
                    c("Pt11_LN1","Pt11_1"),
                    c("Pt11_1","Pt11_14","Pt11_28"),
                    c("Pt11_14","Pt11_28"),
                    c("Pt17_LN1","Pt17_2","Pt17_7","Pt17_12","Pt17_31"),
                    c("Pt17_LN1","Pt17_2"),
                    c("Pt17_2","Pt17_7"),
                    c("Pt17_2","Pt17_12"),
                    c("Pt17_7","Pt17_12"),
                    c("Pt17_7","Pt17_12","Pt17_31"),
                    c("Pt25_SB1","Pt25_1","Pt25_1_8","Pt25_24","Pt25_25Pd","Pt25_AMB25Pd"),
                    c("Pt25_SB1","Pt25_1"),
                    c("Pt25_1","Pt25_1_8"),
                    c("Pt25_1_8","Pt25_24"),
                    c("Pt25_24","Pt25_25Pd"),
                    c("Pt25_25Pd","Pt25_AMB25Pd"),
                    c("Pt25_SB1","Pt25_AMB25Pd"),
                    c("Pt27_1","Pt27_1_8","Pt27_12"),
                    c("Pt27_1","Pt27_1_8"),
                    c("Pt27_1","Pt27_12"),
                    c("Pt27_1_8","Pt27_12"),
                    c("PtB13_Ibp","PtB13_Ib1","PtB13_IbR"),
                    c("PtB13_Ibp","PtB13_Ib1"),
                    c("PtB13_Ib1","PtB13_IbR"),
                    c("PtB13_Ibp","PtB13_IbR"),
                    c("Pt28_LN1","Pt28_1","Pt28_4","Pt28_28"),
                    c("Pt28_LN1","Pt28_1"),
                    c("Pt28_1","Pt28_4"),
                    c("Pt28_4","Pt28_28"),
                    c("Pt28_LN1","Pt28_28"))
(sample = sample_pairs[[i]])
object %<>% subset(idents = sample)

save.path = paste0(path, paste(sample, collapse = "-"),"/")
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

Idents(object) = "X4clusters"
X4clusters_color = ExtractMetaColor(object)
# Store Data in a CellDataSet Object
pd <- new("AnnotatedDataFrame", data = object@meta.data)
fd <- data.frame(gene_short_name = rownames(object),
                 row.names = rownames(object))
fd <- new("AnnotatedDataFrame", data = fd)
expr_matrix <- as.matrix(object@assays$SCT@counts)
cds <- newCellDataSet(expr_matrix, phenoData = pd,featureData = fd,
                      expressionFamily = negbinomial())

#Estimate size factors and dispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
table(pData(cds)$cell.types)
#Filtering low-quality cells
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))
length(expressed_genes)

#Trajectory step 1: choose genes that define a cell's progress
clustering_DEG_genes <- differentialGeneTest(cds[expressed_genes,],
                                             fullModelFormulaStr = "~X4_orig.ident",
                                             cores = detectCores()/2)
saveRDS(clustering_DEG_genes, paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_DE.rds"))

#clustering_DEG_genes = readRDS(paste0(path,sample,"/monocle2_",paste(sample, collapse = "-"),"_DE.rds"))
cds_ordering_genes <-row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][
        1:min(nrow(clustering_DEG_genes),1000)]
cds %<>% setOrderingFilter(ordering_genes = cds_ordering_genes)
cds %<>% reduceDimension(method = 'DDRTree')
cds %<>% orderCells()
group_by <- c("X4_orig.ident", "orig.ident","X4clusters", "Pseudotime")
for(i in seq_along(color_by)){
        jpeg(paste0(save.path,"trajectory_",paste(sample, collapse = "-"),"_",group_by[i],".jpeg"),
             units="in", width=7, height=7,res=600)
        g <- plot_cell_trajectory(cds, color_by = group_by[i],show_branch_points = FALSE)
        if(color_by[i] == "X4clusters") g = g + scale_color_manual(values=X4clusters_color)
        print(g)
        dev.off()
}
