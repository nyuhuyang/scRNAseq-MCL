invisible(lapply(c("Seurat","monocle"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))

# load data
object = readRDS(file = "data/MCL_41_B_20200225.rds")
Idents(object) = "orig.ident"
samples = unique(Idents(object)) %>% as.character()
(sample = samples[i])
object %<>% subset(idents = sample)

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
                                             fullModelFormulaStr = "~X4clusters",
                                             cores = detectCores()/2)
cds_ordering_genes <-row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][
        1:min(nrow(clustering_DEG_genes),1000)]
cds %<>% setOrderingFilter(ordering_genes = cds_ordering_genes)
cds %<>% reduceDimension(method = 'DDRTree')
cds %<>% orderCells()
color_by <- c("X4clusters", "Pseudotime")
if(!dir.exists(paste0(path, sample))) dir.create(paste0(path, sample), recursive = T)
for(i in seq_along(color_by)){
        jpeg(paste0(path, sample,"/trajectory_",sample,"_",color_by[i],".jpeg"), units="in", width=7, height=7,res=600)
        print(plot_cell_trajectory(cds, color_by = "X4clusters"))
        dev.off()
}