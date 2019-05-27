#check package
library("devtools")
library("DDRTree")
library("pheatmap")
library("monocle")
library("reshape")
source("../R/Seurat_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 5.1 Importing data from Seurat object=================
(load(file="data/B_cells_MCL_43_20190521.Rda"))
object <- B_cells_MCL
table(object@ident)
object_Mono <- importCDS(object, import_all = TRUE)

# 5.1.1 Estimate size factors and dispersions
# estimateSizeFactors() and estimateDispersions() will only work,
# and are only needed, if you are working with a CellDataSet 
# with a negbinomial() or negbinomial.size() expression family.
object_Mono <- estimateSizeFactors(object_Mono)
object_Mono <- estimateDispersions(object_Mono)

# 5.1.2 Filtering low-quality cells
object_Mono <- detectGenes(object_Mono, min_expr = 0.1)
print(head(fData(object_Mono)))
print(head(pData(object_Mono)))

# 5.1.3 If you are using RPC values to measure expresion, 
# as we are in this vignette, it's also good to look at the distribution
# of mRNA totals across the cells:
pData(object_Mono)$Total_mRNAs <- Matrix::colSums(exprs(object_Mono))
upper_bound <- 10^(mean(log10(pData(object_Mono)$Total_mRNAs)) + 
                     2*sd(log10(pData(object_Mono)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(object_Mono)$Total_mRNAs)) - 
                     2*sd(log10(pData(object_Mono)$Total_mRNAs)))
qplot(Total_mRNAs, data = pData(object_Mono), color = conditions, geom = "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

# 5.2 Classifying and counting cells of different types
# 5.2.1 Classifying cells with manualHierarchy
object <- SetAllIdent(object, id = "X6_cluster")
all(row.names(pData(object_Mono)) == names(object@ident))
pData(object_Mono)$X6_cluster <- object@ident
table(pData(object_Mono)$X6_cluster)

jpeg(paste0(path,"Pie-cluster.jpeg"), units="in", width=10, height=7,res=600)
pie <- ggplot(pData(object_Mono), aes(x = factor(1), fill = factor(X6_cluster))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()
# 5.3 Constructing Single Cell Trajectories
# 5.3.1 choosing genes that define progress
expressed_genes <- row.names(subset(fData(object_Mono), num_cells_expressed >= 10))
length(expressed_genes)
diff_test_res <- differentialGeneTest(object_Mono[expressed_genes,],
                                      fullModelFormulaStr = "~ X6_clusters",
                                      cores = 1) #takes long time

# get DE genes from Seurat
markers.B = read.csv(paste0("output/20190522/X6_clusters_markers.csv"))

# a. All markers identified by Seurat for object clusters (in B_cells_MCL-only analysis; let’s call it list A)
ordering_genes = unique(markers.B[markers.B$p_val_adj<0.1^10,"gene"])
#b. A subset of most significant differentially expressed markers in the list A: p val adj <10(-10)
ordering_genes = unique(markers.B[markers.B$p_val_adj<0.1^100,"gene"])
#c. Transcription factors (TFs) in the list A
#ordering_genes = FilterGenes(B_cells_MCL,TF)

length(ordering_genes)
object_Mono <- setOrderingFilter(object_Mono, ordering_genes)
jpeg(paste0(path,"a-plot_ordering_genes.jpeg"),units="in", width=10, height=7,res=600)
plot_ordering_genes(object_Mono)
dev.off()

# 5.3.2 reduce data dimensionality
#Now we're ready to try clustering the cells:.
jpeg(paste0(path,"a-plot_pc_variance_explained.jpeg"),units="in", width=10, height=7,res=600)
plot_pc_variance_explained(object_Mono, return_all = F) # norm_method = 'log',
dev.off()
object_Mono <- reduceDimension(object_Mono, max_components = 2,
                                  method = 'DDRTree')
#Trajectory step 3: order cells along the trajectory
(load(file = "data/B_Mono_43_20190523.Rda"))
object_Mono <- orderCells(object_Mono)

g1 <- plot_cell_trajectory(object_Mono, color_by = "X6_clusters",cell_size = 3)
g2 <- plot_cell_trajectory(object_Mono, color_by = "State",cell_size = 3)

jpeg(paste0(path,"a-trajectory_B_cells_MCL.jpeg"), units="in", width=10, height=7,res=600)
g1
dev.off()

jpeg(paste0(path,"a-trajectory_B_cells_MCL_state.jpeg"), units="in", width=10, height=7,res=600)
g2
dev.off()

save(object_Mono, file = "data/B_Mono_43_20190523.Rda")
