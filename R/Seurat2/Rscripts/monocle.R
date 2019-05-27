invisible(sapply(c("DDRTree","pheatmap","monocle","magrittr","reshape"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
source("R/util.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# 5.0 load group number =================
args <- commandArgs(trailingOnly = TRUE)
args[1] <- as.numeric(args[1])
groups <- c("AFT-03","AFT-04","Pt-11","Pt-13","Pt-17",
            "Pt-AA13","Pt-25", "Pt-27")
# 5.1 Importing data from Seurat object=================
(load(file="data/B_cells_MCL_43_20190521.Rda"))
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="groups")
object <- SubsetData(B_cells_MCL, ident.use = groups[args[1]])
object <- SetAllIdent(object, id="orig.ident")
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

# 5.2 Classifying and counting cells of different types
# 5.2.1 Classifying cells with manualHierarchy
object <- SetAllIdent(object, id = "orig.ident")
all(row.names(pData(object_Mono)) == names(object@ident))
#pData(object_Mono)$X8_cluster <- object@ident
table(pData(object_Mono)$orig.ident)

jpeg(paste0(path,groups[args[1]],"_Pie-cluster.jpeg"), units="in", width=10, height=7,res=600)
pie <- ggplot(pData(object_Mono), aes(x = factor(1), fill = factor(orig.ident))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()
# 5.3 Constructing Single Cell Trajectories
# 5.3.1 choosing genes that define progress
expressed_genes <- row.names(subset(fData(object_Mono), num_cells_expressed >= 10))
length(expressed_genes)
diff_test_res <- differentialGeneTest(object_Mono[expressed_genes,],
                                      fullModelFormulaStr = "~ orig.ident",
                                      cores = 1) #takes long time
write.csv(diff_test_res, file = paste0(path,groups[args[1]],"_diff_test_res.csv"))
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
length(ordering_genes)
object_Mono <- setOrderingFilter(object_Mono, ordering_genes)
jpeg(paste0(path,groups[args[1]],"_plot_ordering_genes.jpeg"),units="in", width=10, height=7,res=600)
plot_ordering_genes(object_Mono)
dev.off()

# 5.3.2 reduce data dimensionality
#Now we're ready to try clustering the cells:.
jpeg(paste0(path,groups[args[1]],"_plot_pc_variance_explained.jpeg"),units="in", width=10, height=7,res=600)
plot_pc_variance_explained(object_Mono, return_all = F) # norm_method = 'log',
dev.off()
object_Mono <- reduceDimension(object_Mono, max_components = 2,
                                  method = 'DDRTree')
#Trajectory step 3: order cells along the trajectory
object_Mono <- orderCells(object_Mono)

g1 <- plot_cell_trajectory(object_Mono, color_by = "orig.ident",cell_size = 3)
g2 <- plot_cell_trajectory(object_Mono, color_by = "State",cell_size = 3)

jpeg(paste0(path,groups[args[1]],"_trajectory_B_cells_MCL.jpeg"), units="in", width=10, height=7,res=600)
g1
dev.off()

jpeg(paste0(path,groups[args[1]],"_trajectory_B_cells_MCL_state.jpeg"), units="in", width=10, height=7,res=600)
g2
dev.off()

save(object_Mono, file = paste0("output/",groups[args[1]],"_B_Mono_",gsub("-","",Sys.Date()),"Rda")
