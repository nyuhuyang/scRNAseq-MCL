invisible(lapply(c("Seurat","monocle","dplyr",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- "output/20200413_monocle2/"
if(!dir.exists(path)) dir.create(path, recursive = T)
#SBATCH --mem=32G
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))
run_differentialGeneTest = FALSE
# load data
object = readRDS(file = "data/MCL_41_B_20200225.rds")
object %<>% AddMetaColor(label= "X4clusters", colors = gg_color_hue(4))
object$X4_orig.ident = paste(object$orig.ident,
                             object$X4clusters, sep = "_")
sample_pairs = list(c("N01","N02","N03","N04"),
                    c("PtU01","PtU02","PtU03","PtU04"),
                    c("Pt11_LN1","Pt11_1","Pt11_14","Pt11_28"),
                    c("Pt11_LN1","Pt11_1"),
                    c("Pt11_1","Pt11_14","Pt11_28"),
                    c("Pt11_14","Pt11_28"),
                    c("Pt17_LN1","Pt17_2","Pt17_7","Pt17_12","Pt17_31"),#7
                    c("Pt17_LN1","Pt17_2","Pt17_12"),#8
                    c("Pt17_LN1","Pt17_2","Pt17_7"),#9
                    c("Pt17_LN1","Pt17_2","Pt17_7","Pt17_12"),#10
                    c("Pt17_LN1","Pt17_2"),
                    c("Pt17_2","Pt17_7"),
                    c("Pt17_2","Pt17_12"),
                    c("Pt17_7","Pt17_12"),
                    c("Pt17_2","Pt17_7","Pt17_12"),
                    c("Pt17_7","Pt17_12","Pt17_31"),
                    c("Pt25_SB1","Pt25_1","Pt25_1_8","Pt25_24","Pt25_25Pd","Pt25_AMB25Pd"),
                    c("Pt25_SB1","Pt25_1","Pt25_24","Pt25_25Pd","Pt25_AMB25Pd"),#18
                    c("Pt25_SB1","Pt25_24","Pt25_25Pd","Pt25_AMB25Pd"),#19
                    c("Pt25_24","Pt25_25Pd","Pt25_AMB25Pd"),#20
                    c("Pt25_1","Pt25_1_8","Pt25_24","Pt25_25Pd","Pt25_AMB25Pd"),#21
                    c("Pt25_SB1","Pt25_1"),
                    c("Pt25_1","Pt25_1_8"),
                    c("Pt25_1_8","Pt25_24"),
                    c("Pt25_24","Pt25_25Pd"),
                    c("Pt25_24","Pt25_AMB25Pd"),#26
                    c("Pt25_25Pd","Pt25_AMB25Pd"),
                    c("Pt25_SB1","Pt25_AMB25Pd"),
                    c("Pt25_SB1","Pt25_24"),#29
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
                    c("Pt28_LN1","Pt28_28"),
                    c("Pt25_25Pd"))
(sample = sample_pairs[[i]])
Idents(object) = "orig.ident"
object %<>% subset(idents = sample)
GC()
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
table(pData(cds)$orig.ident)
#Filtering low-quality cells
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))
length(expressed_genes)
#######################
#Trajectory step 1: choose genes that define a cell's progress
if(run_differentialGeneTest){
        clustering_DEG_genes <- differentialGeneTest(cds[expressed_genes,],
                                                     fullModelFormulaStr = "~X4_orig.ident",
                                                     cores = detectCores()/2)
        saveRDS(clustering_DEG_genes, paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_DE.rds"))

} else clustering_DEG_genes = readRDS(paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_DE.rds"))

print("clustering_DEG_genes")
cds_ordering_genes <-row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][
        1:min(nrow(clustering_DEG_genes),1000)]
cds %<>% setOrderingFilter(ordering_genes = cds_ordering_genes)
cds %<>% reduceDimension(method = 'DDRTree')
cds %<>% orderCells()

saveRDS(cds, paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_cds.rds"))
cds = readRDS(paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_cds.rds"))
#Trajectory step 2: generate Trajectory plot for all samples
group_by <- c("orig.ident","X4clusters", "Pseudotime","cell.types")
if(length(unique(cds$orig.ident)) == 5) orig.ident_color <- c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")
if(length(unique(cds$orig.ident)) == 4) orig.ident_color <- c("#F8766D","#A3A500","#00BF7D","#00B0F6")
if(length(unique(cds$orig.ident)) == 3) orig.ident_color <- c("#F8766D","#A3A500","#00BF7D")

for(k in seq_along(group_by)){
        jpeg(paste0(save.path,"trajectory_",paste(sample, collapse = "-"),"_",group_by[k],".jpeg"),
             units="in", width=7, height=7,res=600)
        g <- plot_cell_trajectory(cds, color_by = group_by[k],cell_size = 1,
                                  show_branch_points = FALSE)+
                guides(colour = guide_legend(override.aes = list(size=5))) #larger legend label

        if(group_by[k] == "orig.ident") g = g + scale_color_manual(values=orig.ident_color)
        if(group_by[k] == "X4clusters") g = g + scale_color_manual(values=X4clusters_color)
        if(group_by[k] == "cell.types") g = g +
                scale_color_manual(values=sort(unique(cds$cell.types.colors),decreasing = T))

        print(g)+TitleCenter()
        dev.off()
        Progress(k, length(group_by))
}


#Trajectory step 2: generate Trajectory plot by each sample
save.path.sub = paste0(save.path, "subset/")
if(!dir.exists(save.path.sub)) dir.create(save.path.sub, recursive = T)

samples <- unique(cds$orig.ident)
group_by <- c("X4clusters", "Pseudotime","cell.types")

for(s in seq_along(samples)){
        valid_cells <- row.names(subset(pData(cds), orig.ident == samples[s]))
        subset_cds <- cds[,valid_cells]
        subset_cds@reducedDimS <- subset_cds@reducedDimS[,valid_cells]
        for(k in seq_along(group_by)){
                jpeg(paste0(save.path.sub,"trajectory_",samples[s],"_",group_by[k],".jpeg"),
                     units="in", width=7, height=7,res=600)
                g <- plot_cell_trajectory(subset_cds, cell_size = 1,
                                          color_by = group_by[k],show_branch_points = FALSE)+
                        guides(colour = guide_legend(override.aes = list(size=5))) #larger legend label

                if(group_by[k] == "X4clusters") g = g + scale_color_manual(values=X4clusters_color)
                if(group_by[k] == "cell.types") g = g +
                        scale_color_manual(values=sort(unique(cds$cell.types.colors),decreasing = T))
                g = g + ggtitle(paste(group_by[k], "in", samples[s])) + TitleCenter()
                print(g)
                dev.off()
                Progress((s-1)*length(samples)+k, length(samples)*length(group_by))
        }
}

#Trajectory step 3: generate Trajectory plot by each cluster
X4clusters <- sort(unique(cds$X4clusters))
file_name = paste(c(sample[1],gsub(".*_","",sample[2:length(sample)])),collapse = "_")
sample_name = paste(c(sample[1],gsub(".*_","",sample[2:length(sample)])),collapse = ", ")

for(c in seq_along(X4clusters)){
        valid_cells <- colnames(cds)[cds$X4clusters %in% X4clusters[c]]
        subset_cds <- cds[,valid_cells]
        subset_cds@reducedDimS <- subset_cds@reducedDimS[,valid_cells]
        jpeg(paste0(save.path.sub,"trajectory_",X4clusters[c],"_",file_name,".jpeg"),
             units="in", width=7, height=7,res=600)
        g <- plot_cell_trajectory(subset_cds, cell_size = 1,
                                  color_by = "orig.ident",show_branch_points = FALSE)
        g = g + scale_color_manual(values=orig.ident_color)
        g = g + ggtitle(paste(X4clusters[c], "in", sample_name)) + TitleCenter()
        print(g)
        dev.off()
        Progress(c, length(group_by))
}

my_gene = c("CCND1","CDK4","MCM7","EZH2","EZH1")

# plot genes that are expressed in a branch dependent manner ===============
jpeg(paste0(path,file_name,"_branched_pseudotime.jpeg"),
     units="in", width=7, height=7,res=600)
plot_cell_trajectory(cds["CCND1",], color_by = "CellType") +
        theme(legend.position = "right")
dev.off()
