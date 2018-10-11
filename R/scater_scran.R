########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(DropletUtils)
library(scater)
library(EnsDb.Hsapiens.v86)
library(scran)
library(Matrix)
library(pheatmap)
########################################################################
#
#  0. scater
# 
# ######################################################################
# 0.1. Setting up the data
# 0.1.1 Reading in a sparse matrix
df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
sample_n = which(df_samples$Patients %in% "test1")
df_samples[sample_n,]
samples <- df_samples$Samples[sample_n]
projects <- df_samples$Projects[sample_n]
conditions <- df_samples$Conditions[sample_n]

sce_list <- list()
for(i in 1:length(samples)){
        fname <- paste0("./data/",samples[i],
                        "/outs/raw_gene_bc_matrices/hg19/")
        sce_list[[i]] <- read10xCounts(fname, col.names=TRUE)
}

# 0.1.2 Annotating the rows
for(i in 1:length(samples)){
        rownames(sce_list[[i]]) <- uniquifyFeatureNames(rowData(sce_list[[i]])$ID,
                                                        rowData(sce_list[[i]])$Symbol)
        print(head(rownames(sce_list[[i]]),3))
        print(length(rownames(sce_list[[i]])))
}

# We also identify the chromosomal location for each gene. 
# The mitochondrial percentage is particularly useful for later quality control.
for(i in 1:length(samples)){
        location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce_list[[i]])$ID, 
                           column="SEQNAME", keytype="GENEID")
        rowData(sce_list[[i]])$CHR <- location
        print(summary(location=="MT"))
}

# 0.2 Calling cells from empty droplets
## ----rankplot, "Total UMI count for each barcode in the dataset, 
# plotted against its rank (in decreasing order of total counts). 
# The inferred locations of the inflection and knee points are also shown."----
bcrank <- lapply(sce_list, function(x) barcodeRanks(counts(x)))

# Only showing unique points for plotting speed.
uniq <- lapply(bcrank, function(x) !duplicated(x$rank))

# We call cells at a false discovery rate (FDR) of 1%, 
# meaning that no more than 1% of our called barcodes should be empty droplets on average.
set.seed(100)
e.out <- lapply(sce_list, function(x) emptyDrops(counts(x))) # long time

########################################################################
## --------------------------------------------------------------------------
# generate plots for QC
g <- list()
for(i in 1:length(samples)){
        g[[i]] <- qplot(bcrank[[i]]$rank[uniq[[i]]], bcrank[[i]]$total[uniq[[i]]], log="xy",
             xlab="Rank", ylab="Total UMI count", main = samples[i])+
                theme(text = element_text(size = 20),
                      plot.title = element_text(hjust = 0.5))+
                geom_vline(aes(xintercept = length(which(e.out[[i]]$FDR <= 0.01)),
                               linetype = "FDR <= 0.1"),
                           colour= 'black')+
                geom_hline(aes(yintercept= c(bcrank[[i]]$inflection,#lower
                                             bcrank[[i]]$knee), #higher
                               linetype = c("inflection","knee")),
                               colour = c("red","blue"))+
        scale_linetype_manual(name = "Threshold", values = c(3,2,1))
}
# alternative
par(mfrow =c(3,2))
for(i in 1:length(samples)){
        plot(bcrank[[i]]$rank[uniq[[i]]], bcrank[[i]]$total[uniq[[i]]], log="xy",
             xlab="Rank", ylab="Total UMI count", main=samples[i], cex.lab=1.2)
        
        abline(h=bcrank[[i]]$knee, col="green", lty=1)
        abline(h=bcrank[[i]]$inflection, col="blue", lty=1)
        abline(v=length(which(e.out[[i]]$FDR <= 0.01)), col="red", lty=1)
        
        legend("bottomleft", legend=c("Knee","Inflection","FDR <= 0.01"), 
               col=c("green","blue","red" ), lty=1, cex=1.2)
}
## --------------------------------------------------------------------------
########################################################################
# 0.3
path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
dir.create(path, recursive = T)
# emptyDrops() computes Monte Carlo p-values, 
# so it is important to set the random seed to obtain reproducible results. 
# The number of Monte Carlo iterations also determines the lower bound for the _p_values. 
# If any non-significant barcodes are TRUE for Limited, 
# we may need to increase the number of iterations to ensure that they can be detected.
for(i in 1:length(samples)){
        print(samples[i])
        print(table(Sig=e.out[[i]]$FDR <= 0.01, Limited=e.out[[i]]$Limited))
}
# using which() to automatically remove NAs.
for(i in 1:length(samples)){
        sce_list[[i]] <- sce_list[[i]][,which(e.out[[i]]$FDR <= 0.01)]
}

# 0.4 Quality control on the cells
# It is entirely possible for droplets to contain damaged or dying cells,
# which need to be removed prior to downstream analysis. 
# We compute some QC metrics using  calculateQCMetrics() (McCarthy et al. 2017) 
# and examine their distributions in Figure 2.
sce_list.copy <- sce_list
sce_list <- lapply(sce_list.copy, function(x) calculateQCMetrics(x,compact = TRUE,
                        feature_controls=list(Mito=which(location=="MT"))))
QC = lapply(sce_list,function(x) x$scater_qc)
########################################################################
## --------------------------------------------------------------------------
## ----qchist, Histograms of QC metric distributions in the dataset."----
for(i in 1:length(samples)){
        jpeg(paste0(path,"/0_QC_",samples[i],".jpeg"), units="in", width=10, height=7,
             res=600)
        par(mfrow=c(1,3))
                hist(sce_list[[i]]$log10_total_counts, breaks=20, col="grey80",
             xlab="Log-total UMI count",
             main = paste("nUMI count of",samples[i]))
        hist(sce_list[[i]]$log10_total_features_by_counts, breaks=20, col="grey80",
             xlab="Log-total number of expressed features",
             main = paste("nGene count of",samples[i]))
        hist(sce_list[[i]]$pct_counts_Mito, breaks=20, col="grey80",
             xlab="Proportion of reads in mitochondrial genes",
             main = paste("mitochondrial % of",samples[i]))
        dev.off()
}
## --------------------------------------------------------------------------
########################################################################

# Ideally, we would remove cells with low library sizes or total number of expressed features as described previously.
# However, this would likely remove cell types with low RNA content,
# especially in a heterogeneous population with many different cell types.
# Thus, we use a more relaxed strategy and only remove cells with large mitochondrial proportions,
# using it as a proxy for cell damage. 
# (Keep in mind that droplet-based datasets usually do not have spike-in RNA.)
# Low-quality cells are defined as those with extreme values for these QC metrics and are removed.
for(i in 1:length(samples)){
        high.mito <- isOutlier(QC[[i]]$feature_control_Mito$pct_counts, nmads=3, type="higher")
        low.lib <- isOutlier(QC[[i]]$all$log10_total_counts, type="lower", nmad=3)
        low.genes <- isOutlier(QC[[i]]$all$log10_total_features_by_counts, type="lower", nmad=3)
        discard <- high.mito | low.lib | low.genes
        data.frame(HighMito= sum(high.mito),LowLib=sum(low.lib), 
                   LowNgenes=sum(low.genes),Discard=sum(discard))
        sce_list[[i]] <- sce_list[[i]][,!discard]
        print(summary(discard))
}
# 0.5 Examining gene expression
## Histogram of the log~10~-average counts for each gene in the dataset.----
#par(mfrow = c(3,2))
for(i in 1:length(samples)){
        ave <- calcAverage(sce_list[[i]])
        rowData(sce_list[[i]])$AveCount <- ave
        #hist(log10(ave), col="grey80",main = paste("log10(ave) of",samples[i]))
}
########################################################################
## --------------------------------------------------------------------------
## Percentage of total counts assigned to the top 50 most highly-abundant features in the dataset. 
# For each feature, each bar represents the percentage assigned to that feature for a single cell,
# while the circle represents the average across all cells. 
# Bars are coloured by the total number of expressed features in each cell."----
gg <- list()
for(i in 1:length(samples)){
        gg[[i]] <- plotHighestExprs(sce_list[[i]])
        jpeg(paste0(path,"/0_plotHighestExprs_",samples[i],".jpeg"), units="in", width=10, height=7,
             res=600)
        print(gg[[i]])
        dev.off()
}
## --------------------------------------------------------------------------
########################################################################
# Use natural Log transform to fit Seurat

for(i in 1:length(samples)){
        logcounts(sce_list[[i]]) <- as(log(assay(sce_list[[i]], "counts")+1),"dgCMatrix")
}
# 0.6 Normalizing for cell-specific biases
clusters <- list()
for(i in 1:length(samples)){
        clusters[[i]] <- quickCluster(sce_list[[i]], method="igraph", min.mean=0.1,
                                      assay.type = "logcounts",
                                 irlba.args=list(maxit=1000)) # for convergence.
        print(table(clusters[[i]]))
        sce_list[[i]] <- computeSumFactors(sce_list[[i]], min.mean=0.1, 
                                           cluster=clusters[[i]])
        print(summary(sizeFactors(sce_list[[i]])))
        
}

## ----sfplot, fig.cap="Size factors for all cells in the PBMC dataset, plotted against the library size."----
jpeg(paste0(path,"/0_Size_factors.jpeg"), units="in", width=10, height=7,
     res=600)
par(mfrow=c(3,2))
for(i in 1:length(samples)) plot(sce_list[[i]]$scater_qc$all$total_counts,
                                 sizeFactors(sce_list[[i]]), log="xy")
dev.off()

sce_list <- lapply(sce_list, function(x) normalize(x,exprs_values = "logcounts"))

# 7 Modelling the mean-variance trend
# The lack of spike-in transcripts complicates the modelling of the technical noise.
# One option is to assume that most genes do not exhibit strong biological variation,
# and to fit a trend to the variances of endogenous genes. 
# However, this assumption is generally unreasonable for a heterogeneous population. 
# Instead, we assume that the technical noise is Poisson and create a fitted trend on that basis 
# using the makeTechTrend() function.
#for(i in 1:length(samples))  metadata(sce_list[[i]])$log.exprs.offset =1

new.trend <- lapply(sce_list, function(y) makeTechTrend(x=y))

## ----trendplot, "Variance of normalized log-expression values for each gene in the dataset,
# plotted against the mean log-expression. 
# The blue line represents the mean-dependent trend fitted to the variances, 
# while the red line represents the Poisson noise."----
fit <- lapply(sce_list, function(x) trendVar(x, use.spikes=FALSE, loess.args=list(span=0.05)))
par(mfrow=c(1,1))
plot(fit[[1]]$mean, fit[[1]]$var, pch=16)
curve(fit[[1]]$trend(x), col="dodgerblue", add=TRUE)
curve(new.trend[[1]](x), col="red", add=TRUE)

# decompose the variance for each gene using the Poisson-based trend, 
# and examine the genes with the highest biological components.
## --------------------------------------------------------------------------
fit0 <- fit[[1]]
fit[[1]]$trend <- new.trend[[1]]
dec <- decomposeVar(fit=fit[[1]])
top.dec <- dec[order(dec$bio, decreasing=TRUE),] 
head(top.dec)

## ----hvgplot, "Distributions of normalized log-expression values for the top 10 genes 
# with the largest biological components in the dataset. 
# Each point represents the log-expression value in a single cell."----
plotExpression(sce_list[[1]], features=rownames(top.dec)[1:10])

## --------------------------------------------------------------------------
sce <- denoisePCA(sce, technical=new.trend, approx=TRUE)
ncol(reducedDim(sce, "PCA"))

## ----screeplot, fig.cap="Variance explained by each principal component in the PBMC dataset. The red line represents the chosen number of PCs."----
plot(attr(reducedDim(sce), "percentVar"), xlab="PC",
     ylab="Proportion of variance explained")
abline(v=ncol(reducedDim(sce, "PCA")), lty=2, col="red")

## ----pcaplot-init, fig.cap="Pairwise PCA plots of the first three PCs in the PBMC dataset, constructed from normalized log-expression values of genes with positive biological components. Each point represents a cell, coloured by the log-number of expressed features.", fig.width=9----
plotPCA(sce, ncomponents=3, colour_by="log10_total_features_by_counts")

## ----tsneplot-init, fig.cap="_t_-SNE plots constructed from the denoised PCs of the PBMC dataset. Each point represents a cell and is coloured according to the log-number of expressed features."----
sce <- runTSNE(sce, use_dimred="PCA", perplexity=30, rand_seed=100)
plotTSNE(sce, colour_by="log10_total_features_by_counts")

## --------------------------------------------------------------------------
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)
sce$Cluster <- factor(clusters$membership)
table(sce$Cluster)

## ----clustermod, fig.cap="Heatmap of the log~10~-ratio of the total weight between nodes in the same cluster or in different clusters, relative to the total weight expected under a null model of random links."----
cluster.mod <- clusterModularity(snn.gr, sce$Cluster, get.values=TRUE)
log.ratio <- log2(cluster.mod$observed/cluster.mod$expected + 1)


pheatmap(log.ratio, cluster_rows=FALSE, cluster_cols=FALSE, 
         color=colorRampPalette(c("white", "blue"))(100))

## ----tsneplot-cluster, fig.cap="_t_-SNE plots constructed from the denoised PCs of the PBMC dataset. Each point represents a cell and is coloured according to its cluster identity."----
plotTSNE(sce, colour_by="Cluster")

## --------------------------------------------------------------------------
markers <- findMarkers(sce, clusters=sce$Cluster, direction="up")

## --------------------------------------------------------------------------
marker.set <- markers[["1"]]
head(marker.set[,1:8], 10) # only first 8 columns, for brevity

## ----heatmap, fig.wide=TRUE, fig.cap="Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 1 in the PBMC dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend."----
chosen <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=chosen, exprs_values="logcounts", 
            zlim=5, center=TRUE, symmetric=TRUE, cluster_cols=FALSE,
            colour_columns_by="Cluster", columns=order(sce$Cluster))

## --------------------------------------------------------------------------
saveRDS(sce, file="pbmc_data.rds")

## --------------------------------------------------------------------------
sessionInfo()

