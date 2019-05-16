 invisible(sapply(c("Seurat","dplyr","tidyr","magrittr","dplyr","gplots"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
source("R/util.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 


# 3.1.1 load data
# Rename ident
(load(file="data/MCL_Harmony_43_20190430.Rda"))
object %<>% ScaleData
args <- commandArgs(trailingOnly = TRUE)
args[1] <- as.character(args[1])
# B cells only ================
object <- SetAllIdent(object, id="res.0.6")
table(object@ident)
B_cells_MCL <- SubsetData(object, ident.use = c(0,1,5,6,9,13,14,18,19,20))
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1main")
table(B_cells_MCL@ident)
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B_cells","MCL","HSC"))
table(B_cells_MCL@meta.data$singler1sub) %>% as.data.frame %>%
        .[.[,"Freq"] >0,]
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1sub") 
(keep <- grep("B_cells.|MCL|CLP|HSC",B_cells_MCL@meta.data$singler1sub, value = T) %>% unique)
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use =  keep)


##############################
# re-scale, PCA, harmony and tsne
##############################
if(args[1] == "NormalizeData"){
        B_cells_MCL <- NormalizeData(B_cells_MCL)
        B_cells_MCL <- FindVariableGenes(object = B_cells_MCL, mean.function = ExpMean, 
                                         dispersion.function = LogVMR, 
                                         x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.5)
        length(B_cells_MCL@var.genes)
        
        B_cells_MCL %<>% ScaleData
        B_cells_MCL %<>% RunPCA(pc.genes = B_cells_MCL@var.genes, pcs.compute = 100, do.print = F)
        
        pcs = 1:75
        system.time(B_cells_MCL %<>% RunHarmony("orig.ident", dims.use = pcs,
                                                theta = 2, plot_convergence = TRUE,
                                                nclust = 50, max.iter.cluster = 100))
}

if(args[1] == "ScaleData"){
        B_cells_MCL %<>% ScaleData
        B_cells_MCL %<>% RunPCA(pc.genes = B_cells_MCL@var.genes, pcs.compute = 100, do.print = F)
        
        pcs = 1:75
        system.time(B_cells_MCL %<>% RunHarmony("orig.ident", dims.use = pcs,
                                                theta = 2, plot_convergence = TRUE,
                                                nclust = 50, max.iter.cluster = 100))
}

if(args[1] == "RunPCA"){
        B_cells_MCL %<>% RunPCA(pc.genes = B_cells_MCL@var.genes, pcs.compute = 100, do.print = F)
        pcs = 1:75
        system.time(B_cells_MCL %<>% RunHarmony("orig.ident", dims.use = pcs,
                                                theta = 2, plot_convergence = TRUE,
                                                nclust = 50, max.iter.cluster = 100))
}

if(args[1] == "RunHarmony"){
        pcs = 1:75
        system.time(B_cells_MCL %<>% RunHarmony("orig.ident", dims.use = pcs,
                                                theta = 2, plot_convergence = TRUE,
                                                nclust = 50, max.iter.cluster = 100))
}

system.time(
        B_cells_MCL <- RunTSNE(B_cells_MCL, reduction.use = "harmony", dims.use = 1:75, 
                               perplexity = 30, do.fast = TRUE))
system.time(
        B_cells_MCL %<>% FindClusters(reduction.type = "harmony", resolution = 0.2, dims.use = 1:75,
                                      save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                                      force.recalc = TRUE, print.output = FALSE))

g1 <- TSNEPlot.1(object = B_cells_MCL, do.label = T, group.by = "ident",
                 do.return = TRUE, no.legend = T, 
                 pt.size = 1,label.size = 5 )+
        ggtitle(args[1])+
        theme(plot.title = element_text(hjust = 0.5,size = 18)) 

jpeg(paste0(path,"rerun_",args[1],".jpeg"), units="in", width=10, height=7,res=600)
print(g1)
dev.off()