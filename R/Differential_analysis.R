########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file="data/MCL_Harmony_20181118.Rda"))
MCL <- SetAllIdent(MCL, id="singler1main")
B_cells_MCL <- SubsetData(MCL, ident.use = c("B_cell","Pre-B_cell_CD34-",
                                             "Pro-B_cell_CD34+"))
table(B_cells_MCL@meta.data$singler2main)
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler2main")
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B-cells","HSC"))

# Compare FACS data using GenePlot =======
df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
sample_n = which(df_samples$tests %in% paste0("test",2:4))
table(df_samples$tests)
df_samples[sample_n,] %>% kable() %>% kable_styling()

tests <- paste0("test",2:4)
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$samples[sample_n]
        print(samples)
        
        cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(MCL, cells.use = cell.use)
        
        cells_use_list <- split(row.names(subset.MCL@meta.data), subset.MCL@meta.data[,"orig.ident"])
        G1 <- lapply(cells_use_list, function(cells_use) {
                single.MCL <- SubsetData(subset.MCL, cells_use_list[[2]])
                GenePlot.1(single.MCL, gene1 = "CD5", gene2 = "CD19",use.raw = F, 
                                 title = unique(single.MCL@meta.data$orig.ident))
        })
        G2 <- lapply(cells_use_list, function(cells_use) {
                single.MCL <- SubsetData(subset.MCL, cells_use)
                GenePlot.1(single.MCL, gene1 = "CD3E", gene2 = "CD19",use.raw = F, 
                                 title = unique(single.MCL@meta.data$orig.ident))
        })
        
        jpeg(paste0(path,test,"_CD5_CD3E_CD19.jpeg"), units="in", width=8, height=7,
             res=600)
        #print(plot_grid(G1[[1]],G1[[5]],G1[[4]],G1[[3]],G1[[2]],
        #          G2[[1]],G2[[5]],G2[[4]],G2[[3]],G2[[2]],nrow = 2))
        print(plot_grid(G1[[1]],G1[[2]],G2[[1]],G2[[2]],nrow = 2))
        dev.off()
}
#  count cells in geneplot=============================
CountCells <- function(object = subset.MCL, gene1="CD5", gene2 = "CD19", split.by ="orig.ident"){
        cells_use_list <- split(row.names(subset.MCL@meta.data), subset.MCL@meta.data[,split.by])
        cell_counts <- lapply(cells_use_list, function(cells_use) {
                single.MCL <- SubsetData(subset.MCL, cells_use)
                c(length(which(single.MCL@data[gene1,] == 0 & single.MCL@data[gene2,] == 0)),
                  length(which(single.MCL@data[gene1,] > 0 & single.MCL@data[gene2,] == 0)),
                  length(which(single.MCL@data[gene1,] == 0 & single.MCL@data[gene2,] > 0)),
                  length(which(single.MCL@data[gene1,] > 0 & single.MCL@data[gene2,] > 0)))/
                        length(single.MCL@cell.names)*100
        })
        df_cell_counts <- list2df(cell_counts)
        rownames(df_cell_counts) <- c("l.l","h.l","l.h","h.h")
        return(df_cell_counts)
}

tests <- paste0("test",2:4)
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$samples[sample_n]
        print(samples)
        
        cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(MCL, cells.use = cell.use)
        CD5_CD19 <- CountCells(subset.MCL, "CD5","CD19")
        CD3E_CD19 <- CountCells(subset.MCL, "CD3E","CD19")
        df <- merge(CD5_CD19,CD3E_CD19, by="row.names",all.x=TRUE)
        write.csv(df, paste0(path,test,"_countcell.csv"),row.names = F)
}

#  demonstrate cell cycle in geneplot=============================
tests <- paste0("test",3:4)
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$samples[sample_n]
        print(samples)
        
        cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(MCL, cells.use = cell.use)
        
        cells_use_list <- split(row.names(subset.MCL@meta.data), subset.MCL@meta.data[,"orig.ident"])
        G1 <- lapply(cells_use_list, function(cells_use) {
                single.MCL <- SubsetData(subset.MCL, cells_use)
                GenePlot.1(single.MCL, gene1 = "G2M.Score", gene2 = "S.Score",use.raw = F, 
                           title = unique(single.MCL@meta.data$orig.ident))+
                        xlim(0, 1.25)+ylim(0, 0.8)
        })
        jpeg(paste0(path,test,"_cellcyle_geneplot.jpeg"), units="in", width=8, height=7,
             res=600)
        #print(plot_grid(G1[[1]],G1[[5]],G1[[4]],G1[[3]],G1[[2]]))
        print(do.call(plot_grid,G1))
        dev.off()
}



# geom_density  ===========
markers <-  HumanGenes(MCL,c("CD19","MS4A1","SOX11","ITGA4","CCND1"))
MCL <- B_cells_MCL
tests <- paste0("test",3:4)
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$samples[sample_n]
        print(samples)
        
        cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(MCL, cells.use = cell.use)
        
        g <- split(rownames(subset.MCL@meta.data), subset.MCL@meta.data[,"orig.ident"]) %>% lapply(function(cells_use) {
                single.MCL <- SubsetData(subset.MCL, cells.use = cells_use)
                sample <- unique(single.MCL@meta.data$orig.ident)
                data.use <- single.MCL@data[markers,] %>% as.matrix %>% t %>% as.data.frame %>%
                       gather(key = markers, value = ave.expr)
                ggplot(data.use, aes(x = ave.expr, fill = markers)) +
                        geom_density(alpha = .5) + scale_y_sqrt() + ylim(0, 3)+
                        xlab("Average expression (log nUMI)")+
                        ggtitle(sample)+
                        theme(text = element_text(size=15),
                              legend.position=c(0.4,0.8) ,# legend.position="none", 
                              plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
        })
        jpeg(paste0(path,"density_",test,".jpeg"), units="in", width=10, height=7,res=600)
        print(do.call(plot_grid,g)+ #  plot_grid(g[[1]],g[[5]],g[[4]],g[[3]],g[[2]])
                ggtitle("Density plot for typical markers in MCL B cells")+
                theme(text = element_text(size=15),							
                      plot.title = element_text(hjust = 0.5,size = 15, face = "bold")))
        dev.off()
}

#======== compare threshold ===============
sample = "Pt-1475"
marker = "MS4A1"
thresholds = c(0.01,0.1,0.5,1,2,3,4,5)
cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident == sample]
subset.MCL <- SubsetData(MCL, cells.use = cell.use)

for(threshold in thresholds){
        g <- list()
        g[[1]] <- SingleFeaturePlot.1(object = MD.MCL, threshold=threshold,
                                      feature = marker,title = "MD")
        g[[2]] <- SingleFeaturePlot.1(object = subset.MCL, threshold=threshold,
                                      feature = marker,title = sample)
        jpeg(paste0(path,"Splited_MD_",sample,"_",marker,"_",threshold,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, g)+
                      ggtitle(paste("threshold =",threshold))+
                      theme(text = element_text(size=20),							
                            plot.title = element_text(hjust = 0.5,size = 15, face = "bold")))
        print(paste0(which(threshold == thresholds),":",length(thresholds)))
        dev.off()
}



