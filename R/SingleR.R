library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
lname1 = load(file = "./data/MCL_20181105.Rda");lname1
singler = CreateSinglerObject(MCL@data, annot = NULL, project.name=MCL@project.name,
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL)
# singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(MCL@data)
all.cell = MCL@cell.names;length(all.cell)
know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
MCL = SubsetData(MCL, cells.use = know.cell)
MCL
singler$meta.data$orig.ident = MCL@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = MCL@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = MCL@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_MCL_20181105.RData")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./output/singler_MCL_20181020.RData")
lnames = load(file = "../SingleR/data/Hpca.RData")
lnames
singler$seurat = MCL
SingleR.DrawScatter(sc_data = singler$seurat@data,cell_id = 10, 
                    ref = Hpca, sample_id = 232)

# Step 2: Multiple correlation coefficients per cell types are aggregated 
# to provide a single value per cell type per single-cell. 
# In the examples below we use the 80% percentile of correlation values.
# for visualization purposes we only present a subset of cell types (defined in labels.use)
out = SingleR.DrawBoxPlot(sc_data = singler$seurat@data,cell_id = 10, 
                          ref = immgen,main_types = T,
                          labels.use=c('B cells','T cells','DC','Macrophages','Monocytes','NK cells',
                                       'Mast cells','Neutrophils','Fibroblasts','Endothelial cells'))
print(out$plot)

##############################
# Human Primary Cell Atlas (HPCA)
###############################
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, top.n = Inf)
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"/DrawHeatmap.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50))
dev.off()
jpeg(paste0(path,"/DrawHeatmap_normF.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,normalize = F))
dev.off()
#Next, we can use the fine-tuned labels to color the t-SNE plot:
       
out1 = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single,
                         singler$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
                         label.size = 5, dot.size = 2,do.legend = F,alpha = 1,
                         label.repel = T,force=2)
#label.repel = T,force=2)
jpeg(paste0(path,"/PlotTsne_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
out1+  ggtitle("Supervised sub cell type labeling by HPCA")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle
dev.off()
# main types-------
out2 = SingleR.PlotTsne.1(singler$singler[[2]]$SingleR.single,
                         singler$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
                         label.size = 5, dot.size = 2,do.legend = F,alpha = 1,
                         label.repel = T,force=2)
jpeg(paste0(path,"/PlotTsne_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
out2+  ggtitle("Supervised cell type labeling by Blueprint+ Encode")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$singler[[2]]$SingleR.single$labels,
            singler$meta.data$orig.ident)) %>%
        kable_styling()
kable(table(singler$meta.data$orig.ident,
            singler$singler[[1]]$SingleR.single.main$labels)) %>%
        kable_styling()

##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub"=singler$singler[[1]]$SingleR.single$labels,
                       "singler1main"=singler$singler[[1]]$SingleR.single.main$labels,
                       "singler2sub"=singler$singler[[2]]$SingleR.single$labels,
                       "singler2main"=singler$singler[[2]]$SingleR.single.main$labels,
                       "cell.names" = rownames(singler$singler[[1]]$SingleR.single$labels))
knowDF = data.frame("cell.names"= MCL@cell.names)
ident.DF = full_join(singlerDF,knowDF, by="cell.names")
ident.DF<- apply(ident.DF,2,as.character)
rownames(ident.DF) = ident.DF[,"cell.names"]
ident.DF = ident.DF[,-which(colnames(ident.DF) == "cell.names")]
#ident.DF[is.na(ident.DF)] <- "unknown"
MCL <- AddMetaData(object = MCL,
                   metadata = as.data.frame(ident.DF))
MCL <- SetAllIdent(object = MCL, id = "singler2sub")
##############################
# process color scheme
##############################
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
length(singler_colors1);length(singler_colors2)

apply(MCL@meta.data[,c("singler1sub","singler1main","singler2sub","singler2main")],
      2,function(x) length(unique(x)))

MCL <- AddMetaColor(object = MCL, label= "singler2sub", colors = singler_colors1[3:26])
MCL <- AddMetaColor(object = MCL, label= "singler1sub", colors = singler_colors2[1:42])
MCL <- SetAllIdent(object = MCL, id = "singler1sub")
TSNEPlot.1(MCL, colors.use = ExtractMetaColor(MCL))

# 2018-11-05
singler2sub <- as.data.frame(table(MCL@meta.data$singler2sub,MCL@meta.data$orig.ident))
singler2sub <- singler2sub %>% spread(Var2, Freq)
rownames(singler2sub) <- singler2sub$Var1
singler2sub <- singler2sub[,-1]
rowSums(singler2sub) %>% kable() %>% kable_styling()
write.csv(MCL@meta.data,file = "output/MCL_metadata_20181106.csv")

##############################
# draw tsne plot
##############################

df_colors <- ExtractMetaColor(object = MCL, label = "singler2sub")
#colors_df <- gg_colors(object = MCL, colors.use = singler_colors[3:26])
p3 <- DimPlot.1(object = MCL, reduction.use = "tsne", dim.1 = 1, dim.2 = 2, 
          group.by = "ident", do.return = TRUE,
          cols.use = ExtractMetaColor(object = MCL),no.legend = T,
          do.label =T,label.size=4, label.repel = T,force=1)+
        ggtitle("Supervised sub cell type labeling by HPCA")+
        theme(text = element_text(size=15),
              plot.title = element_text(hjust = 0.5))

jpeg(paste0(path,"PlotTsne_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

save(MCL,file="./output/MCL_20181106.RData")
##############################
# subset Seurat
###############################
table(MCL@meta.data$orig.ident)
table(MCL@ident)

df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
tests <- paste0("test",1:4)
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$samples[sample_n]
        print(samples)
        
        cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(MCL, cells.use = cell.use)
        
        g <- SplitTSNEPlot(subset.MCL,group.by = "ident",split.by = "orig.ident",
                           no.legend = T,do.label =F,label.size=3,
                           return.plots =T, label.repel = T,force=2)
        jpeg(paste0(path,test,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, g))
        dev.off()
}


for(test in tests){
  sample_n = which(df_samples$tests %in% test)
  samples <- df_samples$samples[sample_n]
  print(samples)
  
  cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
  subset.MCL <- SubsetData(MCL, cells.use = cell.use)
  
  subset.MCL <- RunTSNE(object = subset.MCL, reduction.use = "MNN", dims.use = 1:50,
                        do.fast = TRUE, perplexity= 30)
  
  g <- SplitTSNEPlot(subset.MCL,group.by = "ident",split.by = "orig.ident",
                     no.legend = T,do.label =F,label.size=3,
                     return.plots =T, label.repel = T,force=2)
  jpeg(paste0(path,test,"_reshape.jpeg"), units="in", width=10, height=7,
       res=600)
  print(do.call(plot_grid, g))
  dev.off()
}
