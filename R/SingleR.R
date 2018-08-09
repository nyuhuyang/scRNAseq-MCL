library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)

source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
#====== 3.1 Create Singler Object  ==========================================
lname1 = load(file = "./data/MCL_alignment.Rda")
lname1
lname2 = load(file = "../SingleR/data/Hpca.RData")
lname2
singler = CreateSinglerObject(as.matrix(MCL@data), annot = NULL, project.name=MCL@project.name,
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              ref.list = list(Hpca), normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL)
singler$meta.data$orig.ident = MCL@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = MCL@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = MCL@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./data/singler_MCL1.RData")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./data/singler_MCL.RData")
lnames
SingleR.DrawScatter(sc_data = singler$seurat@data,cell_id = 10, 
                    ref = immgen, sample_id = 232)

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
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, top.n = Inf,
                    clusters = singler$meta.data$orig.ident)
#Or by all cell types (showing the top 50 cell types):
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,
                    clusters = singler$meta.data$orig.ident)
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,
                    normalize = F,clusters = singler$meta.data$orig.ident)
#Next, we can use the fine-tuned labels to color the t-SNE plot:
       
out = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single,
                       singler$meta.data$xy,do.label=T,
                       do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
                       label.size = 4, label.repel = T,dot.size = 2,do.legend = F,alpha = 1,
                       force=2)
out$p+  ggtitle("Supervised sub-cell type labeling by HPCA")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle
# main types-------
out = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single.main,
                         singler$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler$singler[[1]]$SingleR.single.main$labels,
                         label.size = 5, dot.size = 2,do.legend = F,alpha = 1,
                         label.repel = T,force=2)
out$p+  ggtitle("Supervised cell type labeling by HPCA")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle

singler$seurat = MCL
SplitSingleR.PlotTsne(singler = singler, split.by = "conditions",do.label=T,alpha = 1,
                      select.plots =c(2,1),label.repel = T, do.legend = FALSE,
                      show.subtype = TRUE,force=20)
output <- SplitSingleR.PlotTsne(singler = singler, split.by = "conditions",main=T,
                              return.plots=T,do.label=T,do.legend = F,alpha = 1,
                              label.repel = T, force=2)
output[[1]]
output[[2]]
#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$singler[[1]]$SingleR.single$labels,
            singler$meta.data$orig.ident)) %>%
        kable_styling()
kable(table(singler$meta.data$orig.ident,
            singler$singler[[1]]$SingleR.single.main$labels)) %>%
        kable_styling()
kable(table(singler$meta.data$orig.ident,singler$seurat@ident)) %>%
        kable_styling()

##############################
# Blueprint+Encode
###############################
SplitSingleR.PlotTsne(singler = singler, split.by = "conditions",do.label=T,alpha = 1,
                      select.plots =c(2,1),label.repel = F, do.legend = FALSE,
                      show.subtype = F,force=3)
output <- SplitSingleR.PlotTsne(singler = singler, split.by = "conditions",show.subtype = F,
                                return.plots=T,do.label=F,do.legend = T,legend.size = 15,
                                alpha = 1,label.repel = T, force=1)
        
output[[1]]
output[[2]]

kable(table(singler$singler[[2]]$SingleR.single$labels,
            singler$meta.data$orig.ident)) %>%
        kable_styling()

kable(table(singler$other,singler$meta.data$orig.ident)) %>%
        kable_styling()
SingleFeaturePlot.1(object = MCL, "FOXP3", threshold = 0.1)

