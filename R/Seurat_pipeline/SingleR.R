library(SingleR)
library(Seurat)
library(dplyr)
library(reshape2)
library(pheatmap)
library(kableExtra)

source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
lname1 = load(file = "./data/MCL_alignment20181031.Rda")

singler = CreateSinglerObject(as.matrix(MCL@data), annot = NULL, project.name=MCL@project.name,
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL)
# Did singler find all cell labels?
length(singler$singler[[1]]$SingleR.single$labels) == ncol(MCL@data)
singler$meta.data$orig.ident = MCL@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = MCL@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = MCL@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./data/singler_MCL20181031.RData")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./data/singler_MCL20181031.RData")
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
jpeg(paste0(path,"/DrawHeatmap.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single, top.n = 50))
dev.off()
jpeg(paste0(path,"/DrawHeatmap_normF.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n = 50,normalize = F))
dev.off()
#Next, we can use the fine-tuned labels to color the t-SNE plot:
       
out = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single,
                       singler$meta.data$xy,do.label=T,
                       do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
                       label.size = 4, label.repel = T,dot.size = 3,do.legend = F,alpha = 1,
                       force=2)
jpeg(paste0(path,"/PlotTsne_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
out+  ggtitle("Supervised sub-cell type labeling by HPCA")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))
dev.off()
# main types-------
out = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single.main,
                         singler$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler$singler[[1]]$SingleR.single.main$labels,
                         label.size = 5, dot.size = 2,do.legend = F,alpha = 1,
                         label.repel = T,force=2)
jpeg(paste0(path,"/PlotTsne_main1.jpeg"), units="in", width=10, height=7,
     res=600)
out+  ggtitle("Supervised main-cell type labeling by HPCA")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))
dev.off()
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
apply(ident.DF, 2, function(x) length(unique(x))) # check number of labels
#ident.DF[is.na(ident.DF)] <- "unknown"
MCL <- AddMetaData(object = MCL,
                   metadata = as.data.frame(ident.DF))
MCL <- SetAllIdent(object = MCL, id = "singler1sub")
##############################
# process singler.color
##############################
singler_colors <- readxl::read_excel("../SingleR/singler.colors.xlsx")
singler_colors = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
length(singler_colors)
length(unique(singler$singler[[1]]$SingleR.single$labels))

#' .AddMetaColor: prepare meta.data data frame to store color code
#' @mat factor at column 1, rownames is cell.names
#' @colors vector of hex colors
#' @df_colors data frame to add to meta data directly
#' @example 
# singlerDF = .AddMetaColor(mat = singler$singler[[2]]$SingleR.single$labels,
#                        colors = singler_colors[3:26])
.AddMetaColor <- function(mat, colors){
        if(class(mat) != "data.frame") mat = as.data.frame(mat)
        mat$cell.names = rownames(mat)
        mat$index <- as.numeric(as.factor(mat[,1]))
        if(length(colors)<length(unique(mat$index))) {
                stop(paste("Not enough colors! Provide at least", 
                           length(unique(mat$index)),"different colors"))}
        df_colors = data.frame(colors[1:length(unique(mat$index))],
                               "index" = 1:length(unique(mat$index)))
        mat_colors <- full_join(mat, df_colors, by = "index")
        # remove NA color if there is any
        mat_colors <- mat_colors[(mat_colors$cell.names %in% rownames(mat)),]
        rownames(mat_colors) = mat_colors$cell.names
        df_colors <- data.frame(mat_colors$colors,
                                row.names = mat_colors$cell.names)
        colnames(df_colors)[1] = paste(colnames(mat)[1],"colors",sep =".")
        
        return(df_colors)
}


#' AddMetaColor: convert one MetaData label into color scheme and store into MetaData
#' @object seurat object
#' @label colname in metadata
#' @colors vector of hex colors
# MCL <- AddMetaColor(object = MCL, label= "singler1sub", colors = singler_colors)
AddMetaColor<- function(object, label = NULL, colors = NULL){
        
        if(is.null(label)) label <- .FindIdentLabel(object)
        if(is.null(colors)) colors <- SingleR:::singler.colors
        mat = data.frame(object@meta.data[,label],
                         row.names = object@cell.names)
        colnames(mat) = get("label")
        newMetaData = .AddMetaColor(mat = mat, colors = colors)
        object <- AddMetaData(object = object, metadata = newMetaData)
        
        return(object)
}

#' FindIdentLabel: Find identical label between ident and metadata
#' @object seurat object
#' @label colname in metadata
.FindIdentLabel <- function(object){
        ident.label <- as.vector(object@ident)
        labels <- apply(object@meta.data, 2, function(x) all(ident.label %in% x))
        return(names(labels[labels]))
        }


MCL <- AddMetaColor(object = MCL, colors = singler_colors)

##############################
# draw tsne plot
##############################

df_colors <- ExtractMetaColor(object = MCL)
df_colors
p3 <- TSNEPlot.1(object = MCL, dim.1 = 1, dim.2 = 2, 
                group.by = "ident", do.return = TRUE,pt.size = 2,
                colors.use = ExtractMetaColor(MCL),no.legend = T,
                do.label =F,label.size=4, label.repel = F,force=1)+
        ggtitle("Supervised sub cell type labeling by HPCA")+
        theme(text = element_text(size=15),
              plot.title = element_text(hjust = 0.5))
jpeg(paste0(path,"PlotTsne_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
p3+xlim(-40, 35)+ylim(-35, 35)
dev.off()

save(MCL,file="./output/MCL_alignment20181031.Rda")
##############################
# subset Seurat
###############################
table(MCL@meta.data$orig.ident)
table(MCL@ident)
g <- SplitTSNEPlot(MCL,group.by = "ident",split.by = "orig.ident",
                   no.legend = T,do.label =F,label.size=3,pt.size = 2,
                   return.plots =T, label.repel = T,force=2)

jpeg(paste0(path,"SplitTSNEPlot~.jpeg"), units="in", width=10, height=7,
     res=600)
print(do.call(plot_grid, g))
dev.off()

jpeg(paste0(path,"SplitTSNEPlot_Pt-MD.jpeg"), units="in", width=10, height=7,
     res=600)
g[[1]]+xlim(-40, 35)+ylim(-35, 35)
dev.off()

jpeg(paste0(path,"SplitTSNEPlot_Pt-RM.jpeg"), units="in", width=10, height=7,
     res=600)
g[[2]]+xlim(-40, 35)+ylim(-35, 35)
dev.off()
