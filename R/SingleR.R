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
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
lname1 = load(file = "./data/MCL_CCA_20181107.Rda");lname1
singler = CreateSinglerObject(MCL@data, annot = NULL, project.name=MCL@project.name,
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              normalize.gene.length = F, variable.genes = "de",
                              fine.tune = T, do.signatures = F, clusters = NULL)
# singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(MCL@data)
all.cell = MCL@cell.names;length(all.cell)
know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
MCL = SubsetData(MCL, cells.use = know.cell)
MCL
singler$meta.data$orig.ident = MCL@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = MCL@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = MCL@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_MCL_20181107.RData")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./output/singler_MCL_20181107.RData")
lnames = load(file = "../SingleR/data/Hpca.RData")
lnames
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
apply(ident.DF,2,function(x) length(unique(x)))
#ident.DF[is.na(ident.DF)] <- "unknown"
MCL <- AddMetaData(object = MCL,
                   metadata = as.data.frame(ident.DF))
MCL <- SetAllIdent(object = MCL, id = "singler1sub")

##############################
# Human Primary Cell Atlas (HPCA)
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"/DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"/DrawHeatmap_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n = 50,normalize = F))
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
                         do.letters = F,labels = singler$singler[[2]]$SingleR.single$labels,
                         label.size = 5, dot.size = 2,do.legend = F,alpha = 1,
                         label.repel = T,force=2)
jpeg(paste0(path,"/PlotTsne_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
out2+  ggtitle("Supervised cell type labeling by Blueprint+ Encode")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$singler[[1]]$SingleR.single$labels,
            singler$meta.data$orig.ident)) %>%
        kable_styling()
kable(table(singler$singler[[2]]$SingleR.single$labels,
            singler$meta.data$orig.ident)) %>%
  kable_styling()

MCL@meta.data$singler1sub %>% table() %>% kable() %>%
  kable_styling()

MCL@meta.data$singler2sub %>% table() %>% kable() %>%
  kable_styling()

##############################
# process color scheme
##############################
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
singler_colors1[duplicated(singler_colors1)];singler_colors2[duplicated(singler_colors2)]
length(singler_colors1);length(singler_colors2)
apply(MCL@meta.data[,c("singler1sub","singler1main","singler2sub","singler2main")],
      2,function(x) length(unique(x)))
MCL@meta.data[,c("singler2sub")] %>% table() %>% kable() %>% kable_styling()
MCL <- AddMetaColor(object = MCL, label= "singler1sub", colors = singler_colors1[1:37])
MCL <- AddMetaColor(object = MCL, label= "singler2sub", colors = singler_colors2[1:21])
MCL <- SetAllIdent(object = MCL, id = "singler1sub")
TSNEPlot.1(MCL, colors.use = ExtractMetaColor(MCL),no.legend = F)

write.csv(MCL@meta.data,file = "output/MCL_metadata_20181107.csv")

##############################
# draw tsne plot
##############################
p3 <- TSNEPlot.1(object = MCL, do.label = F, group.by = "ident", 
                 do.return = TRUE, no.legend = F, 
                 colors.use = ExtractMetaColor(MCL),
                 pt.size = 1,label.size = 3)+
  ggtitle("Supervised sub cell type labeling by HPCA")+
  theme(text = element_text(size=10),							
        plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_sub1~.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

save(MCL,file="./data/MCL_20181107.RData")
##############################
# subset Seurat
###############################
table(MCL@meta.data$orig.ident)
MCL@meta.data$orig.ident = gsub("Pt-MD","MD",MCL@meta.data$orig.ident)
table(MCL@ident)

df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
tests <- paste0("test",2:4)
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
