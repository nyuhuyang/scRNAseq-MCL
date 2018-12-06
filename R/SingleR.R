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
(load(file = "./data/MCL_Harmony_23_20181205.Rda"))
MCL@scale.data = NULL; GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();
(load(file = 'data/ref_MCL_blue_encode.RData'))

singler = CreateSinglerObject(MCL@data, annot = NULL, project.name=MCL@project.name,
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              ref.list = list(ref),normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL)
# if singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(MCL@data)
all.cell = MCL@cell.names;length(all.cell)
know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
MCL = SubsetData(MCL, cells.use = know.cell)
MCL
singler$meta.data$orig.ident = MCL@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = MCL@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = MCL@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_MCL_11T_20181128.RData")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file = "./output/singler_MCL_11T_20181121.RData"))
(load(file = "../SingleR/data/Hpca.RData"))

##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub"=singler$singler[[1]]$SingleR.single$labels,
                       "singler1main"=singler$singler[[1]]$SingleR.single.main$labels,
                       "singler2sub"=singler$singler[[2]]$SingleR.single$labels,
                       "singler2main"=singler$singler[[2]]$SingleR.single.main$labels,
                       "kang" = FineTune(singler$other, main.type = FALSE),
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))

table(rownames(singlerDF) %in% MCL@cell.names)

#knowDF = data.frame("cell.names"= MCL@cell.names)
#ident.DF = full_join(singlerDF,knowDF, by="cell.names")
#ident.DF<- apply(ident.DF,2,as.character)
#rownames(ident.DF) = ident.DF[,"cell.names"]
#ident.DF = ident.DF[,-which(colnames(ident.DF) == "cell.names")]
apply(singlerDF,2,function(x) length(unique(x)))
#ident.DF[is.na(ident.DF)] <- "unknown"
MCL <- AddMetaData(object = MCL,
                   metadata = singlerDF)
MCL <- SetAllIdent(object = MCL, id = "singler2sub")

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"/DrawHeatmap_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler_T$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"/DrawHeatmap_sub2_N.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n = 50,normalize = T))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$singler[[2]]$SingleR.single$labels, singler$meta.data$orig.ident)) %>%
        kable_styling()
MCL@meta.data$singler1sub %>% table() %>% kable() %>% kable_styling()
MCL@meta.data$singler2sub %>% table() %>% kable() %>% kable_styling()

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
MCL <- AddMetaColor(object = MCL, label= "singler1sub", colors = singler_colors1[1:42])
MCL <- AddMetaColor(object = MCL, label= "singler2sub", colors = singler_colors2[1:35])
MCL <- SetAllIdent(object = MCL, id = "singler2sub")
TSNEPlot.1(MCL, colors.use = ExtractMetaColor(MCL),no.legend = F)

write.csv(MCL@meta.data,file = "output/MCL_metadata_20181107.csv")

##############################
# draw tsne plot
##############################
p3 <- TSNEPlot.1(object = MCL, do.label = T, group.by = "ident", 
                 do.return = TRUE, no.legend = T, 
                 colors.use = ExtractMetaColor(MCL),
                 pt.size = 1,label.size = 3,force = 2)+
  ggtitle("Supervised cell type labeling by Blueprint + Encode + MCL")+
  theme(text = element_text(size=10),							
        plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

save(MCL,file="./data/MCL_Harmony_20181127.Rda")
##############################
# subset Seurat
###############################
table(MCL@meta.data$orig.ident)
table(MCL@ident)

df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
tests <- paste0("test",c(1,3:4))
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- unique(df_samples$samples[sample_n])
        print(samples)
        
        cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(MCL, cells.use = cell.use)
        
        g <- SplitTSNEPlot(subset.MCL,group.by = "ident",split.by = "orig.ident",
                           #select.plots = c(1,5,4,3,2),
                           no.legend = T,do.label =F,label.size=3,
                           return.plots =T, label.repel = T,force=2)
        jpeg(paste0(path,test,"_TSNEPlot.jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, g))
        dev.off()
}

###################################
# creat seurat object
(load(file = 'data/ref_MCL_blue_encode.RData'))
MCL_blue_encode <- CreateSeuratObject(raw.data = ref$data) %>%
  NormalizeData %>% ScaleData

MCL_blue_encode@data = ref$data
ident = ref$main_types
ident[!grepl("B_cells|MCL",ident)] = "others"
names(ident) = MCL_blue_encode@cell.names
MCL_blue_encode@ident = factor(ident, levels = c("others","B_cells","MCL"))
MCL_B_markers <- FindAllMarkers.UMI(MCL_blue_encode, logfc.threshold = 0.1, only.pos = T,
                                   test.use = "MAST")
write.csv(MCL_B_markers, paste0(path,"MCL_B_markers.csv"))
g <- DoHeatmap.1(MCL_blue_encode, MCL_B_markers,Top_n = 50,
                 ident.use = paste("B cells vs. MCL vs. other cell types"),
                 group.label.rot = T,cex.row = 4,remove.key =T)
jpeg(paste0(path,"heatmap_B_MD_others.jpeg"), units="in", width=10, height=7,
     res=600)
print(g)
dev.off()
