library(SingleR)
library(Seurat)
library(magrittr)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
source("../R/Seurat3_functions.R")
source("../R/SingleR_functions.R")
source("R/util.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file = "data/MCL_V3_Harmony_43_20190627.Rda"))
(load(file="output/singlerT_MCL_43_20190916.Rda"))

# if singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(object)
if(length(singler$singler[[1]]$SingleR.single$labels) < ncol(object)){
        all.cell = colnames(object);length(all.cell)
        know.cell = names(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = subset(object, cells = know.cell)
}

table(rownames(singler$singler[[1]]$SingleR.single$labels) == colnames(object))
singler$meta.data$orig.ident = object@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = object@reductions$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Idents(object) # the Seurat clusters (if 'clusters' not provided)
save(singler,file="output/singlerT_MCL_43_20190916.Rda")
##############################
# add singleR label to Seurat
###############################

singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       "singler1main" = singler$singler[[1]]$SingleR.single.main$labels,
                       "orig.ident" = gsub("\\_.*","",rownames(singler$singler[[1]]$SingleR.single$labels)),
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))
head(singlerDF)
singlerDF = singlerDF[colnames(object),]
table(rownames(singlerDF) == colnames(object))
apply(singlerDF,2,function(x) length(unique(x)))

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"DrawHeatmap_sub1_N.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,normalize = T))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

table(singlerDF$singler1sub, singlerDF$orig.ident) %>% kable %>%
        kable_styling()
##############################
# adjust cell label
##############################
# combine cell types
singlerDF$singler1sub = gsub("MCL:.*","MCL",singlerDF$singler1sub)
singlerDF$singler1sub = gsub("B_cells:PB","B_cells:Plasma_cells",singlerDF$singler1sub)
singlerDF$manual = gsub("B_cells:.*","B_cells",singlerDF$singler1sub)
singlerDF$manual = gsub("MEP|CLP|HSC|CMP|GMP|MPP","HSC/progenitors",singlerDF$manual)
singlerDF$manual = gsub("T_cells:CD4\\+_.*","T_cells:CD4+",singlerDF$manual)
singlerDF$manual = gsub("T_cells:CD8\\+_.*","T_cells:CD8+",singlerDF$manual)
singlerDF$manual = gsub("T_cells:Tregs","T_cells:CD4+",singlerDF$manual)
singlerDF$manual = gsub("DC|Macrophages|Macrophages:M1","Myeloid cells",singlerDF$manual)
singlerDF$manual = gsub("Erythrocytes","Myeloid cells",singlerDF$manual)
singlerDF$manual = gsub("Eosinophils|Megakaryocytes|Monocytes","Myeloid cells",singlerDF$manual)
singlerDF$manual = gsub("Adipocytes|Fibroblasts|mv_Endothelial_cells","Nonhematopoietic cells",singlerDF$manual)
table(singlerDF$manual, singlerDF$orig.ident) %>% kable() %>% kable_styling()

# reduce false positive results (B cells are labeled as MCL in normal samples)
# and false negative results (MCL cells are labeled as B cells in MCL samples)
# singler1sub false negative results  =========
CCND1 = FetchData(object,"CCND1")
singlerDF$CCND1 = CCND1$CCND1
singlerDF[(singlerDF$CCND1 >0 & singlerDF$manual %in% "B_cells"),"manual"] = "MCL"
# manual false positive results  ========
table(singlerDF$manual, object@meta.data$orig.ident) %>% kable %>% kable_styling()
normal_cells <- object$orig.ident %in% c("BH","DJ","MD","NZ") %>% rownames(singlerDF)[.]
singlerDF[normal_cells,"manual"] = gsub("MCL","B_cells",
                                              singlerDF[normal_cells,"manual"])
# singler1sub false positive results  =========
table(singlerDF$singler1sub, object$orig.ident) %>% kable %>% kable_styling()
singlerDF[normal_cells,"singler1sub"] = gsub("MCL:.*$","B_cells:Memory",
                                              singlerDF[normal_cells,"singler1sub"])
table(singlerDF$manual, object$orig.ident) %>% kable %>% kable_styling()
table(singlerDF$singler1sub, object$orig.ident)%>% kable %>% kable_styling()
table(singlerDF$manual %>% sort)%>% kable %>% kable_styling()
table(singlerDF$manual) %>% kable() %>% kable_styling()
singlerDF = singlerDF[, -5]
##############################
# process color scheme
##############################
singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
singler_colors2[duplicated(singler_colors2)]
length(singler_colors1);length(singler_colors2)
apply(singlerDF[,c("singler1sub","manual")],2,function(x) length(unique(x)))
singlerDF[,c("singler1sub")] %>% table() %>% kable() %>% kable_styling()
(remove = c(grep("singler1sub",colnames(object@meta.data))[1], ncol(object@meta.data)))
object@meta.data = object@meta.data[,-c(remove[1]:remove[2])]
object <- AddMetaData(object = object,metadata = singlerDF)
object <- AddMetaColor(object = object, label= "singler1sub", colors = singler_colors1)
object <- AddMetaColor(object = object, label= "manual", colors = singler_colors2)
Idents(object) <- "manual"

object %<>% sortIdent
TSNEPlot.1(object,cols = ExtractMetaColor(object), label = T, label.repel = T) + NoLegend()
save(object,file="data/MCL_V3_Harmony_43_20190627.Rda")

##############################
# draw tsne plot
##############################
object <- subset(object,idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)
Idents(object) = "manual"
object %<>% sortIdent()
table(Idents(object))
TSNEPlot.1(object = object, label = F, group.by = "manual",
       cols = ExtractMetaColor(object),no.legend = F,
       pt.size = 0.1,label.size = 3, do.print = T,do.return = F,
       title = "Cell type labeling by Blueprint + Encode + MCL")


cell_Freq <- table(Idents(object)) %>% as.data.frame
cell_Freq = cell_Freq[order(cell_Freq$Var1),]
cell_Freq$col = ExtractMetaColor(object)
cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=10, height=7,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "Cell_Type",xlab = "",
          palette = cell_Freq$col,x.text.angle = 90,
          title = "Numbers of major cell types in total 43 samples")+NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,size=15))
dev.off()
     
##############################
# subset Seurat
###############################
(load(file = "data/MCL_V3_Harmony_43_20190627.Rda"))

Idents(object) <-  "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "manual"
object %<>% subset(idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)

table(Idents(object))
object@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",object@meta.data$orig.ident)
Idents(object) = "orig.ident"

df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
tests <- paste0("test",c(2:12))
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df <- as.data.frame(df_samples[sample_n,])
        samples <- unique(df$sample)
        rownames(df) = samples
  
        samples <- c(ifelse(length(samples)>5,NA,"Normal"),df$sample[order(df$tsne)])
        print(samples <- samples[!is.na(samples)])
  
        subset_object <- subset(object, idents = samples)
        subset_object$orig.ident %<>% factor(levels = samples)
  
        Idents(subset_object) = "manual"
  
        subset_object %<>% sortIdent()
        TSNEPlot.1(subset_object, pt.size =0.3, 
                   strip.text.size = min(240/max(stringr::str_length(samples)),30),
                   group.by = "manual",split.by = "orig.ident",
                   cols = ExtractMetaColor(subset_object), ncol = length(samples),
                   unique.name = T, do.print = T,do.return = F,border = T,
                   width=length(samples)*2+2, height=3)
}

subset_object <- subset(object, idents = c("Pt-U01","Pt-U02","Pt-U03","Pt-U04"))
Untreated_exp <- AverageExpression(subset_object)
write.csv(Untreated_exp,file=paste0(path,"Untreated_UMI.csv"))
table(subset_object$manual, subset_object$orig.ident)%>% kable %>% kable_styling()

###################################
# creat seurat object
###################################
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


# correct orig.ident using barcode ==============
meta.data <- object@meta.data
head(meta.data)
for(sample in df_samples$sample){
        idx <- grep(paste0("^",sample,"_"),rownames(object@meta.data))
        if(length(idx)>0) meta.data[idx,"orig.ident"] = sample
}
object@meta.data = meta.data