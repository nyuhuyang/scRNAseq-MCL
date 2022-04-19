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
library(gplots)
library(cowplot)
library(eulerr)
library(openxlsx)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

path <- "Yang/PALIBR/51_samples/Fig. 3/"
if(!dir.exists(path)) dir.create(path, recursive = T)
#============= after running Differential_analysis.R Rsscript ===========
# X4clusters =====
choose = "X4cluster"
(de_file_names = list.files("output/20210802_4clusters",pattern = "MCL_only_51-FC0_cluster"))
idents.all = 1:4
genes.de <- list()
for(i in seq_along(idents.all)){
        genes.de[[i]] <- read.csv(paste0("output/20210802_4clusters/",de_file_names[i]),
                                  stringsAsFactors = F,row.names = 1)
        genes.de[[i]] <- genes.de[[i]][order(genes.de[[i]]$p_val, -genes.de[[i]][, 2]), ]
        genes.de[[i]]$cluster <- idents.all[i]
        genes.de[[i]]$gene <- rownames(x = genes.de[[i]])
}
gde.markers <- bind_rows(genes.de)
rownames(x = gde.markers) <- make.unique(names = as.character(x = gde.markers$gene))
write.csv(gde.markers, file = paste0(path,choose,"/MCL_only_",choose,"_51-FC0.csv"))

#============
choose = "X4cluster_vs_Normal"
(de_file_names = list.files("output/20200226",pattern = "MCL_B_41-FC0.25_"))
idents.all = c("B_cells",paste0("C",1:4))
genes.de <- list()
for(i in seq_along(idents.all)){
        genes.de[[i]] <- read.csv(paste0("output/20200226/",de_file_names[i]),
                                  stringsAsFactors = F,row.names = 1)
        genes.de[[i]] <- genes.de[[i]][order(genes.de[[i]]$p_val, -genes.de[[i]][, 2]), ]

        genes.de[[i]]$cluster <- idents.all[i]
        genes.de[[i]]$gene <- rownames(x = genes.de[[i]])
}
gde.all <- bind_rows(genes.de)
rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
write.csv(gde.all, file = paste0(path,"Figure Sources/",choose,"/",choose,"_41-FC0.25.csv"))
#============
choose = "X4cluster_vs_B"
(de_file_names = list.files("output/20200227",pattern = "MCL_B_41-FC0_"))
idents.all = paste0("C",1:4)
genes.de <- list()
for(i in seq_along(idents.all)){
        genes.de[[i]] <- read.csv(paste0("output/20200227/",de_file_names[i]),
                                  stringsAsFactors = F,row.names = 1)
        genes.de[[i]] <- genes.de[[i]][order(genes.de[[i]]$p_val, -genes.de[[i]][, 2]), ]

        genes.de[[i]]$cluster <- idents.all[i]
        genes.de[[i]]$gene <- rownames(x = genes.de[[i]])
}
gde.all <- bind_rows(genes.de)
rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
write.csv(gde.all, file = paste0(path,choose,"/",choose,"_41-FC0.csv"))
#=================== Doheatmap =================
object = readRDS(file = "data/MCL_SCT_51_20210724.rds")
DefaultAssay(object) = "SCT"

object = subset(object, subset =  Doublets == "Singlet"
                & X4cluster  %in% c("1","2","3","4")
)
object = subset(object, subset =  orig.ident %in% c("Pt2_30Pd","Pt3_BMA72_6","Pt3_72_6",
                                                    "N01","N02","N03","N04"), invert = T)

Idents(object) = "X4cluster"
table(Idents(object))

object$X4cluster %<>% factor(levels=c("1","2","3","4"))
Top_n = 40
markers <- c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11")
markers = markers[markers %in% rownames(object)]
# remove mito
(mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
# remove IGL, IGK, IGK
(IG.gene <- grep(pattern = "^IG", x = gde.markers$gene, value = T))
(IG.gene <- grep(pattern = "^IGL|^IGK|^IGH", x = gde.markers$gene))
if(length(IG.gene)>0) gde.markers = gde.markers[-IG.gene,]
GC()

top = gde.markers %>% group_by(cluster) %>% top_n(Top_n, avg_log2FC)

features = c(as.character(top$gene),
             tail(VariableFeatures(object = object), 2),
             markers)
object %<>% ScaleData(features=features)
featuresNum <- make.unique(features, sep = ".")
object %<>% MakeUniqueGenes(features = features)

DoHeatmap.1(object =object, features = featuresNum, Top_n = Top_n,
            do.print=T, angle = 0, group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0, label=F, cex.row= ifelse(i==2,4,2), legend.size = 0,width=10, height=6.5,
            pal_gsea = FALSE,
            save.path = paste0(path,choose,"/"),file.name = "Heatmap_top40_4cluster.jpeg",
            title = paste("Top",Top_n,"DE genes in 4, B/MCL cells cluster"))

DoHeatmap.2(object =object, group.by = c("X4cluster","orig.ident"),
            features = featuresNum, Top_n = Top_n,
            do.print=T, angle = 45, group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.02, label=T, cex.row= ifelse(i==2,4,2), legend.size = 0,width=10, height=6.5,
            pal_gsea = FALSE,
            save.path = paste0(path,choose,"/"),file.name = "Heatmap_top40_4cluster_legend.jpeg",
            title = paste("Top",Top_n,"DE genes in 4, B/MCL cells cluster"))

# ====================
path <- "Yang/Figure Sources/Log2UMI/"
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
# prepare average expression
# load data
object = readRDS(file = "data/MCL_SCT_51_20210724.rds")
DefaultAssay(object) = "SCT"

B_cells_MCL = subset(object, subset =  Doublets == "Singlet"
                     & X4cluster  %in% c("1","2","3","4")
)
B_cells_MCL$X4_orig.ident = paste(B_cells_MCL$orig.ident,
                                  B_cells_MCL$X4cluster, sep = "_C")
table(B_cells_MCL@meta.data$X4_orig.ident)
Idents(B_cells_MCL) = "X4_orig.ident"

B_cells_MCL_exp <- AverageExpression(B_cells_MCL,assays = "SCT")
exp = (log2(B_cells_MCL_exp$SCT + 1))
write.csv(exp,paste0(path,"B_MCL_log2UMI.csv"))

B_cells_MCL_number = table(B_cells_MCL@meta.data$X4_orig.ident) %>%
        as.data.frame() %>% t
rownames(B_cells_MCL_number) = c("samples","cell.number")
write.csv(B_cells_MCL_number,paste0(path,"B_MCL_cells_number.csv"))

# Doheatmap for MCL longitudinal X5 clusters ================
(load(file = "data/B_cells_MCL_43_20190917.Rda"))
df_samples <- readxl::read_excel("doc/191030_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:12)))
markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
df_samples = df_samples[sample_n,]

table(B_cells_MCL@meta.data$orig.ident)

B_cells_MCL@meta.data$orig.ident %<>% plyr::mapvalues(from = unique(df_samples$sample),
                                                      to = unique(df_samples$publication.id))
NewNames = paste0(B_cells_MCL@meta.data$orig.ident,"_",B_cells_MCL@meta.data$Barcode)
B_cells_MCL %<>% RenameCells(new.names = NewNames)
rownames(B_cells_MCL@reductions$tsne@cell.embeddings) = colnames(B_cells_MCL)
gsub("_.*","",rownames(B_cells_MCL@reductions$tsne@cell.embeddings)) %>% table

Idents(B_cells_MCL) = "groups"
#groups = c("Untreated","Pt-11","Pt-17","AFT-03","AFT-04","Pt-AA13","Pt-25","Pt-27")
clusters =  3:4
groups = c("Untreated","Pt-17","Pt-25")
for(i in 2:3){

        subset_MCL <- subset(B_cells_MCL, idents = groups[i])
        Idents(subset_MCL) = "X5clusters"
        for(k in clusters){
                subset.MCL <- subset(subset_MCL,idents = k)

                print(samples <- unique(subset.MCL$orig.ident))
                df <- df_samples[df_samples$publication.id %in% samples,]
                subset.MCL@meta.data$orig.ident = factor(subset.MCL@meta.data$orig.ident,
                                                         levels = rev(df$publication.id))
                Idents(subset.MCL) = "orig.ident"
                gde.markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.3, only.pos = T,
                                                  test.use = "MAST")
                write.csv(gde.markers,paste0(path,groups[i],"_Clusters_",k,
                                             "_FC0.3_markers.csv"))
                (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
                if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
                GC()
                #DoHeatmap.1======
                Top_n = 40
                top = gde.markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)

                features = c(as.character(top$gene),
                             tail(VariableFeatures(object = subset.MCL), 2),
                             markers)
                subset.MCL %<>% ScaleData(features=features)
                featuresNum <- make.unique(features, sep = ".")
                subset.MCL %<>% MakeUniqueGenes(features = features)
                DoHeatmap.1(object =subset.MCL, features = featuresNum, Top_n = Top_n,
                            do.print=T, angle = 0, group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
                            group.bar.height = 0, label=F, cex.row= ifelse(i==2,4,2), legend.size = 0,width=10, height=6.5,
                            pal_gsea = FALSE,
                            title = paste("Top",Top_n,"DE genes in longitudinal",groups[i],
                                          "B/MCL cells cluster",k))

                # rename file
                v <- UniqueName(object = subset.MCL, fileName = "subset.MCL",unique.name = T)
                v = paste0(v,"_",FindIdentLabel(object))
                old.name = paste0(path,"Heatmap_top",Top_n,"_",v,"_Legend.jpeg")
                file.rename(old.name, paste0(path,"Heatmap_top",Top_n,"_",groups[i],
                                             "_cluster",k, "_rev.jpeg"))
                remove(subset.MCL);GC()
        }
        remove(subset_MCL);GC()
}


# Doheatmap for MCL longitudinal Normal vs MCL ================
df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) =  tolower(colnames(df_samples))
df_samples[df_samples$sample %in% "MD","tsne"] = 0
df_samples[df_samples$sample %in% "MD","sample"] = "Normal"

Idents(B_cells_MCL) = "groups"
#groups = c("Untreated","Pt-11","Pt-17","AFT-03","AFT-04","Pt-AA13","Pt-25","Pt-27")
groups = c("Untreated","Pt-17","Pt-25")
for(i in 1:length(groups)){

        subset.MCL <- subset(B_cells_MCL, idents = c("Normal",groups[i]))

        (samples = unique(subset.MCL$orig.ident))
        df = df_samples[df_samples$sample %in% samples,]
        subset.MCL@meta.data$orig.ident = factor(subset.MCL@meta.data$orig.ident,
                                                 levels = df$sample[order(df$tsne)])
        Idents(subset.MCL) %<>% factor()
        Idents(subset.MCL) = "orig.ident"
        samples = samples[-which(samples %in% "Normal")]
        gde.markers_list <- list()
        for(k in 1:length(samples)){
                gde.markers_list[[k]] <- FindMarkers.UMI(subset.MCL,
                                                         ident.1 = samples[k],
                                                         ident.2 = "Normal",
                                                         logfc.threshold = 0.1, only.pos = F,
                                                         test.use = "MAST")
                gde.markers_list[[k]]$cluster <- samples[k]
                gde.markers_list[[k]]$gene <- rownames(x = gde.markers_list[[k]])
        }
        gde.markers <- do.call(rbind, gde.markers_list)
        write.csv(gde.markers,paste0(path,"B/B_MCL_DE/B_",groups[i],".csv"))

        gde.markers = read.csv(file = paste0("output/20190622/B/B_MCL_DE/B_",groups[i],".csv"))
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        subset.MCL %<>% ScaleData(features= unique(gde.markers$gene))
        Top_n = 20
        DoHeatmap.1(subset.MCL, marker_df = gde.markers, Top_n = Top_n, do.print=T, angle = 0,
                    group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
                    label=T, cex.row=5, legend.size = NULL,width=10, height=7,unique.name = T,
                    title = paste("Top",Top_n,"DE genes of",groups[i],
                                  "MCL cells over Normal B cells"))
}

# venn diagram
path <- "Yang/Figure Sources/51_samples/"
read.path <- "Yang/Figure Sources/51_samples/heatmaps_full/"
save.path <- paste0(path,"VennDiagram/")
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

opts <- data.frame("choose" = rep(c("X4cluster","X4cluster_vs_Normal"), each = 2),
                  "LogFC" = c(0,0.1,0,0.25), stringsAsFactors = F)

for(i in 1){
        choose = opts$choose[i]
        value = opts$LogFC[i]
        gde.all <- read.csv(file = paste0(read.path,choose,"/","_41-FC",value,".csv"))
        gde.all <- gde.all[gde.all$avg_logFC > 0 ,]
        pos_genes <- eulerr(gde.all,shape =  "circle",#key = c("C1","C2","C3","C4","B_cells"),
               cut_off = "p_val_adj", cut_off_value = 0.01,do.print = T,return.raw = T,
               save.path = paste0(save.path, choose,"_FC",value,"_"))

}

euler_df <- eulerr::euler(pos_genes,shape = "circle")
pos_genes_list <- as.list(euler_df$original.values)
names(pos_genes_list) %<>% paste("=",pos_genes_list)
id <- eulerr:::bit_indexr(4)

for (i in nrow(id):1) {
        pos_genes_list[[i]] = Reduce(intersect, pos_genes[id[i,]])  %>%
                setdiff(Reduce(union, pos_genes[!id[i,]]))
}
pos_genes_df <- list2df(pos_genes_list)
write.xlsx(pos_genes_df, asTable = F,
           file = paste0(save.path,"postive_shared_gene_list_",choose,"-FC",value,".xlsx"),
           borders = "surrounding")


#==============
choose = "X4clustes"
gde.all <- read.csv(file = paste0(read.path,choose,"/",choose,"_41-FC0.csv"))
gde.all <- gde.all[gde.all$avg_logFC > 0 ,]

eulerr(gde.all,shape =  "ellipse",key = c("C1","C2"),
       cut_off = "p_val_adj", cut_off_value = 0.01,do.print = T,
       save.path = paste0(path, "C1_C2"))
eulerr(gde.all,shape =  "ellipse",key = c("C1","C2","C3"),
             cut_off = "p_val_adj", cut_off_value = 0.01,do.print = T,
       save.path = paste0(path, "C1_C2_C3"))
eulerr(gde.all,shape =  "ellipse",key = c("C1","C2","C3","C4"),
             cut_off = "p_val_adj", cut_off_value = 0.01,do.print = T,
       save.path = paste0(path, "C1_C2_C3_C4"))

# load data
path <- "Yang/Normal_B/"
object = readRDS(file = "data/MCL_41_B_20200225.rds")
Idents(object) = "orig.ident"
object <- BuildClusterTree(object)
jpeg(paste0(path,"PlotClusterTree_all_B_MCL.jpeg"), units="in", width=10, height=10,res=600)
PlotClusterTree(object)
dev.off()

Normal <- c("N01","N02","N03","N04","Pt25_1","Pt25_24","Pt25_25Pd")
Normal <- subset(object, idents = Normal)
Normal$X4clusters_orig.ident = paste0(Normal$X4clusters,"_",Normal$orig.ident)
Idents(Normal) = "X4clusters_orig.ident"
table(Idents(Normal))
Normal %<>% BuildClusterTree
jpeg(paste0(path,"PlotClusterTree_B.jpeg"), units="in", width=10, height=10,res=600)
PlotClusterTree(Normal)
dev.off()

Idents(Normal) = "X4clusters"
Normal %<>% subset(idents = c("C1","C2"))
Idents(Normal) = "X4clusters_orig.ident"
table(Idents(Normal))
Normal %<>% BuildClusterTree
jpeg(paste0(path,"PlotClusterTree_B_C1_2.jpeg"), units="in", width=10, height=10,res=600)
PlotClusterTree(Normal)
dev.off()

# volcano plots =======================
path <- "Yang/Normal_B/"
save.path <- paste(path, "VolcanoPlots/p_val/")
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

DE_files = c("N01-N02-N03-N04","N01", "N02","N03","N04","Pt25_1","Pt25_24-Pt25_1",
             "Pt25_24","Pt25_25Pd")
for(i in seq_along(DE_files)){
        B_markers = read.csv(paste0(path,"DE files/B_41-FC0_",DE_files[i],".csv"),row.names = 1)
        if(DE_files[i] != "Pt25_24-Pt25_1") B_markers = B_markers[B_markers$cluster %in% "C2",]
        if(DE_files[i] == "Pt25_24-Pt25_1") B_markers = B_markers[B_markers$cluster %in% "C2_Pt25_24",]

        write.csv(B_markers,paste0(path,"DE files/B_41-FC0_",DE_files[i],".csv"))

        g <- VolcanoPlots(B_markers, cut_off_value = 0.05, cut_off = "p_val", cut_off_logFC = 0,top = 20,
                                 cols = c("#2a52be","#d2dae2","#d9321f"),alpha=1, size=2,
                                 legend.size = 12)+ theme(legend.position="bottom")
        if(DE_files[i] != "Pt25_24-Pt25_1") g = g + ggtitle(paste("Cluster 2 / cluster 1 in", DE_files[i]))
        if(DE_files[i] == "Pt25_24-Pt25_1") g = g + ggtitle("Cluster 2 of Pt25_C24 / cluster 2 of Pt25_C1")
        g = g + TitleCenter()#+theme_bw()

        jpeg(paste0(save.path,"VolcanoPlots_",DE_files[i],".jpeg"), units="in", width=10, height=10,res=600)
        print(g)
        dev.off()
        Progress(i, length(DE_files))
}

# venn diagram gene list 20200512 =======================
set.seed(101)
path <- "Yang/PALIBR Figures legends methods/Figure 2/"
save.path <- "Yang/PALIBR Figures legends methods/Figure S3/"
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

choose = "X4clusters"
read.path <- paste0(path, "Figure Sources/",choose,"/")

B_markers <- read.csv(paste0(read.path, "X4clusters_41-FC0.1.csv"))
B_markers <- B_markers[B_markers$avg_logFC > 0 ,]
eulerr(B_markers,shape =  "circle",key = c("C1","C2","C3","C4"),
       cut_off = "p_val_adj", cut_off_value = 0.01,do.print = T,return.raw = F,do.return = T,
       save.path = paste0(save.path, choose,"_"))

pos_genes <- eulerr(B_markers,shape =  "circle", key = c("C1","C2","C3","C4"),
       cut_off = "p_val_adj", cut_off_value = 0.01,do.print = F,return.raw = T,
       save.path = save.path)
euler_df <- eulerr::euler(pos_genes,shape = "circle")
pos_genes_list <- as.list(euler_df$original.values)
names(pos_genes_list) %<>% paste("=",pos_genes_list)
id <- eulerr:::bit_indexr(4)

for (i in nrow(id):1) {
        pos_genes_list[[i]] = Reduce(intersect, pos_genes[id[i,]])  %>%
                setdiff(Reduce(union, pos_genes[!id[i,]]))
}
pos_genes_df <- list2df(pos_genes_list)
write.xlsx(pos_genes_df, asTable = F,
           file = paste0(save.path,"postive_shared_gene_list_",choose,".xlsx"),
           borders = "surrounding")

# choose == "MCL_vs_B_cells" ======================
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# load data
object = readRDS(file = "data/MCL_41_B_20200225.rds")
DefaultAssay(object) = "SCT"
Idents(object) = "orig.ident"
samples = as.character(unique(object$orig.ident))
opts = data.frame(ident.1 = samples[2:length(samples)],
                  ident.2 = rep("N01", length(samples)-1),
                  stringsAsFactors = F)
for(i in 16:nrow(opts)){
        (opt = opts[i,])
        sub_object <- subset(object,idents = c(opt$ident.1,opt$ident.2))
        save.path = paste0(path, "MCL_B_", opt$ident.1, "-",opt$ident.2,"/")
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        file.copy(from = paste0(path,"MCL_B_",opt$ident.1, "-",opt$ident.2,".csv"),
                  to = paste0(save.path,basename(save.path),".csv"))

        MCL_markers = read.csv(file= paste0(save.path,basename(save.path),".csv"),
                             row.names = 1, stringsAsFactors=F)
        table(MCL_markers$cluster)
        markers <- FilterGenes(sub_object,c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B","GZMK"))
        (MT_gene <- grep("^MT-",MCL_markers$gene))
        if(length(MT_gene) >0 ) MCL_markers = MCL_markers[-MT_gene,]
        Top_n = 40

        top = MCL_markers %>% group_by(cluster) %>%
                top_n(40, avg_logFC)
        unique(top$cluster)
        #top = top[order(top$cluster),]
        write.csv(top,paste0(save.path,"Top40_","MCL_B_",
                             opt$ident.1, "-",opt$ident.2,".csv"))
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = sub_object), 2),
                     markers)
        #DoHeatmap.1======
        # raw heatmap
        featuresNum <- make.unique(features, sep = ".")
        exp = AverageExpression(sub_object[features,],
                                assays = "SCT") %>% .$SCT
        exp %<>% MakeUniqueGenes(features = features)
        exp[tail(VariableFeatures(object = sub_object), 2),] =0

        (group.by = c(opt$ident.2, opt$ident.1))
        DoHeatmap.matrix(exp, features = featuresNum,
                         group.by = group.by,
                         size = 6,angle = 0,label =F,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         width=2.5, height=10,res=600,no.legend = F,
                         cex.row=5,
                         group.colors = gg_color_hue(2),
                         do.print = T,
                         unique.name = "cell.types",
                         save.path = paste0(save.path,"Heatmap_top40_","MCL_B_",
                                            opt$ident.1, "-",opt$ident.2,"_raw"))
        # scale heatmap
        sub_object %<>% ScaleData(features = features)
        scale_exp = AverageExpression(sub_object[features,],
                                      assays = "SCT", slot = "scale.data") %>% .$SCT
        scale_exp %<>% MakeUniqueGenes(features = features)
        scale_exp[tail(VariableFeatures(object = sub_object), 2),] =0
        DoHeatmap.matrix(scale_exp, features = featuresNum,
                         group.by = group.by,
                         size = 6,angle = 0,label =F,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         width=2.5, height=10,res=600,no.legend = F,
                         cex.row=5,
                         group.colors = gg_color_hue(2),
                         do.print = T,
                         unique.name = "cell.types",
                         save.path = paste0(save.path,"Heatmap_top40_","MCL_B_",
                                            opt$ident.1, "-",opt$ident.2,"_scale"))
        MCL_markers = MCL_markers[MCL_markers$cluster %in% opt$ident.1,]
        #avg_logFC = MCL_markers[MCL_markers$cluster %in% opt$ident.2,"avg_logFC"]
        #MCL_markers[MCL_markers$cluster %in% opt$ident.2,"avg_logFC"] = avg_logFC * -1
        p <- VolcanoPlots(data = MCL_markers, cut_off_value = 0.05, cut_off = "p_val", cut_off_logFC = 0.1,
                          top = 20, cols = c("#2a52be","#d2dae2","#d9321f"),alpha=1, size=2,
                          legend.size = 12)+
                ggtitle(paste0(opt$ident.1, " \\ ",opt$ident.2, " in ", "MCL and B cells"))+
                theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"),
                      legend.position="bottom")
        jpeg(paste0(save.path,"VolcanoPlots_","MCL_B_",
                    opt$ident.1, "-",opt$ident.2,".jpeg"), units="in", width=10, height=7,res=600)
        print(p)
        dev.off()
        Progress(i, nrow(opts))
}


(load(file = "data/B_cells_MCL_43_20190917.Rda"))
Idents(B_cells_MCL) = "orig.ident"

sub_object = subset(B_cells_MCL, idents = c(grep("Pt-25",unique(Idents(B_cells_MCL)),value = T)))
sub_object = subset(B_cells_MCL, idents = "Pt-25-C25")

Idents(sub_object) = "manual"
TSNEPlot.1(sub_object,split.by ="orig.ident")


monocle2_Pt25_25Pd = readRDS(file = "data/monocle2_Pt25_25Pd_DE.rds")
