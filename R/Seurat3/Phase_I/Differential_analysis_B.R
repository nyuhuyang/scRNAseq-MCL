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
source("../R/Seurat3_functions.R")

path <- "Yang/Figure 2/Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)
#============= after running Differential_analysis.R Rsscript ===========
# X4clusters =====

choose = "X4clusters"
(de_file_names = list.files("output/20200225",pattern = "MCL_only_41-FC0_"))
idents.all = paste0("C",1:4)
genes.de <- list()
for(i in seq_along(idents.all)){
        genes.de[[i]] <- read.csv(paste0("output/20200225/",de_file_names[i]),
                                  stringsAsFactors = F,row.names = 1)
        genes.de[[i]] <- genes.de[[i]][order(genes.de[[i]]$p_val, -genes.de[[i]][, 2]), ]
        
        genes.de[[i]]$cluster <- idents.all[i]
        genes.de[[i]]$gene <- rownames(x = genes.de[[i]])
}
gde.all <- bind_rows(genes.de)
rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
write.csv(gde.all, file = paste0(path,choose,"/",choose,"_41-FC0.csv"))
#============
choose = "X4cluster_vs_Normal"
(de_file_names = list.files("output/20200226",pattern = "MCL_B_41-FC0_"))
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
write.csv(gde.all, file = paste0(path,choose,"/",choose,"_41-FC0.csv"))
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
#===================
path <- "Yang/Figure Sources/Log2UMI/"
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
# prepare average expression
# load data
B_cells_MCL = readRDS(file = "data/MCL_41_B_20200207.rds")

B_cells_MCL$orig.ident %<>% gsub("N02|N01|N03","Normal",.)
B_cells_MCL$X4_orig.ident = paste(B_cells_MCL$orig.ident,
                                  B_cells_MCL$X4clusters, sep = "_")
B_cells_MCL@meta.data$X4_orig.ident = gsub('^Normal_.*', 'Normal', B_cells_MCL@meta.data$X4_orig.ident)

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
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
save.path <- "Yang/Figure 2/Figure Sources/"
choose = "X4cluster_vs_Normal"
gde.all <- read.csv(file = paste0(save.path,choose,"/",choose,"_41-FC0.csv"))
gde.all <- gde.all[gde.all$avg_logFC > 0 ,]
eulerr(gde.all,shape =  "ellipse",key = c("C1","C2","C3","C4","B_cells"),
       cut_off = "p_val_adj", cut_off_value = 0.01,do.print = T,
       save.path = paste0(path, "All_"))

B_markers <- gde.all[gde.all$cluster %in% "B_cells",]
gde.all <- gde.all[!(gde.all$gene %in% B_markers$gene),]
gde.all %<>% rbind(B_markers)
pos.share_genes <- eulerr(gde.all,shape =  "ellipse",key = c("C1","C2","C3","C4","B_cells"),
       cut_off = "p_val_adj", cut_off_value = 0.01,do.print = T,
       save.path = paste0(path, "B_exclusive"), return.raw = T)
core_MCL_genes <- Reduce(intersect, pos.share_genes[2:4])
write.csv(core_MCL_genes, paste0(path, "core_MCL_genes.csv"))

core_MCL_genes_df <- gde.all[gde.all$gene %in% core_MCL_genes,] %>% group_by(cluster) %>%
        top_n(20, avg_logFC) %>%
        arrange(p_val_adj, 
                desc(avg_UMI.1)
        )
write.csv(core_MCL_genes_df, paste0(path, "core_MCL_genes_table.csv"))

#==============
choose = "X4clustes"
gde.all <- read.csv(file = paste0(save.path,choose,"/",choose,"_41-FC0.csv"))
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

