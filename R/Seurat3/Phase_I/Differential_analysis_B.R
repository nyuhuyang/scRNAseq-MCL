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
source("../R/Seurat3_functions.R")

path <- "Yang/PALIBR Figures legends methods/Figure 2/"
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
path <- "Yang/PALIBR Figures legends methods/Figure 2/"
choose = "X4clusters"
save.path <- paste0(path, "Figure Sources/",choose,"/")
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
B_markers <- read.csv(paste0(save.path, "X4clusters_41-FC0.csv"))
B_markers <- B_markers[B_markers$avg_logFC > 0 ,]
eulerr(B_markers,shape =  "circle",key = c("C1","C2","C3","C4"),
       cut_off = "p_val_adj", cut_off_value = 0.01,do.print = T,return.raw = F,do.return = T,
       save.path = save.path)
       
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
opts = data.frame(ident.1 = c("PtU01","PtU02","PtU03","PtU04"),
                  ident.2 = c("N01", "N01", "N01","N01"),
                  stringsAsFactors = F)
for(i in 1:nrow(opts)){
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
        #MCL_markers = MCL_markers[MCL_markers$cluster %in% opt$ident.1,]
        avg_logFC = MCL_markers[MCL_markers$cluster %in% opt$ident.2,"avg_logFC"]
        MCL_markers[MCL_markers$cluster %in% opt$ident.2,"avg_logFC"] = avg_logFC * -1
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
