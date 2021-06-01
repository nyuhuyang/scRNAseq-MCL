########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
library(Seurat)
library(dplyr)
library(kableExtra)
library(magrittr)
library(ggplot2)
library(cowplot)
library(fgsea)
library(tibble)
library(ggsci)
library(openxlsx)
library(gplots)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

# load data
(load(file="data/MCL_41_harmony_20200225.Rda"))
df_samples <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
table(object$orig.ident)

(df_samples = df_samples[df_samples$`sample name` %in% object$orig.ident,])
(samples = df_samples$`sample name`[df_samples$`sample name` %in% object$orig.ident])

# preprocess
object$orig.ident %<>% factor(levels = samples)
Idents(object) = "orig.ident"
object %<>% subset(idents = "Pt2_30Pd", invert = T)
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "cell.types"
object %<>% subset(idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)
table(Idents(object))

#==== Figure 3-A ===========
path <- "Yang/PALIBR/Archive/Figure 3/Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)

object$cell.types <- plyr::mapvalues(object@meta.data$cell.types,
                                     from = c("B_cells",
                                              "Monocytes:CD14+",
                                              "Monocytes:CD16+",
                                              "NK_cells",
                                              "T_cells:CD4+",
                                              "T_cells:CD8+"),
                                             to = c("B",
                                                    "CD14 Monocytes",
                                                    "CD16 Monocytes",
                                                    "NK",
                                                    "CD4 T",
                                                    "CD8 T"))
object$cell.types %<>% as.character()
object$cell.types.colors = object$cell.types.colors
Idents(object) = "cell.types"
#bject %<>% sortIdent()
TSNEPlot.1(object = object, label = F, label.repel = F, group.by = "cell.types",
           cols = ExtractMetaColor(object),no.legend = F,border = T,
           pt.size = 0.1, do.print = T,do.return = F,legend.size = 25,
           title.size = 20,title = "tSNE plots for cell types of 41 samples",
           units= "in",width=9, height=7,hjust =0.5, save.path = path)
#==== Figure 3-B ===========
features <- FilterGenes(object,c("CD19","CCND1","SOX11",
                                 "CD3D","CD4","CD8A",
                                 "CD14","FCGR3A","FCGR1A",
                                 "GNLY","KLRC1","NCAM1"))
FeaturePlot.1(object,features = features, pt.size = 0.005, cols = c("gray90", "red"),
              alpha = 1,reduction = "tsne",
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 3,
              units = "in",width=9, height=12, no.legend = T,
              save.path = path, file.name = "FeaturePlot_label.jpeg")
FeaturePlot.1(object,features = features, pt.size = 0.005, cols = c("gray90", "red"),
              alpha = 1,reduction = "tsne",
              threshold = 1, text.size = 0, border = T,do.print = T, do.return = F,ncol = 3,
              units = "in",width=9, height=12, no.legend = T,
              save.path = path, file.name = "FeaturePlot_nolabel.jpeg")


#==== Figure 3-C ===========
path <- "Yang/PALIBR/Archive/Figure 3/Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)

B_cells_MCL = readRDS(file = "data/MCL_41_B_20200225.rds")
Idents(B_cells_MCL) = "orig.ident"
B_cells_MCL %<>% subset(idents = "Pt2_30Pd", invert = T)
markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
group.colors = c("#181ea4","#5f66ec","#f46072","#e6001c")
choose = c("X4clusters","X4cluster_vs_Normal","X4cluster_vs_B")[2]
filter_by = c("top40","logFC0.01")[1]
B_cells_MCL[["patient"]] = gsub("_.*","",B_cells_MCL$orig.ident)
B_cells_MCL[["samples"]] = droplevels(B_cells_MCL[["orig.ident"]])
Idents(B_cells_MCL) = "patient"
#B_cells_MCL %<>% subset(idents = "Pt25")
#cells <- sample(colnames(B_cells_MCL), size = 300)
#B_cells_MCL %<>% subset(cells = cells)
if(choose == "X4clusters"){
        B_cells_MCL %<>% subset(idents = c("N01","N02","N03","N04"), invert = T)
        B_cells_MCL@meta.data$samples %<>% droplevels()
        Idents(B_cells_MCL) = "X4clusters"

        B_cells_MCL %<>% sortIdent()
        Idents(B_cells_MCL) %<>% factor(levels = paste0("C",1:4))
        table(Idents(B_cells_MCL))
        X4clusters_markers = read.csv(file= paste0(path,choose,"/X4clusters_41-FC0.csv"),
                                      row.names = 1, stringsAsFactors=F)
        table(X4clusters_markers$cluster)
        X4clusters_markers$cluster %<>% factor(levels = paste0("C",1:4))
        markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
        (MT_gene <- grep("^MT-",X4clusters_markers$gene))
        X4clusters_markers = X4clusters_markers[-MT_gene,]
        Top_n = 40
        top = switch (filter_by,
                      "top40"= X4clusters_markers %>% group_by(cluster) %>%
                              top_n(Top_n, cluster) %>% top_n(Top_n, avg_logFC),
                      "logFC0.01" = X4clusters_markers %>% group_by(cluster) %>%
                              filter(avg_logFC < 0.05 & avg_logFC > 0.01)
        )
        table(top$cluster)
        top = top[order(top$cluster),]
        write.csv(top,paste0(path,choose,"/",filter_by,"_genes_heatmap.csv"))
        #write.csv(top,paste0(path,choose,"/avg_logFC0.01_genes_heatmap.csv"))

        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)

        B_cells_MCL %<>% ScaleData(features=features)
        featuresNum <- make.unique(features, sep = ".")
        B_cells_MCL %<>% MakeUniqueGenes(features = features)
        DoHeatmap.2(B_cells_MCL, features = featuresNum, #Top_n = Top_n,
                    do.print=T, do.return=F,
                    angle = 0,
                    group.by = c("X4clusters","tissues"),group.bar = T,
                    group.colors = group.colors,
                    group2.colors= Singler.colors,
                    title.size = 16, no.legend = F,size=5,hjust = 0.5,
                    group.bar.height = 0.02, label=F, cex.row= 5,
                    width=12, height=14,
                    pal_gsea = FALSE,
                    file.name = "Heatmap_",filter_by,"_X4clusters_tissue.jpeg",
                    title = "Top 40 DE genes in 4 B/MCL clusters",
                    save.path = paste0(path,choose))
}


if(choose == "X4cluster_vs_Normal"){
        B_cells_MCL$X4clusters_normal = as.character(B_cells_MCL$X4clusters)
        B_cells_MCL@meta.data[B_cells_MCL$orig.ident %in% c("N01","N02","N03"),
                                "X4clusters_normal"] = "Normal"
        Idents(B_cells_MCL) = "X4clusters_normal"
        #B_cells_MCL %<>% subset(idents = "B_cells", invert = T)
        B_cells_MCL %<>% sortIdent()
        table(Idents(B_cells_MCL))
        csv_files <- list.files("output/20210531",full.names = T,pattern = ".csv")
        deg_list <- lapply(csv_files,function(x) {
                tmp = read.csv(x,row.names = 1, stringsAsFactors=F)
                tmp$cluster = gsub(".*FC0.25_","",x) %>%
                        gsub("\\.csv","",.) %>%
                        gsub("[1-9]_","",.)
                tmp$gene = rownames(tmp)
                tmp
        })
        X4clusters_markers = bind_rows(deg_list)

        #X4clusters_markers = read.csv(file=paste0(path, choose,"/X4cluster_vs_Normal_41-FC0.csv"),
        #                               row.names = 1, stringsAsFactors=F)
        #X4clusters_markers %<>% .[!(.$cluster %in% "B_cells"),]
        table(X4clusters_markers$cluster)
        markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
        (MT_gene <- grep("^MT-",X4clusters_markers$gene))
        X4clusters_markers = X4clusters_markers[-MT_gene,]
        Top_n = 30
        top = X4clusters_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
        table(top$cluster)
        #top = top[top$cluster %in% c("C1","C2","C3","C4"),]
        write.csv(top,paste0("output/20210531/X4cluster_vs_Normal/","top40_4clusters_over_normal_genes_heatmap.csv"))
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        B_cells_MCL %<>% ScaleData(features=features)
        featuresNum <- make.unique(features, sep = ".")
        B_cells_MCL = MakeUniqueGenes(object = B_cells_MCL, features = features)

        B_cells_MCL$X4clusters_normal %<>% factor(levels = c("Normal", paste0("C",1:4)))
        Idents(B_cells_MCL) = "X4clusters_normal"
        table(Idents(B_cells_MCL))
        DoHeatmap.2(B_cells_MCL, features = featuresNum, Top_n = Top_n,
                    do.print=T, angle = 0,
                    group.by = c("X4clusters_normal","orig.ident"),group.bar = T,
                    group.colors = c("#31aa3a",group.colors,Singler.colors[1:4]),
                    title.size = 16, no.legend = F,size=5,hjust = 0.5,
                    group.bar.height = 0.02, label=F, cex.row= 5,
                    width=14, height=19,
                    pal_gsea = FALSE,
                    file.name = "Heatmap_top30_X4clusters_vs_normal.jpeg",
                    title = "Top 40 DE genes in 4 B/MCL clusters and healthy donors",
                    save.path = paste0(path,choose))
        }

if(choose == "X4cluster_vs_B"){
        Idents(B_cells_MCL) = "orig.ident"
        B_cells_MCL %<>% subset(idents = c("N01","N02","N03","N04"),invert = T)
        B_cells_MCL$X4clusters_B = as.character(B_cells_MCL$X4clusters)
        B_cells_MCL$X4clusters_B %<>% paste(B_cells_MCL$cell.types, sep = "_")
        B_cells_MCL$X4clusters_B %<>% gsub(".*_B_cells","B_cells",.)
        B_cells_MCL$X4clusters_B %<>% gsub("_MCL","",.)
        Idents(B_cells_MCL) = "X4clusters_B"
        B_cells_MCL %<>% sortIdent()
        table(Idents(B_cells_MCL))

        X4clusters_markers = read.csv(file=paste0(path,choose,"/",choose,"_41-FC0.csv"),
                                      row.names = 1, stringsAsFactors=F)
        table(X4clusters_markers$cluster)
        markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
        (MT_gene <- grep("^MT-",X4clusters_markers$gene))
        X4clusters_markers = X4clusters_markers[-MT_gene,]
        Top_n = 40
        top = X4clusters_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
        table(top$cluster)
        top = top[top$cluster %in% c("C1","C2","C3","C4"),]
        write.csv(top,paste0(path,choose,"/Fig3D_4clusters_over_normal_genes_heatmap.csv"))
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        B_cells_MCL %<>% ScaleData(features=features)
        featuresNum <- make.unique(features, sep = ".")
        B_cells_MCL = MakeUniqueGenes(object = B_cells_MCL, features = features)

        B_cells_MCL@meta.data[,"X4clusters_B"] %<>% factor(levels = c("B_cells", paste0("C",1:4)))
        Idents(B_cells_MCL) = "X4clusters_B"
        table(Idents(B_cells_MCL))
        DoHeatmap.2(B_cells_MCL, features = featuresNum, Top_n = Top_n,
                    do.print=T, angle = 0, group.bar = T,
                    group.by = c("X4clusters_B","samples"),
                    group.colors = c("#31aa3a",group.colors),
                    title.size = 16, no.legend = F,size=5,hjust = 0.5,
                    group.bar.height = 0.02, label=F, cex.row= 5,
                    width=13, height=14,
                    pal_gsea = FALSE,
                    file.name = "Heatmap_top40_X4clusters_vs_B.jpeg",
                    title = "Top 40 DE genes in 4 MCL clusters and B cells",
                    save.path = paste0(path,choose))
}
#==== Figure 3-D ===========
choose <- c("X4clusters","X4cluster_vs_Normal")[2]
res = read.csv(file = paste0("Yang/PALIBR Figures legends methods/Figure 3/Figure Sources/",
                             choose,"/",choose,"_41-FC0.csv"),
               row.names = 1, stringsAsFactors=F)
table(res$cluster)
head(res)
res = res[order(res["p_val_adj"]),]
head(res, 20)
(clusters <- unique(res$cluster))
hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v6.2.symbols.gmt")
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(hallmark) = gsub("\\_"," ",names(hallmark))
hallmark$`NF-kB signaling` =  read.delim("data/200222 NFKB pathway gene list.txt") %>%
        pull %>% as.character()
hallmark$`MYC TARGETS` = c(hallmark$`MYC TARGETS V1`,hallmark$`MYC TARGETS V2`)
select = c("TNFA SIGNALING VIA NFKB","MYC TARGETS","E2F TARGETS","OXIDATIVE PHOSPHORYLATION",
           "G2M CHECKPOINT","DNA REPAIR","REACTIVE OXIGEN SPECIES PATHWAY",
           "MTORC1 SIGNALING","GLYCOLYSIS","UNFOLDED PROTEIN RESPONSE","FATTY ACID METABOLISM",
           "CHOLESTEROL HOMEOSTASIS","HYPOXIA", "PI3K AKT MTOR SIGNALING","INTERFERON ALPHA RESPONSE",
           "TGF BETA SIGNALING","APOPTOSIS","IL6 JAK STAT3 SIGNALING", "INTERFERON GAMMA RESPONSE",
           "P53 PATHWAY","IL2 STAT5 SIGNALING", "INFLAMMATORY RESPONSE","NF-kB signaling","KRAS SIGNALING DN",
           "KRAS SIGNALING UP")
hallmark =  hallmark[select]
Names = c("TNFa signaling via NF-kB","Myc targets","E2F targets",
          "Oxidative phosphorylation","G2M checkpoint","DNA repair",
          "Reactive oxygen species pathway","mTORC1 signaling","Glycolysis",
          "Unfolded protein response","Fatty acid metabolism","Cholesterol homeostasis",
          "Hypoxia","PI3K/AKT/mTOR signaling","Interferon alpha response","TFG b signaling",
          "Apoptosis", "IL-6/JAK/STAT3 signaling", "Interferon gamma response",
          "p53 pathway","IL-2/STAT5 signaling","Inflammatory response",
          "NF-kB signaling","KRAS signaling dn","KRAS signaling up")
names(hallmark) = Names
# Now, run the fgsea algorithm with 1000 permutations:
fgseaRes = FgseaDotPlot(stats=res, pathways=hallmark,
                        padj = 0.25,pval = 0.05,
                        order.yaxis.by = c("C4","NES"),
                        order.xaxis = if(choose == "X4cluster_vs_Normal") {
                                c("B_cells", paste0("C",1:4))} else paste0("C",1:4),
                        decreasing = F,
                        Rowv = F,Colv = F,
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = "Hallmark",rotate.x.text = T,
                        title = "Cluster",
                        font.xtickslab=10, font.main=14, font.ytickslab = 10,
                        font.legend = list(size = 10),font.label = list(size = 10),
                        do.return = T, do.print = T,
                        width = 4.3,height = 4,save.path = paste0(path,choose))
file.rename(paste0(path,choose,"/Dotplot_Cluster_Hallmark_0.25_0.05.jpeg"),
            paste0(path,choose,"/Dotplot_",choose,"_Hallmark_0.25_0.05.jpeg"))
write.csv(fgseaRes, file = paste0(path,choose,"/Dotplot",choose,"_FDR0.25_pval0.05.csv"))

##########################################
#==== Figure 3C ===========
##########################################
path <- "Yang/Figure 3/Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)

Idents(object) = "groups"
table(Idents(object))
sub_object <- subset(object, idents = c("Normal", "Untreated"))
sub_object$orig.ident %<>% gsub("N02|N01|N03","Normal",.)
keep_normal = which(sub_object$orig.ident %in% "Normal")
keep_normal = sample(keep_normal,size = length(keep_normal)/3, replace = F)
keep_untreated = which(!sub_object$orig.ident %in% "Normal")
sub_object = sub_object[,c(keep_normal,keep_untreated)]
Idents(sub_object) = "orig.ident"
table(Idents(sub_object))

Idents(sub_object) = "cell.types"
TSNEPlot.1(sub_object, pt.size =0.3,
           text.size = 14,no.legend = T,
           group.by = "cell.types",split.by = "orig.ident",legend.size = 0,
           cols = ExtractMetaColor(sub_object), ncol = length(unique(sub_object$orig.ident)),
           unique.name = "groups", do.print = T,do.return = F,border = T,
           width=8.5, height=2, save.path = path)

###############################
# All pairwise heatmaps:
###############################
path <- "Yang/B_pairwise_heatmaps/"
if(!dir.exists(path)) dir.create(path, recursive = T)
group.colors = c("#181ea4","#5f66ec","#f46072","#e6001c")
# =========Doheatmap for Normal / MCL ============
B_cells_MCL = readRDS(file = "data/MCL_41_B_20200225.rds")
# remove B cell in MCL, and keep the B cells in Normal
B_cells_MCL$tests %<>% gsub("test.*","test",.)
B_cells_MCL$tests_cell.types = paste(B_cells_MCL$tests,
                                     B_cells_MCL$cell.types, sep = "_")
Idents(B_cells_MCL) = "tests_cell.types"
B_cells_MCL %<>% subset(idents = c("control_B_cells","test_MCL"))

B_cells_MCL$orig.ident %<>% gsub("N01|N02|N03","Normal",.)
B_cells_MCL$X4_orig.ident = paste(B_cells_MCL$orig.ident,
                                   B_cells_MCL$X4clusters, sep = "_")
B_cells_MCL@meta.data$X4_orig.ident = gsub('^Normal_.*', 'Normal', B_cells_MCL@meta.data$X4_orig.ident)
table(B_cells_MCL@meta.data$X4_orig.ident)
Idents(B_cells_MCL) = "X4_orig.ident"
table(Idents(B_cells_MCL))
df_samples <- readxl::read_excel("doc/191001_scRNAseq_info.xlsx",sheet = "heatmap")
list_samples <- df2list(df_samples)
scRNAseq_info <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
colnames(scRNAseq_info) <- colnames(scRNAseq_info) %>% tolower
list_samples <- lapply(list_samples[c("MCL","MCL.1","MCL.2")],
                       function(x) plyr::mapvalues(x,
                                               from = scRNAseq_info$sample,
                                               to = scRNAseq_info$`sample name`))
all(list_samples %>% unlist %>% as.vector %>% unique %in%
            B_cells_MCL$orig.ident)
Idents(B_cells_MCL) = "orig.ident"
B_cells_MCL %<>% subset(idents = "N04", invert = T)
markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
(block <- VariableFeatures(B_cells_MCL) %>% tail(2))

choose = c("X4cluster_vs_Normal","X4clusters")[1]
run_DE = F
for(sample in list_samples$MCL[1]){
        subset.MCL <- subset(B_cells_MCL, idents = c("Normal",sample))

        # SplitTSNEPlot======
        Idents(subset.MCL) = "X4_orig.ident"
        subset.MCL %<>% sortIdent()
        TSNEPlot.1(subset.MCL, split.by = "orig.ident",pt.size = 1,label = F,
                   do.print = F, unique.name = T)

        # remove cluster with less than 3 cells======
        table_subset.MCL <- table(subset.MCL$X4_orig.ident) %>% as.data.frame
        keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character()

        (X4_cluster <- keep.MCL %>% unique %>%
                        gsub('.*\\_C',"",.) %>%
                        sub("Normal","",.) %>%
                        as.numeric %>% sort )

        print(ident.1 <- rep("Normal",length(X4_cluster)))
        print(ident.2 <- paste(sample,X4_cluster,sep="_C"))
        subset.MCL <- subset(subset.MCL, idents = c(ident.1,ident.2))

        # FindAllMarkers.UMI======
        if(run_DE) {
                gde.markers <- FindPairMarkers(subset.MCL, ident.1 = c(ident.1,ident.2),
                                               ident.2 = c(ident.2,ident.1), only.pos = T,
                                               logfc.threshold = 0.1,min.cells.group =3,
                                               min.pct = 0.1,return.thresh = 0.05,
                                               latent.vars = "nCount_SCT")
                write.csv(gde.markers, paste0(path,"DE_analysis_files/",sample,"_vs_Normal.csv"))
        }

        gde.markers = read.csv(paste0(path,"DE_analysis_files/",sample,"_vs_Normal.csv"),row.names = 1)
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        Top_n = 40
        top <-  gde.markers %>% group_by(cluster1.vs.cluster2) %>%
                top_n(Top_n, avg_logFC) %>% as.data.frame()
        write.csv(top, paste0(path,"DE_analysis_files/","top40_",sample,"_vs_Normal.csv"))
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        featuresNum <- make.unique(features, sep = ".")
        exp = AverageExpression(subset.MCL[features,],
                                assays = "SCT") %>% .$SCT
        exp = MakeUniqueGenes(object = exp, features = features)
        exp[tail(VariableFeatures(object = B_cells_MCL), 2),] =0
        scale_exp <- exp %>% t %>% scale %>% t
        colnames(scale_exp)
        group.by = factor(c("Normal",ident.2), levels = c("Normal",ident.2))
        DoHeatmap.matrix(scale_exp, features = featuresNum,
                         group.by = group.by,size = 6,angle =90,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         group.colors = c("#31aa3a",group.colors),
                         width=length(group.by)+0.75/2, height=22,res=600,no.legend = T,
                         cex.row=5,group.bar.height = 0.03,
                         do.print = T,
                         file.name = paste0("Heatmap_top",Top_n,"_",sample,"_",choose,".jpeg"),
                         title = paste("40 DEGs in",sample,"_Normal"),
                         save.path = paste0(path,choose))
}

# =========Doheatmap for MCL.1 / MCL.2 ============
df_samples <- readxl::read_excel("doc/191001_scRNAseq_info.xlsx",sheet = "heatmap")
list_samples <- df2list(df_samples)
scRNAseq_info <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
colnames(scRNAseq_info) <- colnames(scRNAseq_info) %>% tolower
list_samples <- lapply(list_samples[c("MCL","MCL.1","MCL.2")],
                       function(x) plyr::mapvalues(x,
                                                   from = scRNAseq_info$sample,
                                                   to = scRNAseq_info$`sample name`))
all(list_samples %>% unlist %>% as.vector %>% unique %in%
            B_cells_MCL$orig.ident)
(block <- VariableFeatures(B_cells_MCL) %>% tail(2))
markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))

Idents(B_cells_MCL) = "orig.ident"
choose = c("X4cluster_vs_Normal","X4clusters")[1]
run_DE = F
for(i in c(1,4)){ #1:10

        (samples1 = list_samples$MCL.1[i])
        (samples2 = list_samples$MCL.2[i])

        subset.MCL <- subset(B_cells_MCL, idents = c(samples1,samples2))
        subset.MCL@meta.data$orig.ident %<>% factor(levels = c(samples1,samples2))
        # remove cluster with less than 3 cells======

        table_subset.MCL <- table(subset.MCL@meta.data$X4_orig.ident) %>% as.data.frame
        (keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 3,"Var1"] %>% as.character())
        (X4_cluster <- keep.MCL %>% unique %>%
                        gsub('.*\\_C',"",.) %>% as.numeric %>% sort %>% .[duplicated(.)])

        print(ident.1 <- paste(samples1,X4_cluster,sep="_C"))
        print(ident.2 <- paste(samples2,X4_cluster,sep="_C"))

        #---SplitTSNEPlot----
        Idents(subset.MCL) = "X4_orig.ident"
        subset.MCL <- subset(subset.MCL, idents = c(ident.1,ident.2))

        Idents(subset.MCL) %<>% factor(levels = c(ident.1,ident.2))
        TSNEPlot.1(subset.MCL, split.by = "orig.ident",pt.size = 1,label = F,
                   do.return = F,do.print = F, unique.name = T)

        if(run_DE) {
                gde.markers <- FindPairMarkers(subset.MCL, ident.1 = c(ident.1,ident.2),
                                       ident.2 = c(ident.2,ident.1), only.pos = T,
                                       logfc.threshold = 0.05,min.cells.group =3,
                                       min.pct = 0.1,
                                       latent.vars = "nCount_SCT")
                write.csv(gde.markers, paste0(path,"DE_analysis_files/",samples1,"_vs_",samples2,".csv"))
        }
        gde.markers = read.csv(paste0(path,"DE_analysis_files/",samples1,"_vs_",samples2,".csv"),row.names = 1)
        print(table(gde.markers$cluster1.vs.cluster2))
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        Top_n = 40
        top <-  gde.markers %>% group_by(cluster1.vs.cluster2) %>%
                top_n(Top_n, avg_logFC) %>% as.data.frame()
        write.csv(top, paste0(path,"DE_analysis_files/","top40_",samples1,"_vs_",samples2,".csv"))
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        featuresNum <- make.unique(features, sep = ".")
        exp = AverageExpression(subset.MCL[features,],
                                assays = "SCT") %>% .$SCT
        exp = MakeUniqueGenes(object = exp, features = features)
        exp[tail(VariableFeatures(object = B_cells_MCL), 2),] =0
        scale_exp <- exp %>% t %>% scale %>% t
        colnames(scale_exp)
        group.by = factor(c(ident.1, ident.2), levels = c(ident.1, ident.2))

        DoHeatmap.matrix(scale_exp, features = featuresNum,
                         group.by = group.by,size = 6,angle =90,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,label = T,
                         group.colors = rep(group.colors[X4_cluster], 2),
                         width=length(group.by)+0.75/2, height=22,res=600,no.legend = T,
                         cex.row=5,
                         do.print = T,
                         file.name = paste0("Heatmap_top",Top_n,"_",samples1,"_vs_",samples2,"_",choose,".jpeg"),
                         title = paste0("40 DEGs in ",samples1,"_",samples2),
                         save.path = paste0(path,choose))
}
###############################
#### Fig 4
###############################
path = "Yang/PALIBR Figures legends methods/Figure 4/"
(load(file="data/MCL_41_harmony_20200225.Rda"))
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "cell.types"
object %<>% subset(idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)
table(Idents(object))

meta.data = object@meta.data[,c("orig.ident", "cell.types")]
(mytable <- table(meta.data$orig.ident,meta.data$cell.types))
df <- prop.table(mytable, 1) %>% as.data.frame.matrix
cell_list <- list(as.data.frame.matrix(mytable),df*100)
names(cell_list) = c("cell number", "percentage % in each specimen")
write.xlsx(cell_list, file = paste0(path,"cell.types-orig.ident.xlsx"),rowNames =TRUE,
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
# subset
Idents(object) = "cell.types"
cell_types <- c("T_cells:CD4+", "T_cels:CD8+","NK_cells","Myeloid cells")
exp <- list()
for(i in seq_along(cell_types)) {
        sub_object <- subset(object, idents = cell_types[i])
        Idents(sub_object) = "orig.ident"
        exp[[i]] = AverageExpression(sub_object, assays = "SCT")
        exp[[i]] = exp[[i]]$SCT
}
names(exp) = gsub(":","_",cell_types)
write.xlsx(exp, file = paste0(path,"UMI-cell.types-orig.ident.xlsx"),rowNames =TRUE,
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
#### Heatmap T ===========
# select pair
opts = data.frame(cell.types = rep(c("T_cells:CD8+","T_cells:CD4+"),  each = 23),
                  ident.1 = rep(c("Pt17_7","Pt17_12","Pt17_31","Pt28_4","Pt28_28",
                                  "Pt25_1_8","Pt25_24","Pt25_25Pd","Pt25_25Pd","Pt11_28",
                                  "Pt25_1","PtU01","PtU02","PtU03","PtU04",
                                  "Pt17_LN1","Pt17_2","Pt17_7","Pt17_12","Pt17_31",
                                  "PtB13_Ibp","PtB13_Ib1","PtB13_IbR"),2),
                  ident.2 = rep(c("Pt17_2","Pt17_2","Pt17_2","Pt28_1","Pt28_1",
                                  "Pt25_1","Pt25_1","Pt25_1","Pt25_24","Pt11_14",
                                  rep("N01",13)),2),
                  stringsAsFactors = F)
(samples <- sort(unique(c(opts$ident.1,opts$ident.2))))
samples = c("N01","PtU01","PtU02","PtU03","PtU04",
            "Pt17_LN1","Pt17_2","Pt17_7","Pt17_31",
            "Pt25_SB1","Pt25_1","Pt25_1_8","Pt25_24","Pt25_25Pd","Pt25_AMB25Pd",
            "Pt28_LN1","Pt28_1","Pt28_4","Pt28_28",
            "Pt11_LN1", "Pt11_1","Pt11_14","Pt11_28",
            "PtB13_Ibp","PtB13_Ib1","PtB13_IbR")
Idents(object) = "cell.types"

T_cells <- c("CD8+","CD4+")
for(i in seq_along(T_cells)){
        df <- readxl::read_excel(paste0(path, "200514_DE genes CD8_CD4 T cell.xlsx"),sheet = T_cells[i])
        features <- unique(df$Up, df$Down) %>% na.omit

        sub_object <- subset(object, idents = paste0("T_cells:",T_cells[i]))
        Idents(sub_object) = "orig.ident"
        sub_object %<>% subset(idents = samples[samples %in% Idents(sub_object)])
        sub_object %<>% ScaleData(features = features)
        exp = AverageExpression(sub_object[features,],
                                      assays = "SCT") %>% .$SCT
        #if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

        samples = samples[samples %in% colnames(exp)]
        jpeg(paste0(path,"heatmap2_",T_cells[i],"_short.jpeg"), units="in", width=10, height=7,res=600)
        heatmap.2(as.matrix(exp[,samples]),
                  breaks = seq(-3,3,length.out = 300),
                  cexRow = 0.5,
                  #col = sns.RdBu_r,
                  dendrogram = "both",
                  margins = c(5,5),
                  col = bluered(299),
                  key.xlab = "scale log nUMI",
                  Colv = F,
                  Rowv = T,
                  scale= "row",
                  trace = "none",
                  density.info="none",
                  main = paste(T_cells[i],"T cells DE genes in MCL and normal blood"))
        dev.off()
}

###############################
#### Fig 4E #####
###############################
path <- "Yang/Figure 4/"
if(!dir.exists(path)) dir.create(path, recursive = T)
B_cells_MCL = readRDS(file = "data/MCL_41_B_20200225.rds")
B_cells_MCL$X4_orig.ident = paste(B_cells_MCL$orig.ident,
                                  B_cells_MCL$X4clusters, sep = "_")
Idents(B_cells_MCL) = "cell.types"
B_cells_MCL %<>% subset(idents = "MCL")
(samples1 = list_samples$MCL.1[1])
(samples2 = list_samples$MCL.2[1])
Idents(B_cells_MCL) = "orig.ident"
subset.MCL <- subset(B_cells_MCL, idents = c(samples1,samples2))
subset.MCL$orig.ident %<>% factor(levels = c(samples1,samples2))
# merge C1+C2 and C3+C4

subset.MCL$X4_orig.ident %<>% gsub("C1|C2","C1+2",.)
subset.MCL$X4_orig.ident %<>% gsub("C3|C4","C3+4",.)

print(ident.1 <- paste0(samples1,c("_C1+2","_C3+4")))
print(ident.2 <- paste0(samples2,c("_C1+2","_C3+4")))

Idents(subset.MCL) = "X4_orig.ident"
TSNEPlot.1(subset.MCL, split.by = "orig.ident",pt.size = 1,label = F,
           do.return = T,do.print = F, unique.name = T)

gde.markers <- FindPairMarkers(subset.MCL, ident.1 = ident.2,
                               ident.2 = ident.1, only.pos = F,
                               logfc.threshold = 0.1,min.cells.group =3,
                               min.pct = 0.1,return.thresh = 0.05,
                               latent.vars = "nCount_SCT")
write.csv(gde.markers, paste0(path,samples2,"_vs_",samples1,".csv"))
gde.markers = read.csv(paste0(path,samples2,"_vs_",samples1,".csv"),row.names = 1)
print(table(gde.markers$cluster1.vs.cluster2))
(mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
GC()
gde.markers$cluster1.vs.cluster2 %<>% as.character()
Clusters <- c("C1","C3")
for(i in seq_along(Clusters)){
        gde <- gde.markers[grepl(Clusters[i], gde.markers$cluster1.vs.cluster2),]
        p <- VolcanoPlots(data = gde, cut_off_pvalue = 0.0000001, cut_off_logFC = 0.25,
                          top = 20, cols = c("#0000ff","#d2dae2","#ff0000"),alpha=0.8, size=2,
                          legend.size = 12)+
                ggtitle(paste0(Clusters[i],"/",i*2," AMB/SB"))+
                theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))
        jpeg(paste0(path,"VolcanoPlots_",Clusters[i],"_",i*2,"_AMB_SB.jpeg"), units="in", width=10, height=7,res=600)
        print(p)
        dev.off()
}



### Fig. 3F ==========
path <- "Yang/PALIBR Figures legends methods/Figure 7/Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)
group.colors = c("#181ea4","#5f66ec","#f46072","#e6001c")

B_cells_MCL = readRDS(file = "data/MCL_41_B_20200225.rds")
choose = c("X4clusters","X4cluster_vs_Normal")[1]
if(choose == "X4clusters"){
        Idents(B_cells_MCL) = "cell.types"
        B_cells_MCL %<>% subset(idents = "MCL")
        Idents(B_cells_MCL) = "orig.ident"
        B_cells_MCL %<>% subset(idents = "Pt10_LN2Pd")
        Idents(B_cells_MCL) = "X4clusters"
        B_cells_MCL %<>% sortIdent()
        table(Idents(B_cells_MCL))
        system.time(MCL_markers <- FindAllMarkers.UMI(B_cells_MCL,
                                                      only.pos = T,
                                                      test.use = "MAST",
                                                      logfc.threshold = 0.05,
                                                      min.pct = 0.1,return.thresh = 0.05,
                                                      latent.vars = "nCount_SCT"))
        write.csv(MCL_markers,paste0(path,"Pt10_LN2Pd_X4clusters_FC0.05_markers.csv"))
        X4clusters_markers = read.csv(file= paste0(path,"Pt10_LN2Pd_X4clusters_FC0.05_markers.csv"),
                                      row.names = 1, stringsAsFactors=F)
        table(X4clusters_markers$cluster)
        X4clusters_markers$cluster %<>% factor(levels = paste0("C",1:4))
        markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
        (MT_gene <- grep("^MT-",X4clusters_markers$gene))
        X4clusters_markers = X4clusters_markers[-MT_gene,]
        Top_n = 40
        top = X4clusters_markers %>% group_by(cluster) %>%
                top_n(Top_n, cluster) %>% top_n(Top_n, avg_logFC)
        unique(top$cluster)
        top = top[order(top$cluster),]
        write.csv(top,paste0(path,"top40_genes_heatmap.csv"))
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        #DoHeatmap.1======
        featuresNum <- make.unique(features, sep = ".")
        exp = AverageExpression(B_cells_MCL[features,],
                                assays = "SCT") %>% .$SCT
        exp = MakeUniqueGenes(object = exp, features = features)
        exp[tail(VariableFeatures(object = B_cells_MCL), 2),] =0
        scale_exp <- exp %>% t %>% scale %>% t
        colnames(scale_exp)
        (group.by = unique(top$cluster))
        DoHeatmap.matrix(scale_exp, features = featuresNum,
                         group.by = 1:4,size = 6,angle = 0,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         width=1.5, height=10,res=600,no.legend = T,
                         cex.row=5,
                         group.colors = group.colors,
                         do.print = T,
                         file.name = paste0("Heatmap_top",Top_n,"_",sample,"_",choose,".jpeg"),
                         title = "40 DEGs",
                         save.path = path)
        file.rename(paste0(path,"Heatmap_top40_object_MCL_B_cells_X4clusters_Legend.jpeg"),
                    paste0(path,"Heatmap_top40_Pt10_LN2Pd_X4clusters.jpeg"))
}

#==== Figure ===========
path <- "Yang/Figure Sources/tSNE plots/Groups/"
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data
(load(file="data/MCL_41_harmony_20200225.Rda"))
# reduce normal sample
object$orig.ident %<>% gsub("N01|N02|N03","Normal",.)
all_normal = which(object$orig.ident %in% "Normal")
rm_normal = sample(all_normal,size = length(all_normal)/3*2, replace = F)
object = object[,-rm_normal]

# order the sample
df_samples <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
table(object$orig.ident)

(df_samples = df_samples[df_samples$`sample name` %in% object$orig.ident,])
(samples = c("Normal",df_samples$`sample name`[df_samples$`sample name` %in% object$orig.ident]))
object$orig.ident %<>% factor(levels = samples)

# preprocess
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "cell.types"
object %<>% subset(idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)
table(Idents(object))

object$cell.types <- plyr::mapvalues(object@meta.data$cell.types,
                                     from = c("B_cells","MCL",
                                              "Myeloid cells",
                                              "NK_cells","T_cells:CD4+",
                                              "T_cells:CD8+"),
                                     to = c("B","MCL",
                                            "Monocytes",
                                            "NK","CD4 T",
                                            "CD8 T"))
object$cell.types %<>% as.character()
object$cell.types.colors = object$cell.types.colors


Idents(object) = "groups"
(groups = Idents(object) %>% unique %>% as.character %>% sort)
for(i in 1:length(groups)){
        sub_object <- subset(object, idents = groups[i])
        Idents(sub_object) = "cell.types"
        TSNEPlot.1(object = sub_object, label = F, label.repel = F, group.by = "cell.types",
                   cols = ExtractMetaColor(sub_object),no.legend = F,border = T,
                   pt.size = 1, do.print = T,do.return = F,legend.size = 25,
                   unique.name = "groups",
                   title.size = 20,title = paste("tSNE plot of",groups[i],"samples"),
                   units= "in",width=9, height=7,hjust =0.5, save.path = path)
        Progress(i-2,length(groups)-2)
}
path <- "Yang/Figure Sources/tSNE plots/Groups_split/"
if(!dir.exists(path)) dir.create(path, recursive = T)

Idents(object) = "groups"
for(i in 1:length(groups)){
        sub_object <- subset(object, idents = c("Normal",groups[i]))
        Idents(sub_object) = "cell.types"
        TSNEPlot.1(object = sub_object, label = F, label.repel = F,
                   group.by = "cell.types",split.by = "orig.ident",
                   cols = ExtractMetaColor(sub_object),no.legend = T,border = T,
                   pt.size = 0.2, do.print = T,do.return = F,legend.size = 12,
                   unique.name = "groups",
                   title.size = 14,title = paste("tSNE plot of Normal and",groups[i],"samples"),
                   units= "in",
                   width=length(unique(sub_object$orig.ident))*2.5,
                   height = 3,hjust =0.5, save.path = path)
        Progress(i-2,length(groups)-2)
}

path <- "Yang/Figure Sources/tSNE plots/B_Groups_split/"
if(!dir.exists(path)) dir.create(path, recursive = T)
Idents(object) = "cell.types"
object <- subset(object, idents = c("B","MCL"))
Idents(object) = "groups"

for(i in 1:length(groups)){
        sub_object <- subset(object, idents = c("Normal",groups[i]))
        Idents(sub_object) = "cell.types"
        TSNEPlot.1(object = sub_object, label = F, label.repel = F,
                   group.by = "cell.types",split.by = "orig.ident",
                   cols = ExtractMetaColor(sub_object),no.legend = T,border = T,
                   pt.size = 0.2, do.print = T,do.return = F,legend.size = 12,
                   unique.name = "groups",
                   title.size = 14,title = paste("tSNE plot of Normal and",groups[i],"samples"),
                   units= "in",
                   width=length(unique(sub_object$orig.ident))*2.5,
                   height = 3,hjust =0.5, save.path = path)
        Progress(i-2,length(groups)-2)
}

###############################
#### Fig 4
###############################
path = "Yang/PALIBR/Fig. 6/"
(load(file="data/MCL_41_harmony_20200225.Rda"))
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
object@meta.data[object$SCT_snn_res.0.8 %in% c(5,16), "cell.types"] = "Monocytes"
Idents(object) = "cell.types"
object %<>% subset(idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)
table(Idents(object))

samples = c("N01","PtU01","PtU02","PtU04",#"PtU03",
            "Pt17_2","Pt17_7","Pt17_31",#"Pt17_LN1",
            "Pt25_1","Pt25_24","Pt25_25Pd"#,"Pt25_SB1","Pt25_1_8","Pt25_AMB25Pd",
            #"Pt28_LN1","Pt28_1","Pt28_4","Pt28_28",
            #"Pt11_LN1", "Pt11_1","Pt11_14","Pt11_28",
            #"PtB13_Ibp","PtB13_Ib1","PtB13_IbR"
)
samples %<>% factor(levels = as.character(samples))

Idents(object) = "cell.types"
gene_list <- list("T_cells_CD4+" = c("IL32", "CD52", "TXNIP", "EEF1G", "EML4","TSC22D3",
                                     "KLRG1","TOX","PDCD1","PCNA"),#"TNFAIP3",
                  "T_cells_CD8+" = c("IL32", "CD52", "LTB", "IFITM1", "JUN", "JUNB", "RPS26",
                                     "EEF1G", "EML4", "TSC22D3","KLRG1","TOX","PDCD1","PCNA"))#"PIK3R1", "TNFAIP3",
cell.types <- c("T_cells_CD4+","T_cells_CD8+")

#' adjust hex color scheme acorrding to the range, put white color at middle
#' @example skewed_corlor(sns.RdBu_r, Min = -1, Max = 2)
skewed_corlor <- function(hex_color, Min, Max){
        #if(Min > 0) return(hex_color[7:12])
        #if(Max < 0) return(hex_color[1:6])
        if(Min < 0 & Max > 0) {
                l = length(hex_color)
                Range = max(Max, -Min)*2
                remove <- (Min + Max) / Range * l
                if(remove > 0) return(hex_color[ceiling(remove):l])
                if(remove < 0) return(hex_color[1:(l+min(-1,remove)+1)])
        }
}
#sns.color_palette("RdBu_r", 15).as_hex()

# average scale before
for(i in seq_along(cell.types)){
        #sns.RdBu <- if(i ==1) sns.RdBu_r_199 else sns.RdBu_r_15
        features = gene_list[[cell.types[i]]]
        sub_object <- subset(object, idents = gsub("_CD",":CD",cell.types[i]))
        Idents(sub_object) = "orig.ident"
        sub_object %<>% subset(idents = samples[samples %in% Idents(sub_object)])
        sub_object %<>% ScaleData(features = features)
        exp = AverageExpression(sub_object[features,],
                                assays = "SCT") %>% .$SCT
        #if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        exp %<>% t %>% scale %>% t
        S = samples[samples %in% colnames(exp)]
        DoHeatmap.matrix(exp[,as.character(S)], features = features,
                         group.by = S,size = 4,angle = 45,
                         draw.lines =F, raster = FALSE,
                         colors = skewed_corlor(sns.RdBu_r_199,
                                                Min = min(exp[,as.character(S)]),
                                                Max = max(exp[,as.character(S)])),
                         width=ifelse(i == 1,4,4), height=ifelse(i == 1,3,4),res=600,
                         no.legend = T,
                         cex.row=10,
                         group.bar.height = 0,
                         do.print = T,
                         file.name = paste0("Heatmap_",cell.types[i],"_","average_scale.jpeg"),
                         save.path = paste0(path,"Heatmaps"))
        print(min(exp[,as.character(S)]))
        print(max(exp[,as.character(S)]))
}

