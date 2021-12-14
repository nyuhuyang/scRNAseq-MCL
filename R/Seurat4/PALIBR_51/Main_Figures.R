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
library(stringr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")

# load data
object = readRDS(file = "data/MCL_SCT_51_20210724.rds")
df_samples <- readxl::read_excel("doc/20210715_scRNAseq_info.xlsx",sheet = "fastq")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples %<>% filter(sequence %in% "GEX") %>% filter(phase %in% "PALIBR_I") %>%
        filter(sample != "Pt11_31")
table(object$orig.ident %in% df_samples$sample)
table(df_samples$sample %in% object$orig.ident)

(df_samples = df_samples[df_samples$sample %in% object$orig.ident,])
(samples = df_samples$sample[df_samples$sample %in% object$orig.ident])

# preprocess
object$orig.ident %<>% factor(levels = samples)
Idents(object) = "orig.ident"
object %<>% subset(idents = "Pt2_30Pd", invert = T)
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "cell.types"
object %<>% subset(idents = c("HSC/progenitors","Nonhematopoietic cells","unknown"), invert = TRUE)
table(Idents(object))

#==== Figure 3-A ===========
path <- "Yang/PALIBR/51_samples/Fig. 3/"
if(!dir.exists(path)) dir.create(path, recursive = T)

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
read_path = "Yang/Figure Sources/51_samples/heatmaps_full/"
path <- "Yang/PALIBR/51_samples/Fig. 3/"
if(!dir.exists(path)) dir.create(path, recursive = T)

B_cells_MCL = readRDS(file = "data/MCL_SCT_51_20210724.rds")
DefaultAssay(B_cells_MCL) = "SCT"

B_cells_MCL = subset(B_cells_MCL, subset =  Doublets == "Singlet"
                & X4cluster  %in% c("1","2","3","4")
)
B_cells_MCL = subset(B_cells_MCL, subset =  orig.ident %in% c("Pt3_BMA72_6","Pt3_72_6"), invert = T)
Idents(B_cells_MCL) = "X4cluster"
B_cells_MCL$X4cluster %<>% factor(levels=c("1","2","3","4"))

markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
group.colors = c('#40A635','#FE8205','#8861AC','#E83C2D') #(green, orange, light purple,red)
choose = c("X4cluster","X4cluster_vs_Normal","X4cluster_vs_B")[1]
filter_by = c("top40","logFC0.01")[1]
B_cells_MCL[["patient"]] = gsub("_.*","",B_cells_MCL$orig.ident)
B_cells_MCL[["samples"]] = droplevels(B_cells_MCL[["orig.ident"]])
B_cells_MCL$PB_or_not = plyr::mapvalues(B_cells_MCL$tissue,
                                             from = c("AMB","BM","BMA","LN","SB"),
                                             to = rep("non-PB",5))
B_cells_MCL$PB_or_not %<>% factor(levels = c("PB","non-PB"))
Idents(B_cells_MCL) = "PB_or_not"

B_MCL_list <- pbapply::pblapply(c("PB","non-PB"),function(x) subset(B_cells_MCL,idents = x))

B_cells_MCL <- Reduce(function(x, y) merge(x, y, do.normalize = F), B_MCL_list)

choose = "X4cluster"
#B_cells_MCL_copy <- B_cells_MCL
B_cells_MCL <- B_cells_MCL_copy
#B_cells_MCL %<>% subset(idents = "Pt25")
#cells <- sample(colnames(B_cells_MCL), size = 1000)
#B_cells_MCL %<>% subset(cells = cells)
B_cells_MCL@meta.data$samples %<>% factor(levels = df_samples$sample)
if(choose == "X4cluster"){
        Idents(B_cells_MCL) = "orig.ident"

        B_cells_MCL %<>% subset(idents = c("N01","N02","N03","N04"), invert = T)
        B_cells_MCL@meta.data$samples %<>% droplevels()
        Idents(B_cells_MCL) = "X4cluster"

        Idents(B_cells_MCL) %<>% factor(levels = 1:4)
        table(Idents(B_cells_MCL))
        X4clusters_markers = read.csv(file= paste0(read_path,choose,"/MCL_only_",choose,"_51-FC0.csv"),
                                      row.names = 1, stringsAsFactors=F)
        table(X4clusters_markers$cluster)
        X4clusters_markers$cluster %<>% factor(levels = 1:4)
        markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
        # remove mito
        (MT_gene <- grep("^MT-",X4clusters_markers$gene))
        if(length(MT_gene)>0) X4clusters_markers = X4clusters_markers[-MT_gene,]
        # remove IGL, IGK, IGK
        (IG.gene <- grep(pattern = "^IGL|^IGK|^IGH", x = X4clusters_markers$gene))
        if(length(IG.gene)>0) X4clusters_markers = X4clusters_markers[-IG.gene,]
        # remove RPL, RPS
        (RP.gene <- grep(pattern = "^RPL|^RPS", x = X4clusters_markers$gene))
        if(length(RP.gene)>0) X4clusters_markers = X4clusters_markers[-RP.gene,]

        Top_n = 40
        top = switch (filter_by,
                      "top40"= X4clusters_markers %>% group_by(cluster) %>%
                              top_n(Top_n, cluster) %>% top_n(Top_n, avg_log2FC),
                      "logFC0.01" = X4clusters_markers %>% group_by(cluster) %>%
                              filter(avg_log2FC < 0.05 & avg_logFC > 0.01)
        )
        table(top$cluster)
        top = top[order(top$cluster),]
        write.csv(top,paste0(path,choose,"_",filter_by,"_genes_heatmap_.csv"))
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
                    group.by = c("X4cluster","samples"),group.bar = T,
                    group1.colors = group.colors,
                    group2.colors= Singler.colors,
                    title.size = 16, no.legend = F,size=5,hjust = 0.5,
                    group.bar.height = 0.02, label=F, cex.row= 5,
                    width=14, height=12,
                    pal_gsea = FALSE,
                    file.name = paste0("Fig3D_Heatmap_",filter_by,"_",Top_n,"_X4clusters_samples.jpeg"),
                    title =  paste("Top",Top_n, "DE genes in 4 B/MCL clusters"),
                    save.path = path,
                    nrow = 5, ncol = 8, design = c(patchwork::area(1, 1, 5, 5),
                                                   patchwork::area(3, 6, 3, 8)))
        DoHeatmap.2(B_cells_MCL, features = featuresNum, #Top_n = Top_n,
                    do.print=T, do.return=F,
                    angle = 0,
                    group.by = c("X4cluster","tissue"),group.bar = T,
                    group1.colors = group.colors,
                    group2.colors= c("#ee7576","#e94749","#E78AC3","#ff0000","#33A02C","#d80172"),
                    title.size = 16, no.legend = F,size=5,hjust = 0.5,
                    group.bar.height = 0.02, label=F, cex.row= 5,
                    width=12, height=12,
                    pal_gsea = FALSE,
                    file.name = paste0("Fig3D_Heatmap_",filter_by,"_",Top_n,"_X4clusters_tissue.jpeg"),
                    title =  paste("Top",Top_n, "DE genes in 4 B/MCL clusters"),
                    save.path = path,
                    nrow = 5, ncol = 8, design = c(patchwork::area(1, 1, 5, 5),
                                                   patchwork::area(3, 6, 3, 7)))
}


#==== Figure 3-E ===========
path <- "Yang/PALIBR/51_samples/Fig. 3/"
choose <- c("X4cluster","X4cluster_vs_Normal")[1]
res = read.csv(file = paste0(path,"/MCL_only_",choose,"_51-FC0.csv"),
               row.names = 1, stringsAsFactors=F)
table(res$cluster)
head(res, 20)
(clusters <- unique(res$cluster))
hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v7.4.symbols.gmt")
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
                        order.yaxis.by = c("1","NES"),
                        order.xaxis = if(choose == "X4cluster_vs_Normal") {
                                c("B_cells", 1:4)} else paste0(1:4),
                        decreasing = F,
                        Rowv = F,Colv = F,
                        size = " -log10(padj)", fill = "NES",
                        title = "Hallmark in MCL 4 cluster",
                        file.name = "Fig3E_fgsea",
                        font.xtickslab=10, font.main=14, font.ytickslab = 10,
                        font.legend = list(size = 10),font.label = list(size = 10),
                        do.return = T, do.print = T,
                        width = 4.3,height = 4,save.path = path)
write.csv(fgseaRes, file = paste0(path,"/Fig3E_Dotplot",choose,"_FDR0.25_pval0.05.csv"))

#==== Figure 3-F venn diagram ===========
path <- "Yang/PALIBR/51_samples/Fig. 3/"
read_path = "output/20210812_MCL_vs_N01/"
PD_sampples = c("Pt10_LN2","Pt13_BMA1","Pt25_SB1","Pt25_AMB25",
         "Pt15_BMA1","Pt18_5_1","Pt18_5_8","Pt25_25")
deg_list <- lapply(PD_sampples, function(x) {
        csv_name = paste0(read_path,"MCL_B_",x,"-vs-N01.csv")
        tmp = read.csv(csv_name,stringsAsFactors = F)
        tmp$sample = x
        tmp
})

gde.all <- bind_rows(deg_list)
gde.all %<>% filter(avg_log2FC >0) %>% filter(cluster != "N01")

for(value in c(1,0.05,0.01,0.001)){
        pos_genes <- eulerr(gde.all,group.by = "cluster", shape =  "circle",#key = c("C1","C2","C3","C4","B_cells"),
                            cut_off = "p_val_adj", cut_off_value = value,do.print = T,return.raw = T,
                            save.path = path, file.name =paste0("Fig3F_Venn_padj_",value,".jpeg"))
        euler_df <- eulerr::euler(pos_genes,shape = "circle")
        pos_genes_list <- as.list(euler_df$original.values)
        names(pos_genes_list) %<>% paste(":",pos_genes_list)
        id <- eulerr:::bit_indexr(8)
        for (i in nrow(id):1) {
                pos_genes_list[[i]] = Reduce(intersect, pos_genes[id[i,]])  %>%
                        setdiff(Reduce(union, pos_genes[!id[i,]]))
        }
        pos_genes_df <- list2df(pos_genes_list)
        pos_genes_df = pos_genes_df[,sapply(pos_genes_list,length) != 0]
        write.xlsx(pos_genes_df, asTable = F,
                   file = paste0(path,"Fig3F_Venn_postive_shared_gene_list_padj",value,".xlsx"),
                   borders = "surrounding")
        print(value)
}

gde.all <- bind_rows(deg_list)
gde.all %<>% filter(p_val_adj <0.01) %>% filter(cluster != "N01")
for(value in c(0.5,0.25,0.1, 0)){
        pos_genes <- eulerr(gde.all,group.by = "cluster", shape =  "circle",#key = c("C1","C2","C3","C4","B_cells"),
                              cut_off = "avg_log2FC", cut_off_value = value,do.print = T,return.raw = T,
                              save.path = path, file.name =paste0("Fig3F_Venn_log2FC_",value,".jpeg"))
        euler_df <- eulerr::euler(pos_genes,shape = "circle")
        pos_genes_list <- as.list(euler_df$original.values)
        names(pos_genes_list) %<>% paste(":",pos_genes_list)
        id <- eulerr:::bit_indexr(8)
        for (i in nrow(id):1) {
                pos_genes_list[[i]] = Reduce(intersect, pos_genes[id[i,]])  %>%
                        setdiff(Reduce(union, pos_genes[!id[i,]]))
        }
        pos_genes_df <- list2df(pos_genes_list)
        pos_genes_df = pos_genes_df[,sapply(pos_genes_list,length) != 0]
        write.xlsx(pos_genes_df, asTable = F,
                   file = paste0(path,"Fig3F_Venn_postive_shared_gene_list_FC",value,".xlsx"),
                   borders = "surrounding")
        print(value)
}

#=====3G heatmap Pt10 ======================
path <- "Yang/PALIBR/51_samples/Fig. 3/"
if(!dir.exists(path)) dir.create(path, recursive = T)

B_cells_MCL = readRDS(file = "data/MCL_SCT_51_20210724.rds")
DefaultAssay(B_cells_MCL) = "SCT"

B_cells_MCL = subset(B_cells_MCL, subset =  Doublets == "Singlet"
                     & X4cluster  %in% c("1","2","3","4")
)
sub_object <- subset(B_cells_MCL, subset = orig.ident == "Pt10_LN2")
Idents(sub_object) = "X4cluster"
MCL_markers <- FindAllMarkers_UMI(sub_object,
                               logfc.threshold = 0,
                               return.thresh = 1,
                               only.pos = T,
                               test.use = "MAST",
                               latent.vars = "nFeature_SCT")
write.csv(MCL_markers,paste0(path,"Pt10_LN2_X4Cluster_FC0_markers.csv"))
table(MCL_markers$cluster)
markers <- FilterGenes(sub_object,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
(MT_gene <- grep("^MT-",MCL_markers$gene))
if(length(MT_gene) >0 ) MCL_markers = MCL_markers[-MT_gene,]
Top_n = 40

top = MCL_markers %>% group_by(cluster) %>%
        top_n(40, avg_log2FC)
#top = top[order(top$cluster),]
write.csv(top,paste0(path,"Top40_Pt10_LN2_X4Cluster_FC0_markers.csv"))
features = c(as.character(top$gene),
             tail(VariableFeatures(object = sub_object), 2),
             markers)
#DoHeatmap.1======
# heatmap
featuresNum <- make.unique(features, sep = ".")
exp = AverageExpression(sub_object,features = features,
                        assays = "SCT") %>% .$SCT
exp %<>% MakeUniqueGenes(features = features)
exp %<>% t %>% scale %>% t
exp[tail(VariableFeatures(object = sub_object), 2),] =0

(group.by = as.character(1:4))
DoHeatmap.matrix(exp, features = featuresNum,
                 group.by = group.by,
                 size = 6,angle = 0,label =F,
                 draw.lines =F, raster = FALSE,
                 pal_gsea = FALSE,
                 width=2.5, height=11,res=600,no.legend = F,
                 cex.row=5,
                 group.colors = c('#40A635','#FE8205','#8861AC','#E83C2D'),
                 do.print = T,
                 save.path = path,
                 file.name = paste0("Heatmap_top40_Pt10_LN2_X4Cluster.jpeg")
)

#  remove RP
top = MCL_markers %>% filter(!str_detect(gene, "^MT-")) %>%
        filter(!str_detect(gene, "^RPL")) %>%
        filter(!str_detect(gene, "^RPS")) %>%
        group_by(cluster) %>%
        top_n(40, avg_log2FC)
unique(top$cluster)
#top = top[order(top$cluster),]
write.csv(top,paste0(path,"Top40_Pt10_LN2_X4Cluster_FC0_markers_rmPR.csv"))
features = c(as.character(top$gene),
             tail(VariableFeatures(object = sub_object), 2),
             markers)
featuresNum <- make.unique(features, sep = ".")
exp = AverageExpression(sub_object,features = features,
                        assays = "SCT") %>% .$SCT
exp %<>% MakeUniqueGenes(features = features)
exp %<>% t %>% scale %>% t
exp[tail(VariableFeatures(object = sub_object), 2),] =0

(group.by = as.character(1:4))
DoHeatmap.matrix(exp, features = featuresNum,
                 group.by = group.by,
                 size = 6,angle = 0,label =F,
                 draw.lines =F, raster = FALSE,
                 pal_gsea = FALSE,
                 width=2.5, height=10,res=600,no.legend = F,
                 cex.row=5,
                 group.colors = c('#40A635','#FE8205','#8861AC','#E83C2D'),
                 do.print = T,
                 save.path = path,
                 file.name = paste0("Heatmap_top40_Pt10_LN2_X4Cluster_rmPR.jpeg")
)

