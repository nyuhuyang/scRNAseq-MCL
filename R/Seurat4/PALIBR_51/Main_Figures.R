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

meta.data = object@meta.data
meta.data$patient = plyr::mapvalues(meta.data$orig.ident,
                                    from = df_samples$sample,
                                    to = df_samples$patient)
meta.data$response =  plyr::mapvalues(meta.data$orig.ident,
                                      from = df_samples$sample,
                                      to = df_samples$response)
meta.data$treatment =  plyr::mapvalues(meta.data$orig.ident,
                                       from = df_samples$sample,
                                       to = df_samples$treatment)
meta.data[meta.data$Doublets %in% c("Doublet-High Confidence","Doublet-Low Confidence"),"discard"] = TRUE
meta.data[meta.data$cell.types %in% c("Erythrocytes","unknown","Nonhematopoietic cells"),"discard"] = TRUE
meta.data[meta.data$orig.ident %in% "Pt2_30" & meta.data$cell.types %in% c("B_cells","MCL"),"discard"] = TRUE
meta.data$Mean.Reads.per.Cell %<>% gsub(",","",.) %>% as.integer()
meta.data$Number.of.Reads %<>% gsub(",","",.) %>% as.integer()
meta.data$Sequencing.Saturation %<>% gsub("%","",.) %>% as.numeric()

object$response %<>% factor(levels = c("Normal","Untreated","CR","PR","PD"))
object$treatment %<>% factor(levels = c("Normal","Untreated","PALIBR+Ibrutinib"))


colnames(object@meta.data)[grep("Frequency",colnames(object@meta.data))]="tcr.frequency"

saveRDS(object@meta.data, file = "output/MCL_51_20210724_meta.data_v1.rds")

table(colnames(object) == rownames(meta.data))


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
path <- "Yang/PALIBR/51_samples/Fig. 2/"
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


#==== Figure 3-D ===========
path <- "Yang/PALIBR/51_samples/Fig. 3/"
if(!dir.exists(path)) dir.create(path, recursive = T)

object = readRDS(file = "data/MCL_51_20210724.rds")
meta.data <- readRDS(file = "output/MCL_B_51_20230113_meta.data_v2.rds")
if(all(rownames(meta.data) %in% colnames(object))){
    print("all cellID within!")
    object %<>% subset(cells = rownames(meta.data))
    object@meta.data = meta.data
}

object %<>% subset(subset = Doublets == "Singlet"
                   & X4cluster  %in% c("1","2","3","4"))

#B_cells_MCL = subset(B_cells_MCL, subset =  orig.ident %in% c("Pt3_BMA72_6","Pt3_72_6"), invert = T)
Idents(object) = "X4cluster"
object$X4cluster %<>% factor(levels=c("1","2","3","4"))

markers <- FilterGenes(object,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
group1.colors = c('#40A635','#FE8205','#8861AC','#E83C2D') #(green, orange, light purple,red)
choose = "X4cluster"

excel_files <- list.files(path = "output/20230203",pattern = "X4cluster_1_2_3_4",full.names = TRUE)
degs <- readxl::read_excel(path = excel_files)
Top_n = 150
degs <- degs[grep("^MT-",degs$genes,invert = TRUE),]

top <- degs %>%
    group_by(group) %>%
    arrange(desc(avg_log2FC),.by_group = TRUE) %>%
    distinct() %>%
    top_n(Top_n, avg_log2FC) %>%
    ungroup()
table(top$group)

markers <- c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11")
markers = markers[markers %in% rownames(object)]

features = c(as.character(top$genes),
             tail(VariableFeatures(object = object), 2),
             markers)
table(duplicated(features))
sub_obj <- object[features,object$response != "Normal"]


ident <- "X4cluster"
Idents(sub_obj) <- ident
sub_obj %<>% ScaleData(features=features)
featuresNum <- make.unique(features, sep = ".")
sub_obj %<>% MakeUniqueGenes(features = features)


DoHeatmap.2(object =sub_obj, group.by = c(ident,"orig.ident"),
            features = features,
            do.print=TRUE, angle = 0, group.bar = TRUE, title.size = 20, no.legend = FALSE,size=20,hjust = 0.5,
            group.bar.height = 0.02, label=TRUE, cex.row= 1, legend.size = 12,width=14, height=12,
            group1.colors = c('#40A635','#FE8205','#8861AC','#E83C2D'),
            save.path = path,file.name = paste0("Fig3D_Heatmap_top",Top_n,"_",ident,"_orig.ident_legend.jpeg"),
            title = paste("Top",Top_n,"DE genes in 4 B/MCL cells clusters"),
            nrow = 5, ncol = 8, design = c(patchwork::area(1, 1, 5, 5),
                                           patchwork::area(3, 6, 3, 8)))

DoHeatmap.2(object =sub_obj, group.by = c(ident,"patient"),
            features = features,
            do.print=T, angle = 45, group.bar = TRUE, title.size = 20, no.legend = FALSE,size=20,hjust = 0.5,
            group.bar.height = 0.02, label=TRUE, cex.row= 1, legend.size = 12,width=14, height=13,
            group1.colors = c('#40A635','#FE8205','#8861AC','#E83C2D'),
            save.path = path,file.name = paste("Fig3D_Heatmap_top",Top_n,"_",ident,"_patient_legend.jpeg"),
            title = paste("Top",Top_n,"DE genes in 4 B/MCL cells clusters"),
            nrow = 5, ncol = 8, design = c(patchwork::area(1, 1, 5, 5),
                                           patchwork::area(3, 6, 3, 7)))

Top_n = 40
markers <- c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11")
markers = markers[markers %in% rownames(object)]
degs <- degs[grep("^MT-",degs$genes,invert = TRUE),]
degs$group %<>% factor(levels = c("1","2","3","4"))
top1 <- degs %>%
    group_by(group) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    top_n(Top_n, avg_log2FC)

table(top1$group)
openxlsx::write.xlsx(list("top40" = top1,"top150" = top), file =  paste0(path,"DEGs_X4cluster.xlsx"),
                          colNames = TRUE,rownames = T,borders = "surrounding")
features = c(as.character(top1$genes),
             tail(VariableFeatures(object = object), 2),
             markers)
table(duplicated(features))
sub_obj <- object[features,object$response != "Normal"]


ident <- "X4cluster"
Idents(sub_obj) <- ident
sub_obj %<>% ScaleData(features=features)
featuresNum <- make.unique(features, sep = ".")
sub_obj %<>% MakeUniqueGenes(features = features)
dim(sub_obj[["SCT"]]@scale.data)
table(duplicated(rownames(sub_obj[["SCT"]]@scale.data)))
table(featuresNum %in% rownames(sub_obj[["SCT"]]@scale.data))

DoHeatmap.2(object =sub_obj, group.by = c(ident,"orig.ident"),
            features = features,
            do.print=TRUE, angle = 0, group.bar = TRUE, title.size = 20, no.legend = FALSE,size=20,hjust = 0.5,
            group.bar.height = 0.02, label=TRUE, cex.row= 4.5, legend.size = 12,width=14, height=12,
            group1.colors = c('#40A635','#FE8205','#8861AC','#E83C2D'),
            save.path = path,file.name = paste0("Fig3D_Heatmap_top",Top_n,"_X4cluster_v1_orig.ident_legend.jpeg"),
            title = paste("Top",Top_n,"DE genes in 4 B/MCL cells clusters"),
            nrow = 5, ncol = 8, design = c(patchwork::area(1, 1, 5, 5),
                                           patchwork::area(3, 6, 3, 8)))

DoHeatmap.2(object =sub_obj, group.by = c(ident,"patient"),
            features = features,
            do.print=T, angle = 45, group.bar = TRUE, title.size = 20, no.legend = FALSE,size=20,hjust = 0.5,
            group.bar.height = 0.02, label=TRUE, cex.row= 4, legend.size = 12,width=14, height=13,
            group1.colors = c('#40A635','#FE8205','#8861AC','#E83C2D'),
            save.path = path,file.name = paste("Fig3D_Heatmap_top",Top_n,"_X4cluster_v1_patient_legend.jpeg"),
            title = paste("Top",Top_n,"DE genes in 4 B/MCL cells clusters"),
            nrow = 5, ncol = 8, design = c(patchwork::area(1, 1, 5, 5),
                                           patchwork::area(3, 6, 3, 7)))
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

