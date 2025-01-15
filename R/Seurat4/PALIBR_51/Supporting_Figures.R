########################################################################
#
#  07 setup environment, install libraries if necessary, load libraries
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
library(ggpubr)
library(pbapply)
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

object$orig.ident %<>% factor(levels = df_samples$sample)

object$Mean.Reads.per.Cell %>% unique %>% gsub(",","",.) %>%
        as.integer() %>% mean
object$Sequencing.Saturation %>% unique %>%
        gsub("\\.","",.) %>% gsub("%","",.) %>%
        as.integer() %>% mean

table(object$Doublets)

object %<>% subset(subset = Doublets == "Singlet")
mean(object$nCount_SCT)
mean(object$nFeature_SCT)

# Extend Data

object %<>% subset(subset = cell.types %in% c("unknown","Nonhematopoietic cells"), invert =T)
table(object$cell.types)

path <- "Yang/PALIBR/51_samples/Fig.S3/"
if(!dir.exists(path)) dir.create(path, recursive = T)

Idents(object) = "cell.types"
object %<>% sortIdent()

cell_Freq <- table(Idents(object)) %>% as.data.frame
cell_Freq$Percent <- prop.table(cell_Freq$Freq) %>% scales::percent(accuracy = 0.1)

cell_Freq = cell_Freq[order(cell_Freq$Var1),]
cell_Freq$col = ExtractMetaColor(object)
cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")
cell_Freq$Cell_Type %<>% gsub("_"," ",.)

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=7, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,
          x.text.angle = 45,
          ylab = "Cell Number",
          label = cell_Freq$Percent,
          sort.val = "desc",
          width = 1, size = 0.5,
          title = "Numbers of major cell types in total 51 samples")+NoLegend()+
        theme(plot.title = element_text(hjust = 0.5,size=15))+
        scale_y_continuous(expand = c(0, 0), limits = c(0,max(cell_Freq$Cell_Number)+2000))
dev.off()

#=======
meta.data = object@meta.data
meta.data = meta.data[!duplicated(meta.data$orig.ident),]
meta.data$Mean.Reads.per.Cell %<>% gsub(",","",.) %>% as.integer()
meta.data$Mean.Reads.per.Cell = meta.data$Mean.Reads.per.Cell/1000
meta.data$Number.of.Reads %<>% gsub(",","",.) %>% as.integer()
meta.data$Sequencing.Saturation %<>% gsub("%","",.) %>% as.numeric()

jpeg(paste0(path,"mean.Reads.per.Cell.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(meta.data, x = "Doublets",
         y= "Mean.Reads.per.Cell",
         #title = "Mean reads per cell",
         xlab = "",ylab = expression(Mean~Reads~per~Cell~x10^{"3"}),
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5,
         yscale = "log10"
         )+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,350000))+
                 theme(plot.title = element_text(hjust = 0.5,size=15),
                       axis.title.y=element_text(size =20),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank())
dev.off()

jpeg(paste0(path,"Sequencing.Saturation.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(meta.data, x = "Doublets", y= "Sequencing.Saturation",
         #title = "Sequencing saturation in each scRNA-seq",
         xlab = "",ylab = "Sequencing Saturation %",
         add = c("jitter","mean_sd"),
         #ylim = c(0, max(QC_list$cell.number)+1000),
         draw_quantiles = 0.5,
         yscale = "log10"
         )+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,6000))+
        theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.y=element_text(size =20),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()


meta.data = object@meta.data
meta.data_summary <- group_by(meta.data, orig.ident) %>%
        summarize(cell.number = n()/1000,
                  mean.nCount = mean(nCount_SCT)/1000,
                  mean.nFeature = mean(nFeature_SCT)/1000,
                  mean.percent.mt = mean(percent.mt)
                  )
meta.data_summary$Doublets = "Singlet"

jpeg(paste0(path,"S3_cell.number.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(meta.data_summary, x = "Doublets", y= "cell.number",
         #title = "Cell number in each scRNA-seq",
         xlab = "",ylab = expression(Cell~Number~x10^{"3"}),
         add = c("jitter","mean_sd"),
         #ylim = c(0, max(QC_list$cell.number)+1000),
         draw_quantiles = 0.5,
         yscale = "log10"
)+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,6000))+
        theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.y=element_text(size =20),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()


jpeg(paste0(path,"S3_mean.nCount.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(meta.data_summary, x = "Doublets", y= "mean.nCount",
         #title = "Mean transcripts per cell",
         xlab = "",ylab = expression(Mean~UMI~per~Cell~x10^{"3"}),
         add = c("jitter","mean_sd"),
         #ylim = c(0, max(QC_list$cell.number)+1000),
         draw_quantiles = 0.5,
         yscale = "log10"
)+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,6000))+
        theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.y=element_text(size =20),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()

jpeg(paste0(path,"S3_mean.nGene.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(meta.data_summary, x = "Doublets", y= "mean.nFeature",
         #title = "Mean genes per cell",
         xlab = "",ylab = expression(Mean~genes~number~per~Cell~x10^{"3"}),
         add = c("jitter","mean_sd"),
         #ylim = c(0, max(QC_list$cell.number)+1000),
         draw_quantiles = 0.5,
         yscale = "log10"
)+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,6000))+
        theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.y=element_text(size =20),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()


jpeg(paste0(path,"S3_mean.percent.mt.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(meta.data_summary, x = "Doublets", y= "mean.percent.mt",
         #title = "Mean mitochondrial gene % in each scRNA-seq",
         xlab = "",ylab = "Mean mitochondrial gene %",
         add = c("jitter","mean_sd"),
         #ylim = c(0, max(QC_list$cell.number)+1000),
         draw_quantiles = 0.5,
         yscale = "log10"
)+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,6000))+
        theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.y=element_text(size =20),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()
#==== Figure ===========
read.path <- "output/20201030/"
csv_list <- list.files(paste0(read.path,"3. CD4 T vs CD8 T"),full.names = T)
res_list = lapply(csv_list, function(fileName) {
                x = read.csv(file = fileName,row.names = 1, stringsAsFactors=F)
                x$cluster = sub("-.*","",basename(fileName))
                x$gene = rownames(x)
                x
        })
res = rbind(res_list[[1]],res_list[[2]])
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
                        padj = 0.25,pval = 0.25,
                        order.yaxis.by = c("T_cells:CD4+","NES"),
                        decreasing = F,
                        Rowv = F,Colv = F,
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = "Hallmark",rotate.x.text = T,
                        title = "Cluster",
                        font.xtickslab=10, font.main=14, font.ytickslab = 10,
                        font.legend = list(size = 10),font.label = list(size = 10),
                        do.return = T, do.print = T,
                        width = 4.3,height = 4,save.path = paste0(read.path,"/3. CD4 T vs CD8 T/"))
write.csv(fgseaRes, file = paste0(read.path,"/3. CD4 T vs CD8 T/","Dotplot_FDR1_pval1.csv"))

# support Figure 3
save_path = "Yang/PALIBR/51_samples/Fig. S2/"
if(!dir.exists(save_path)) dir.create(save_path, recursive = T)

object = readRDS(file = "data/MCL_51_20210724.rds")
meta.data <- readRDS(file = "output/MCL_B_51_20230113_meta.data_v2.rds")


if(all(rownames(meta.data) %in% colnames(object))){
    print("all cellID within!")
    object %<>% subset(cells = rownames(meta.data))
    object@meta.data = meta.data
}

object %<>% subset(subset =  Doublets == "Singlet"
                & X4cluster  %in% c("1","2","3","4"))
object$X4cluster_normal = as.character(object$X4cluster)
object@meta.data[object$orig.ident %in% c("N01","N02","N03"),
                 "X4cluster_normal"] = "Normal"
object@meta.data[(!object$X4cluster_normal %in% "Normal") & object$cell.types == "B_cells",
                 "X4cluster_normal"] = "B cells"
object %<>% subset(subset =  orig.ident == "N04" | X4cluster_normal == "B cells", invert = T)
Idents(object) = "X4cluster_normal"
object$X4cluster_normal %<>% factor(c("Normal","1","2","3","4"))
table(object$X4cluster_normal)

csv_file_names <- list.files("output/20230203",pattern = "MCL_Normal_51-FC0.25",full.names = TRUE)
degs <- pbapply::pblapply(csv_file_names, function(x) {
    tmp <- read.csv(x,header = TRUE,row.names = 1)
    tmp$cluster <- sub("\\.csv","",x) %>% sub(".*_","",.)
    tmp$genes <- rownames(tmp)
    tmp
    }) %>%
    bind_rows()
degs <- degs[grep("^MT-",degs$gene,invert = TRUE),]
Top_n = 40
top1 <- degs %>% group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    distinct() %>%
    top_n(Top_n, avg_log2FC)


markers <- c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11")
markers = markers[markers %in% rownames(object)]

features = c(as.character(top1$genes),
             tail(VariableFeatures(object = object), 2),
             markers)
table(duplicated(features))
sub_obj <- object[features,]
Idents(sub_obj) = "X4cluster_normal"
table(Idents(sub_obj))
sub_obj %<>% ScaleData(features=features)
featuresNum <- make.unique(features, sep = ".")
sub_obj %<>% MakeUniqueGenes(features = features)
dim(sub_obj[["SCT"]]@scale.data)
table(duplicated(rownames(sub_obj[["SCT"]]@scale.data)))
table(featuresNum %in% rownames(sub_obj[["SCT"]]@scale.data))
DoHeatmap.2(sub_obj, features = featuresNum, #Top_n = Top_n,
            do.print=T,
            angle = 45,
            group.by = c("X4cluster_normal","orig.ident"), group.bar = T,
            group1.colors = c('#B2DF8A','#40A635','#FE8205','#8861AC','#E83C2D'),
            group2.colors= Singler.colors,
            title.size = 16, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.02, label=F, cex.row= 5,
            width=14, height=12,
            file.name = paste0("Fig S2. top",Top_n,"_X4cluster_normal_orig.ident.jpeg"),
            title =  paste("Top",Top_n, "DE genes in 4 B/MCL clusters"),
            save.path = save_path,
            nrow = 5, ncol = 8, design = c(patchwork::area(1, 1, 5, 5),
                                           patchwork::area(3, 6, 3, 8)))

DoHeatmap.2(sub_obj, features = featuresNum, #Top_n = Top_n,
            do.print=T,
            angle = 45,
            group.by = c("X4cluster_normal","patient"), group.bar = T,
            group1.colors = c('#B2DF8A','#40A635','#FE8205','#8861AC','#E83C2D'),
            group2.colors= Singler.colors,
            title.size = 16, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.02, label=F, cex.row= 5,
            width=14, height=12,
            file.name = paste0("Fig S2. top",Top_n,"_X4cluster_normal_patient.jpeg"),
            title =  paste("Top",Top_n, "DE genes in 4 B/MCL clusters"),
            save.path = save_path,
            nrow = 5, ncol = 8, design = c(patchwork::area(1, 1, 5, 5),
                                           patchwork::area(3, 6, 3, 7)))


Top_n = 150
top <- degs %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC),.by_group = TRUE) %>%
    distinct() %>%
    top_n(Top_n, avg_log2FC) %>%
    ungroup()
table(top$cluster)

#write.table(top,file = paste0("output/20230131/top",Top_n,"_degs_X4cluster_v1.csv"),row.names = FALSE,quote = FALSE, col.names = FALSE)

features = c(as.character(top$genes),
             tail(VariableFeatures(object = object), 2),
             markers)
table(duplicated(features))
sub_obj <- object[features,]
Idents(sub_obj) = "X4cluster_normal"
sub_obj %<>% ScaleData(features=features)
featuresNum <- make.unique(features, sep = ".")
sub_obj %<>% MakeUniqueGenes(features = features)
dim(sub_obj[["SCT"]]@scale.data)
table(duplicated(rownames(sub_obj[["SCT"]]@scale.data)))
table(featuresNum %in% rownames(sub_obj[["SCT"]]@scale.data))
DoHeatmap.2(sub_obj, features = featuresNum, #Top_n = Top_n,
            do.print=T,
            angle = 45,
            group.by = c("X4cluster_normal","orig.ident"), group.bar = T,
            group1.colors = c('#B2DF8A','#40A635','#FE8205','#8861AC','#E83C2D'),
            group2.colors= Singler.colors,
            title.size = 16, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.02, label=F, cex.row= 1,
            width=14, height=12,
            file.name = paste0("Fig S2. top",Top_n,"_X4cluster_normal_orig.ident.jpeg"),
            title =  paste("Top",Top_n, "DE genes in 4 B/MCL clusters"),
            save.path = save_path,
            nrow = 5, ncol = 8, design = c(patchwork::area(1, 1, 5, 5),
                                           patchwork::area(3, 6, 3, 8)))

DoHeatmap.2(sub_obj, features = featuresNum, #Top_n = Top_n,
            do.print=T,
            angle = 45,
            group.by = c("X4cluster_normal","patient"), group.bar = T,
            group1.colors = c('#B2DF8A','#40A635','#FE8205','#8861AC','#E83C2D'),
            group2.colors= Singler.colors,
            title.size = 16, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.02, label=F, cex.row= 1,
            width=14, height=12,
            file.name = paste0("Fig S2. top",Top_n,"_X4cluster_normal_patient.jpeg"),
            title =  paste("Top",Top_n, "DE genes in 4 B/MCL clusters"),
            save.path = save_path,
            nrow = 5, ncol = 8, design = c(patchwork::area(1, 1, 5, 5),
                                           patchwork::area(3, 6, 3, 7)))

openxlsx::write.xlsx(list("top40" = top1,"top150" = top), file =  paste0(path,"DEGs_X4cluster_normal.xlsx"),
                     colNames = TRUE,rownames = T,borders = "surrounding")

#write.csv(top,file = paste0("output/20230131/top",Top_n,"_degs_X4cluster_normal_v1.csv"),row.names = FALSE,quote = FALSE)
openxlsx::write.xlsx(top, file =  paste0("output/20230131/top",Top_n,"_degs_X4cluster_normal_v1.xlsx"),
                     colNames = TRUE,rownames = T,borders = "surrounding")

#==== Figure 3S-F venn diagram ===========
path <- "Yang/PALIBR/51_samples/Fig. 3/"
opts <- data.frame("choose" = rep(c("X4clusters","X4cluster_vs_Normal"), each = 2),
                   "LogFC" = c(0,0.1,0,0.25), stringsAsFactors = F)
save_path <- "Yang/PALIBR/51_samples/Fig. S3/"

choose = "X4cluster"
gde.all = read.csv(file = paste0(path,"/MCL_only_",choose,"_51-FC0.csv"),
                   row.names = 1, stringsAsFactors=F)
gde.all = gde.all[gde.all$p_val_adj <0.01 ,]

for(value in c(0.5,0.25,0.1)){
        pos_genes <- eulerr(gde.all,shape =  "circle",#key = c("C1","C2","C3","C4","B_cells"),
                            cut_off = "avg_log2FC", cut_off_value = value,do.print = T,return.raw = T,
                            save.path = save_path, file.name =paste0("Fig3F_Venn_log2FC_",value,".jpeg"))
        euler_df <- eulerr::euler(pos_genes,shape = "circle")
        pos_genes_list <- as.list(euler_df$original.values)
        names(pos_genes_list) %<>% paste(":",pos_genes_list)
        id <- eulerr:::bit_indexr(4)
        for (i in nrow(id):1) {
                pos_genes_list[[i]] = Reduce(intersect, pos_genes[id[i,]])  %>%
                        setdiff(Reduce(union, pos_genes[!id[i,]]))
        }
        pos_genes_df <- list2df(pos_genes_list)
        write.xlsx(pos_genes_df, asTable = F,
                   file = paste0(save_path,"Fig3F_Venn_postive_shared_gene_list_FC",value,".xlsx"),
                   borders = "surrounding")
}

#======= dendrogram ============
df_list <- split(DEG,DEG[,"sample"])
comm_shared_genes <- lapply(df_list, function(df) df[(df$avg_logFC > 0),]) %>%
        lapply(function(df) df[(abs(df[,"avg_logFC"]) > 0.25),"gene"]) %>%
        Reduce(intersect, .)
length(comm_shared_genes)

B_cells_MCL = readRDS(file = "data/MCL_41_B_20200225.rds")
Samples = c("Pt10_LN2Pd","Pt13_BMA1","Pt25_SB1","Pt25_AMB25Pd")

Idents(B_cells_MCL) = "orig.ident"
B_cells_MCL %<>% subset(idents = Samples)
exp = AverageExpression(B_cells_MCL,assays = "SCT",features = comm_shared_genes)
exp = exp$SCT


require(gtools)
require(RColorBrewer)
require(made4)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
jpeg(paste0(read_path, "/F3_dendrogram_spearman.jpeg"),units="in", width=4.5, height=7.5,res=600)
heatplot.1(exp,cor.method = "spearman",cexCol=0.6)
dev.off()

jpeg(paste0(read_path, "/F3_dendrogram_pearson.jpeg"),units="in", width=4.5, height=7.5,res=600)
heatplot.1(exp,cor.method = "pearson",cexCol=0.6)
dev.off()


ccgenes_shared_mtr =  exp %>% tibble::rownames_to_column("gene") %>%
        filter(gene %in% c(unlist(cc.genes,use.names = F),"CCND1")) %>%
        tibble::column_to_rownames("gene")

cc_genes = c(unlist(cc.genes,use.names = F),"CCND1","CDK4") %>% .[. %in% rownames(B_cells_MCL)]
exp = AverageExpression(B_cells_MCL,assays = "SCT",features = cc_genes)
exp = exp$SCT
write.csv(exp,file = paste0(read_path, "/FS3_402_genes.csv"))
write.csv(exp,file = paste0(read_path, "/FS3_96_genes.csv"))

cols.gentleman <- function() {
        library(RColorBrewer)
        hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
        return(rev(hmcol))
}


#=====S4B heatmap Pt25 ======================
path <- "Yang/PALIBR/51_samples/Fig. S4/"
if(!dir.exists(path)) dir.create(path, recursive = T)

B_cells_MCL = readRDS(file = "data/MCL_SCT_51_20210724.rds")
DefaultAssay(B_cells_MCL) = "SCT"
B_cells_MCL <- subset(B_cells_MCL, subset = orig.ident %in% c("Pt25_SB1","Pt25_AMB25","N01"))
B_cells_MCL@meta.data[B_cells_MCL$orig.ident %in% "N01","X4cluster"] = "Normal"


B_cells_MCL = subset(B_cells_MCL, subset =  Doublets == "Singlet"
                     & X4cluster  %in% c("1","2","3","4","Normal")
)
B_cells_MCL$orig.ident %<>% droplevels()

sub_ob1 <- subset(B_cells_MCL, subset = orig.ident %in% c("Pt25_SB1","N01"))
(group.by = c("Normal","1","2","3","4"))

sub_ob1@meta.data$X4cluster_normal %<>% factor(levels = group.by)

Idents(sub_ob1) = "X4cluster"
Idents(sub_ob1) %<>% factor(levels = group.by)
table(Idents(sub_ob1))

ident_list = list(list("Normal","1"),
                  list("Normal","2"),
                  list("Normal","3"),
                  list("Normal","4"),
                  list("1","Normal"),
                  list("2","Normal"),
                  list("3","Normal"),
                  list("4","Normal")
                  )
MCL_markers <- pblapply(ident_list, function(ident){
        mark = FindMarkers_UMI(sub_ob1,ident.1 = ident[[1]],
                        ident.2 = ident[[2]],
                        logfc.threshold = 0.75,
                        only.pos = T,
                        test.use = "MAST",
                        latent.vars = "nFeature_SCT")
        mark$cluster = paste0(ident[[1]],"_vs_",ident[[2]])
        mark$gene = rownames(mark)

        return(mark)
        }
)

MCL_markers %<>% bind_rows()
write.csv(MCL_markers,paste0(path,"Pt25_SB1_X4Cluster_normal_markers.csv"))
table(MCL_markers$cluster)
markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
(MT_gene <- grep("^MT-",MCL_markers$gene))
if(length(MT_gene) >0 ) MCL_markers = MCL_markers[-MT_gene,]
Top_n = 40

top = MCL_markers %>%   filter(!stringr::str_detect(gene, "^MT-")) %>%
        group_by(cluster) %>%
        top_n(50, avg_log2FC)
#top = top[order(top$cluster),]
write.csv(top,paste0(path,"Top40_Pt25_SB1_X4Cluster_normal_markers.csv"))
features = top$gene#c(as.character(top$gene),
             #tail(VariableFeatures(object = sub_ob1), 2),
            # markers)
#DoHeatmap.1======
# heatmap
featuresNum <- features# make.unique(features, sep = ".")
exp = AverageExpression(sub_ob1,features = features,
                        assays = "SCT") %>% .$SCT
#exp %<>% MakeUniqueGenes(features = features)
exp %<>% t %>% scale %>% t
#exp[tail(VariableFeatures(object = B_cells_MCL), 2),] =0

(group.by = c("Normal","1","2","3","4"))
group.by %<>% factor(levels = group.by)
DoHeatmap.matrix(exp, features = featuresNum,
                 group.by = group.by,
                 size = 6,angle = 0,label =F,
                 draw.lines =F, raster = FALSE,
                 pal_gsea = FALSE,
                 width=4, height=9,res=600,no.legend = F,
                 cex.row=5,
                 group.colors = c('#B2DF8A','#40A635','#FE8205','#8861AC','#E83C2D'),
                 do.print = T,
                 save.path = path,
                 file.name = paste0("Heatmap_top40_Pt25_SB1_X4Cluster.jpeg")
)

#  remove RP
top = MCL_markers %>%   filter(!stringr::str_detect(gene, "^MT-")) %>%
        filter(!stringr::str_detect(gene, "^RPL")) %>%
        filter(!stringr::str_detect(gene, "^RPS")) %>%
        group_by(cluster) %>%
        top_n(50, avg_log2FC)
unique(top$cluster)
#top = top[order(top$cluster),]
write.csv(top,paste0(path,"Top40_Pt25_SB1_X4Cluster_normal_markers_rmRP.csv"))
features = as.character(top$gene)
featuresNum <- features#make.unique(features, sep = ".")
exp = AverageExpression(sub_ob1,features = features,
                        assays = "SCT") %>% .$SCT
#exp %<>% MakeUniqueGenes(features = features)
exp %<>% t %>% scale %>% t
#exp[tail(VariableFeatures(object = sub_object), 2),] =0

(group.by = c("Normal","1","2","3","4"))
group.by %<>% factor(levels = group.by)
DoHeatmap.matrix(exp, features = featuresNum,
                 group.by = group.by,
                 size = 6,angle = 0,label =F,
                 draw.lines =F, raster = FALSE,
                 pal_gsea = FALSE,
                 width=4, height=9,res=600,no.legend = F,
                 cex.row=5,
                 group.colors = c('#B2DF8A','#40A635','#FE8205','#8861AC','#E83C2D'),
                 do.print = T,
                 save.path = path,
                 file.name = paste0("Heatmap_top40_Pt25_SB1_X4Cluster_rmRP.jpeg")
)

#-------Pt25_SB1 vs Pt25_AMB25----------------------------
sub_ob2 <- subset(B_cells_MCL, subset = orig.ident %in% c("Pt25_SB1","Pt25_AMB25"))

(group.by <- mapply(paste, rep(c("Pt25_SB1","Pt25_AMB25"),each = 4),
       rep(c("1","2","3","4"),2)))
sub_ob2$orig.ident_cluster = paste(sub_ob2$orig.ident,sub_ob2$X4cluster)
sub_ob2$orig.ident_cluster %<>% factor(levels = group.by)


Idents(sub_ob2) = "orig.ident_cluster"
table(Idents(sub_ob2))

ident_list = list(list("Pt25_SB1 1","Pt25_AMB25 1"),
                  list("Pt25_SB1 2","Pt25_AMB25 2"),
                  list("Pt25_SB1 3","Pt25_AMB25 3"),
                  list("Pt25_SB1 4","Pt25_AMB25 4"),
                  list("Pt25_AMB25 1","Pt25_SB1 1"),
                  list("Pt25_AMB25 2","Pt25_SB1 2"),
                  list("Pt25_AMB25 3","Pt25_SB1 3"),
                  list("Pt25_AMB25 4","Pt25_SB1 4")
)
MCL_markers <- pblapply(ident_list, function(ident){
        mark = FindMarkers_UMI(sub_ob2,ident.1 = ident[[1]],
                               ident.2 = ident[[2]],
                               logfc.threshold = 0.75,
                               only.pos = T,
                               test.use = "MAST",
                               latent.vars = "nFeature_SCT")
        mark$cluster = paste0(ident[[1]],"_vs_",ident[[2]])
        mark$gene = rownames(mark)

        return(mark)
}
)

MCL_markers %<>% bind_rows()
write.csv(MCL_markers,paste0(path,"Pt25_SB1_AMB25_markers.csv"))
table(MCL_markers$cluster)
markers <- FilterGenes(sub_ob2,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
(MT_gene <- grep("^MT-",MCL_markers$gene))
if(length(MT_gene) >0 ) MCL_markers = MCL_markers[-MT_gene,]
Top_n = 40

top = MCL_markers %>% group_by(cluster) %>%
        top_n(40, avg_log2FC)
#top = top[order(top$cluster),]
write.csv(top,paste0(path,"Top40_Pt25_SB1_AMB25_markers.csv"))
features = as.character(top$gene)
#DoHeatmap.1======
# heatmap
featuresNum <- features# make.unique(features, sep = ".")
exp = AverageExpression(sub_ob2,features = features,
                        assays = "SCT") %>% .$SCT
#exp %<>% MakeUniqueGenes(features = features)
exp %<>% t %>% scale %>% t
#exp[tail(VariableFeatures(object = sub_object), 2),] =0

(group.by %<>% factor(levels = group.by))
DoHeatmap.matrix(exp, features = featuresNum,
                 group.by = group.by,
                 size = 6,angle = 0,label =F,
                 draw.lines =F, raster = FALSE,
                 pal_gsea = FALSE,
                 width=5, height=9,res=600,no.legend = F,
                 cex.row=5,
                 group.colors = rep(c('#40A635','#FE8205','#8861AC','#E83C2D'),time =2),
                 do.print = T,
                 save.path = path,
                 file.name = paste0("Heatmap_top40_Pt25_SB1_AMB25.jpeg")
)

#  remove RP
top = MCL_markers %>%
        filter(!stringr::str_detect(gene, "^RPL")) %>%
        filter(!stringr::str_detect(gene, "^RPS")) %>%
        group_by(cluster) %>%
        top_n(40, avg_log2FC)
unique(top$cluster)
#top = top[order(top$cluster),]
write.csv(top,paste0(path,"Top40_Pt25_SB1_AMB25_markers.csv"))
features = c(as.character(top$gene))#,
             #tail(VariableFeatures(object = sub_object), 2),
             #markers)
featuresNum <- features#make.unique(features, sep = ".")
exp = AverageExpression(sub_ob2,features = features,
                        assays = "SCT") %>% .$SCT
#exp %<>% MakeUniqueGenes(features = features)
exp %<>% t %>% scale %>% t
#exp[tail(VariableFeatures(object = sub_object), 2),] =0

group.by %<>% factor(levels = group.by)
DoHeatmap.matrix(exp, features = featuresNum,
                 group.by = group.by,
                 size = 6,angle = 0,label =F,
                 draw.lines =F, raster = FALSE,
                 pal_gsea = FALSE,
                 width=5, height=9,res=600,no.legend = F,
                 cex.row=5,
                 group.colors = rep(c('#40A635','#FE8205','#8861AC','#E83C2D'),time =2),
                 do.print = T,
                 save.path = path,
                 file.name = paste0("Heatmap_top40_Pt25_SB1_X4Cluster_rmRP.jpeg")
)
#==============================================





jpeg(paste0(read_path, "/FS3_ccgenes_dendrogram_spearman~.jpeg"),units="in", width=4.5, height=14,res=600)
heatplot.1(exp,cor.method = "spearman",cexCol=0.7)
dev.off()

jpeg(paste0(read_path, "/FS3_ccgenes_dendrogram_spearman.jpeg"),units="in", width=4.5, height=7.5,res=600)
heatplot.1(exp,cor.method = "spearman",cexCol=0.6)
dev.off()

jpeg(paste0(read_path, "/FS3_ccgenes_dendrogram_pearson~.jpeg"),units="in", width=4.5, height=7.5,res=600)
heatplot.1(exp,cor.method = "pearson",cexCol=0.6)
dev.off()

write.csv(ccgenes_shared_mtr,file = paste0(read_path, "/FS3_ccgenes.csv"))
