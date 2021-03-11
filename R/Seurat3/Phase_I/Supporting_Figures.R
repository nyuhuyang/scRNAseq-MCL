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
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- "Yang/Figure 2S/Supplementary Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data

(load(file="data/MCL_41_harmony_20200225.Rda"))
df_samples <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
object$orig.ident %<>% factor(levels = df_samples$`sample name`)

mean(object$nCount_SCT)
mean(object$nFeature_SCT)

# Extend Data
Idents(object) = "cell.types"
object %<>% sortIdent()

cell_Freq <- table(Idents(object)) %>% as.data.frame
cell_Freq$Percent <- prop.table(cell_Freq$Freq) %>% scales::percent()

cell_Freq = cell_Freq[order(cell_Freq$Var1),]
cell_Freq$col = ExtractMetaColor(object)
cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")
cell_Freq$Cell_Type %<>% gsub("_"," ",.)

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=6, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,x.text.angle = 45,
          ylab = "Cell Number",
          label = "Percent",
          sort.val = "desc",
          width = 1, size = 0.5,
          title = "Numbers of major cell types in total 41 samples")+NoLegend()+
        theme(plot.title = element_text(hjust = 0.5,size=15))+
        scale_y_continuous(expand = c(0, 0), limits = c(0,max(cell_Freq$Cell_Number)+1500))
dev.off()

#=======
Seurat_list <- SplitObject(object, split.by = "orig.ident")

QC_list <- as.data.frame(df_samples)
QC_list["cell.number"] <- sapply(Seurat_list, function(x) length(colnames(x)))
QC_list["mean.nUMI"] <- sapply(Seurat_list, function(x) mean((x$nCount_SCT)))
QC_list["mean.nGene"] <- sapply(Seurat_list, function(x) mean((x$nFeature_SCT)))
QC_list["mean.percent.mt"] <- sapply(Seurat_list, function(x) mean((x$percent.mt)))

write.csv(QC_list,paste0(path,"QC_list.csv"))
#QC.list %>% kable() %>% kable_styling()

remove(Seurat_list);GC()

QC_list <- read.csv(paste0(path,"QC_list.csv"), stringsAsFactors = F)
jpeg(paste0(path,"mean.Reads.per.Cell.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(QC_list, x = "submitter", y= "mean.reads.per.cell",
         title = "Mean reads per cell in each scRNA-seq",
         xlab = "",ylab = "Mean Reads per Cell",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5,
         yscale = "log10"
         )+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,350000))+
                 theme(plot.title = element_text(hjust = 0.5,size=15),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank())
dev.off()

jpeg(paste0(path,"cell.number.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(QC_list, x = "submitter", y= "cell.number",
         title = "Cell number in each scRNA-seq",
         xlab = "",ylab = "Cell Number",
         add = c("jitter","mean_sd"),
         #ylim = c(0, max(QC_list$cell.number)+1000),
         draw_quantiles = 0.5,
         yscale = "log10"
         )+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,6000))+
        theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()

jpeg(paste0(path,"UMI.per.Cell.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(QC_list, x = "submitter", y= "mean.nUMI",
         title = "Mean transcripts per cell in each scRNA-seq",
         xlab = "",ylab = "Mean UMI per Cell",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5,
         yscale = "log10"
         )+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,12000))+
        theme(plot.title = element_text(hjust = 0.5,size=13),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()

jpeg(paste0(path,"nGene.per.Cell.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(QC_list, x = "submitter", y= "mean.nGene",
         title = "Mean genes per cell in each scRNA-seq",
         xlab = "",ylab = "Mean genes per Cell",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5,
         yscale = "log10"
         )+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,4000))+
        theme(plot.title = element_text(hjust = 0.5,size=13),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()

jpeg(paste0(path,"percent.mt.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(QC_list, x = "submitter", y= "mean.percent.mt",
         title = "Mean mitochondrial gene % in each scRNA-seq",
         xlab = "",ylab = "Mean mitochondrial gene percentages",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5,
         #yscale = "log10"
         )+
        scale_y_continuous(expand = c(0, 0), limits = c(0,30))+
        theme(plot.title = element_text(hjust = 0.5,size=13),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()
# ================
Idents(object) = "orig.ident"

(mito.features <- grep(pattern = "^MT-", x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = "^MT-")
g2 <- lapply(c("nFeature_SCT", "nCount_SCT", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=8),legend.position="none")
})

jpeg(paste0(path,"S2_nGene.jpeg"), units="in", width=7, height=5,res=600)
g2[[1]]+ggtitle("Distribution of Gene Number per Cells")+
        xlab("scRNA-seq samples")+
        ylab("Gene Number")+
        scale_y_log10()+
        #ylim(0,max(object$nFeature_SCT)+100)+
        theme(plot.title = element_text(face = 'plain'))
dev.off()

jpeg(paste0(path,"S2_nUMI.jpeg"), units="in", width=7, height=5,res=600)
g2[[2]]+ggtitle("Distribution of transcripts per Cells")+
        xlab("scRNA-seq samples")+
        ylab("UMI per Cell")+
        scale_y_log10()+
        #ylim(0,max(object$nCount_SCT)+1000)+
        theme(plot.title = element_text(face = 'plain'))
dev.off()

jpeg(paste0(path,"S2_mito.jpeg"), units="in", width=7, height=5,res=600)
g2[[3]]+ggtitle("Distribution of mitochondrial gene percentage per Cells")+
        xlab("scRNA-seq samples")+
        ylab("Mitochondrial gene percentage %")+
        theme(plot.title = element_text(face = 'plain'))
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

