library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
library(readxl)
library(pheatmap)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
source("R/util.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
# 3.1.1 load bulk RNA-seq data
df_samples <- readxl::read_excel("doc/190126_scRNAseq_info.xlsx",
                                 sheet = "bulk")
colnames(df_samples) = tolower(colnames(df_samples))
df_samples$bulk.rnaseq.id = gsub(" ","_",df_samples$bulk.rnaseq.id)

MCL_6669 <- read_excel("data/RNAseq/GRCF-Elemento-MDL-6669-expression.genes.xls")
MCL_6392 <- read_excel("data/RNAseq/GRCF-Elemento-MDL-6392-expression.genes.xls")

(sample_6669 <- colnames(MCL_6669)[colnames(MCL_6669) %in% df_samples$bulk.rnaseq.id])
MCL_6669 = MCL_6669[,c("Gene_name",sample_6669)]

(sample_6392 <- colnames(MCL_6392)[colnames(MCL_6392) %in% df_samples$bulk.rnaseq.id])
MCL_6392 = MCL_6392[,c("Gene_name",sample_6392)]

colnames(MCL_6392)[-1] = paste0(colnames(MCL_6392)[-1],"_bulk")

MCL_bulk <- merge(MCL_6669, MCL_6392, by = "Gene_name")
MCL_bulk <- RemoveDup(MCL_bulk)
log_MCL_bulk <- log1p(MCL_bulk)
testMMM(log_MCL_bulk)
# 3.1.2 load sc RNA-seq data
(load(file="data/MCL_Harmony_24_20190128.Rda"))
object %<>% SetAllIdent("orig.ident")
MCL_sc <- AverageExpression(object,use.raw = T)
testMMM(MCL_scRNA)

# 3.1.3 merge bulk and sc
MCL_bulk_sc <- merge(log1p(MCL_bulk), log1p(MCL_scRNA), by="row.names")
table(duplicated(MCL_bulk_sc$Row.names))
rownames(MCL_bulk_sc) = MCL_bulk_sc$Row.names
MCL_bulk_sc = MCL_bulk_sc[,-which(colnames(MCL_bulk_sc) %in% "Row.names")]
testMMM(MCL_bulk_sc)
dim(MCL_bulk_sc)

var.genes <- rownames(MCL_bulk_sc)[rownames(MCL_bulk_sc) %in% B_cells_MCL@var.genes]
length(var.genes)
MCL_bulk_sc = MCL_bulk_sc[var.genes,]

MCL_bulk_sc <- apply(MCL_bulk_sc, 2, function(x) (x- mean(x))/sd(x))
boxplot(MCL_bulk_sc)

# 3.1.4 Spearman correlation  ==================
y <- cor(MCL_bulk_sc, method="spearman")
diag(y) <-NA
ident_num <- ncol(as.matrix(MCL_bulk))
short_y <- y[(ident_num+1):nrow(y),1:ident_num]


hc_c <- hclust(as.dist(1-cor(short_y, method="spearman")), method="complete")
hc_r <- hclust(as.dist(1-cor(t(short_y), method="spearman")), method="complete")
jpeg(paste0(path,"cluster_heatmap.jpeg"), units="in", width=10, height=7,res=600)
print(pheatmap(short_y,cex=.9,
         cluster_rows= hc_r,
         cluster_cols = hc_c,
         fontsize_row = 15,
         fontsize_col = 12,
         fontsize = 20,
         main = ""))
dev.off()

# 3.2.1 CD19+ sort cells ===========================
df_samples <- readxl::read_excel("doc/190126_scRNAseq_info.xlsx",
                                 sheet = "bulk_MCL")
colnames(df_samples) = tolower(colnames(df_samples))
MCL_181120 <- read_excel("data/RNAseq/181120 MCL WTS.xlsx")
MCL_190311 <- read_excel("data/RNAseq/190311 PALIBR samples WTS for scRNAseq.xlsx")
MCL_AA13 <- read_excel("data/RNAseq/190311 Pt AA13 WTS .xlsx")
MCL_AFT <- read_excel("data/RNAseq/190311 AFT-03 AFT-04 WTS.xlsx")

MCL_181120 = CleanUp(MCL_181120)
MCL_190311 = CleanUp(MCL_190311)
MCL_AA13 = CleanUp(MCL_AA13)
MCL_AFT = CleanUp(MCL_AFT)

MCL_list <- lapply(list(MCL_181120,MCL_190311,MCL_AA13,MCL_AFT), 
                   function(x) {
                       x = as.data.frame(x)
                       x$rn = rownames(x)
                       x})
MCL_bulk <- plyr::join_all(MCL_list, by = 'rn', type = 'inner')
rownames(MCL_bulk) = MCL_bulk$rn
MCL_bulk = MCL_bulk[,-which(colnames(MCL_bulk) == "rn")]
colnames(MCL_bulk) = df_samples$sample[match(colnames(MCL_bulk),df_samples$sample.id)]
MCL_bulk = MCL_bulk[,!duplicated(colnames(MCL_bulk))]
write.csv(MCL_bulk, file="data/RNAseq/MCL_bulk_20190313.csv")

# 3.2.2 AverageExpression scRNA ===========================
# generate B_cells_MCL from Differential_analysis_B.R
B_cells_MCL %<>% SetAllIdent(id="orig.ident")
MCL_scRNA <- AverageExpression(B_cells_MCL,use.raw = T)

# 3.3 run Seruat pipeline=================
bulk <- CreateSeuratObject(raw.data = MCL_bulk, project = "bulk", min.cells = 1)
bulk@meta.data$type = bulk@meta.data$orig.ident
bulk@meta.data$orig.ident = rownames(bulk@meta.data)
bulk@meta.data$group = df_samples$group[match(rownames(bulk@meta.data),df_samples$sample)]
bulk@meta.data$label = df_samples$label[match(rownames(bulk@meta.data),df_samples$sample)]

sc <- CreateSeuratObject(raw.data = MCL_scRNA, project = "scRNAseq", min.cells = 1)
sc@meta.data$type = sc@meta.data$orig.ident
sc@meta.data$orig.ident = rownames(sc@meta.data)
sc@meta.data$label = df_samples$label[match(rownames(sc@meta.data),df_samples$sample)]

bulk_sc <- MergeSeurat(bulk, sc, add.cell.id1 = "bulk", add.cell.id2 = "sc")
#=======
bulk %<>% NormalizeData %>% 
    ScaleData %>% 
    FindVariableGenes(do.plot = FALSE) %>%
    RunPCA(pcs.compute = 45)
PCElbowPlot(bulk, num.pc = 45)
bulk %<>% FindClusters(resolution = 1,dims.use = 1:40) %>%
    RunTSNE(perplexity=15,dims.use = 1:40)
bulk %<>% SetAllIdent("res.1")
TSNEPlot.1(bulk,do.label =T,no.legend = T,label.repel = T,pt.size = 4,
           text.repel = F, do.print=T)
bulk %<>% SetAllIdent("orig.ident")
TSNEPlot.1(bulk,do.label =F,no.legend = T,label.repel = T,pt.size = 3,
           text.repel = F, do.print=T)
FeaturePlot(bulk,features.plot = "CCND1",pt.size = 4)

marker.list <- HumanGenes(bulk,c("CD19","CD5","CCND1","MS4A1","CDK4","SOX11"))
p <- lapply(marker.list, function(marker) {
    SingleFeaturePlot.1(object = bulk, feature = marker,pt.size = 5,
                        threshold=NULL,legend.title ="")+ #log(FPKM)
        ggtitle(paste0(marker))+
        theme(plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
})

jpeg(paste0(path,"B_cells_bulk_RNAseq.jpeg"), units="in", width=10, height=7,res=600)
print(do.call(plot_grid, p)+ ggtitle("MCL markers")+
          theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold")))
dev.off()

jpeg(paste0(path,"bulk_CCND1.jpeg"), units="in", width=10, height=7,res=600)
CCND1  <- bulk@data["CCND1",]
hist(CCND1,breaks = 35)
dev.off()

sort(CCND1) %>% .[20:39] %>% kable %>% kable_styling()

keep <- names(CCND1)[CCND1 > 0.5]
(keep = keep[!(keep %in% "6PB1")])
write.csv(MCL_bulk[,keep],file="data/RNAseq/MCL_bulk.csv")

# prepare MCL-stage
bulk@meta.data$stage = bulk@meta.data$res.1
CCND1_high <- names(CCND1)[CCND1 > 2]
CCND1_median <- names(CCND1)[CCND1 > 0.5 & CCND1 < 2]
CCND1_low <- names(CCND1)[CCND1 < 0.5]

bulk@meta.data$stage[(bulk@meta.data$orig.ident %in% CCND1_high)] = "CCND1_high"
bulk@meta.data$stage[(bulk@meta.data$orig.ident %in% CCND1_median)] = "CCND1_median"
bulk@meta.data$stage[(bulk@meta.data$orig.ident %in% CCND1_low)] = "CCND1_low"

bulk %<>% SetAllIdent("stage")
TSNEPlot.1(bulk,do.label =T,no.legend = T,label.repel = T,pt.size = 3,
           text.repel = F, do.print=T)
