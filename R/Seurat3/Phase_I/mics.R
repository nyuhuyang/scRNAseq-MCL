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

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

(csv_list = list.files("Yang/B_pairwise_heatmaps/DE_analysis_files"))
for(i in seq_along(csv_list)){
        gde.markers = read.csv(paste0(path,"DE_analysis_files/", csv_list[i]),
                               row.names = 1, stringsAsFactors = F)
        print(table(gde.markers$cluster1.vs.cluster2))
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        Top_n = 40
        top <-  gde.markers %>% group_by(cluster1.vs.cluster2) %>%
                top_n(Top_n, avg_logFC) %>% as.data.frame()
        write.csv(top,file = paste0(path,"DE_analysis_files/top40_", csv_list[i]))
}

read.path = "Yang/PALIBR Figures legends methods/"
X4clusters = read.csv(file = paste0(read.path,"Figure 2/Figure Sources/X4clusters/X4clusters_41-FC0.1.csv"))
test_genes = X4clusters[X4clusters$p_val_adj < 0.05, "gene"] %>% as.character %>% unique
length(test_genes)

Pt25 = read.csv(file = paste0(read.path,"Figure 5/X4clusters_vs_X4clusters/Pt25/DE_FC0_Pt25_AMB25Pd_Pt25_SB1.csv"))
keep = order(Pt25$p_val_adj) %>% head(80)
test_genes = unique(c(test_genes,as.character(Pt25[keep,"gene"])))
length(test_genes)

Pt10 = read.csv(file = paste0(read.path,"Figure 7/Figure Sources/top40_Pt10_LN2Pd_X4clusters_genes_heatmap.csv"))
test_genes = unique(c(test_genes,as.character(Pt10[,"gene"])))
length(test_genes)

csv_list <- list.files(path = paste0(read.path,"Figure 5/X4clusters_vs_X4clusters/PtB13/"),pattern = ".csv",
                       full.names = T)
PtB13 <- lapply(csv_list, read.csv) %>% bind_rows
keep = order(PtB13$p_val_adj) %>% head(160)
test_genes = unique(c(test_genes,as.character(PtB13[keep,"gene"])))
length(test_genes)

tested_genes <- c("IL32", "CD52", "LTB", "TXNIP", "PIK3R1", "TNFAIP3","IL10", "IL10RA", "IL10RB",
                  "EZH2","EZH1","CCND1","E2F1","PCNA","IRF4","PIK3IP1","HLA-DPA1","MCM7",
                  "HLA-DPB1","HLA-DRA","HLA-A","HLA-B","BCL6","MYC","MEF2B","CDKN1A","NFKB2","MAP3K8",
                  "FOXM1","RELB","POLR2M","CRBN","IKZF1","IKZF3","MBOAT7","MIF")

test_genes <- c(test_genes,"CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11",tested_genes) %>% unique()


length(test_genes)
write.csv(test_genes, file = paste0(path, "test_genes.csv"))
Top_genes <- c()
cell.types <- c("B","Myeloid","CD4T","CD8T","Monocytes")
i= 1
for(m in seq_along(cell.types)){
        cell.type = cell.types[m]
        folder_list <- list.dirs(path = paste0("output/20200612-Correlation/",cell.type,"-Correlation-pvalues"),
                               full.names = T)
        for(n in 2:length(folder_list)){
                xlsx <- list.files(path = folder_list[n],pattern = ".xlsx",full.names = T)
                for(s in 1:length(tested_genes)){
                        df = read.xlsx(xlsxFile = xlsx,sheet = tested_genes[s])
                        df = df[order(df$correlation),]
                        Top_genes = unique(c(Top_genes,head(df$genes,10),tail(df$genes,10)))
                        Progress(i, length(tested_genes)*(length(folder_list)-1)*length(cell.types))
                        i = i+1
                }

        }
}
