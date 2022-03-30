########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
library(Seurat)
library(dplyr)
library(tidyr)
#library(kableExtra)
library(magrittr)
library(ggplot2)
library(cowplot)
library(fgsea)
library(tibble)
library(ggsci)
library(progress)
set.seed(101)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

read_path = "output/20211214/response_X4clusters_vs_rest"
csv_list <- list.files(pattern = "FC0.1",path = read_path,full.names = T)
deg_list <- pbapply::pblapply(csv_list, function(x){
    tmp = read.csv(x,row.names = 1)
    tmp = tmp[order(tmp$avg_log2FC,decreasing = T), ]
    tmp
})
res =  bind_rows(deg_list)

hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v7.4.symbols.gmt")
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(hallmark) = gsub("\\_"," ",names(hallmark))
nfkb = read.table("data/200222 NFKB pathway gene list.txt",header = F) %>% pull
hallmark[["NF-kB pathway"]] =nfkb

# hallmark
Fgsea_res <- FgseaDotPlot(stats=res, pathways=hallmark,Rowv = T,
                          size = " -log10(padj)",
                          title = "enriched hallmark pathways",
                          plot.title = element_text(hjust = 1,size = 15),
                          axis.text.x = element_text(angle = 45, hjust = 1,size = 10),
                          width = 6,do.return = T)
pthy = c("MYC TARGETS V1", "MYC TARGETS V2","TNFA SIGNALING VIA NFKB","NF-kB pathway",
             "IL2 STAT5 SIGNALING","INTERFERON ALPHA RESPONSE","INTERFERON GAMMA RESPONSE")
for(pthy in pathways){
    pathway_res = Fgsea_res %>% filter(pathway %in% pthy)
    Score_df = pathway_res[,c("NES","cluster")]
    Score_df$response = gsub("_.*","",Score_df$cluster)
    Score_df$clusters = gsub(".*_","C",Score_df$cluster)
    Score_df %<>% pivot_wider(id_cols = -cluster,
                              values_from = "NES",
                              names_from = "clusters") %>%
        tibble::column_to_rownames("response")
    Score_df[is.na(Score_df)] =0

    prop_df = pathway_res[,c(" -log10(padj)","cluster")]
    prop_df$response = gsub("_.*","",prop_df$cluster)
    prop_df$clusters = gsub(".*_","C",prop_df$cluster)
    prop_df %<>% pivot_wider(id_cols = -cluster,
                             values_from = " -log10(padj)",
                             names_from = "clusters") %>%
        tibble::column_to_rownames("response")
    prop_df[is.na(prop_df)] =0
    features = paste0("C",1:4)[paste0("C",1:4) %in% colnames(prop_df)]
    g <- DotPlot.2(Score_df,prop_df,scale = FALSE,log.data = NULL,features = features,
                        score_title = "NES", pct.title = " -log10(padj)",
              scale.by = "size", scale.max  = 5,
              col.min = 0, exp.max = 5,dot.scale = 10)+
              ggtitle(pthy)+
        scale_x_discrete(position = "top")
    jpeg(paste0(path,pthy,"_fgsea_dotplot.jpeg"), units="in", width=6, height=5,res=600)
    print(g)
    dev.off()
    svMisc::progress(which(pathways %in% pthy)/length(pathways)*100)
}


openxlsx::write.xlsx(Fgsea_res,
                     file =  paste0(path,"hallmark_gsea.xlsx"),
                     colNames = TRUE,rowNames = F,borders = "surrounding")
write.csv(Fgsea_res,file = paste0(path,"hallmark_gsea.xlsx"))
