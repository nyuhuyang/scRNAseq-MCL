# ######################################################################
invisible(lapply(c("Seurat","dplyr","ggpubr","openxlsx","Hmisc",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

#======1.2 load  Seurat =========================
B_cells_MCL = readRDS(file = "data/MCL_41_B_20200225.rds")
B_cells_MCL$orig.ident %<>% gsub("N01|N02|N03","Normal",.)
Idents(B_cells_MCL) = "orig.ident"
samples = c("Pt25_SB1","Pt25_1","Pt25_1_8","Pt25_24","Pt25_25Pd","Pt25_AMB25Pd",
            "Pt10_LN2Pd","All_samples","PtB13_Ibp","PtB13_Ib1","PtB13_IbR","Normal",
            "Pt11_LN1", "Pt17_LN1","PtU01","PtU02","PtU03","PtU04")

test_genes <- c("EZH2","EZH1","CCND1","E2F1","PCNA","IRF4","PIK3IP1","HLA-DPA1","MCM7",
                   "HLA-DPB1","HLA-A","HLA-B","BCL6","MYC","MEF2B","CDKN1A","NFKB2","MAP3K8",
                   "FOXM1","RELB","POLR2M","CRBN","IKZF1","IKZF3","MBOAT7")
test_genes = test_genes[test_genes %in% rownames(B_cells_MCL)]

s = samples[args]
sub_object <- subset(B_cells_MCL, idents = s)
## Column clustering (adjust here distance/linkage methods to what you need!)
sub_object <- FindVariableFeatures(object = sub_object, selection.method = "vst",
                               num.bin = 20,nfeatures = 3000,
                               mean.cutoff = c(0.1, 8), 
                               dispersion.cutoff = c(1, Inf))
test_genes <- unique(c(VariableFeatures(sub_object),test_genes))
y = sub_object[["SCT"]]@data[test_genes,]
system.time(cor_res <- Hmisc::rcorr(t(as.matrix(y)), type="spearman"))
cor_res$r[is.na(cor_res$r)] = 0

save.path <- paste0(path, "Cor-",s)
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

jpeg(paste0(save.path,"heatmap-cor-",s,".jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(cor_res$r,trace="none")
dev.off()
df_list <- list()
for(i in seq_along(test_genes)){
        gene <- test_genes[i]
        df <- data.frame("correlation" = cor_res$r[gene,],
                         "p.value" = cor_res$P[gene,])
        df = df[!is.na(df$correlation),]
        df = df[(rownames(df) != gene),]
        df[,"log10.p.value"] = -log10(df$`p.value`)
        df = df[order(df$correlation),]
        df[,"genes"] = rownames(df)
        n = 10
        low_cor <- head(df$correlation,n)[n]
        high_cor <- tail(df$correlation,n)[n]
        jpeg(paste0(save.path,"cor-pvalue-",gene,"-",s,".jpeg"), units="in", width=10, height=7,res=600)
        g <- ggline(df, x = "correlation", y = "log10.p.value",
                    numeric.x.axis = TRUE,
                    ylab = "-log10(p-value)",
                    xlab = "Spearman Correlation",
                    label = "genes",             # Add point labels
                    label.select = list(criteria = paste("`x` <=",low_cor,"| `x` >=",high_cor)),           # show only labels for the top 2 points
                    repel = TRUE, 
                    title = paste("Top 3000 variable correlated with",gene,"in",s))+
                theme_minimal() + TitleCenter()
        print(g)
        dev.off()
        colnames(df)[3] = " -log10(p-value)"
        df_list[[i]] = df
        Progress(i,length(test_genes))
}
names(df_list) = test_genes
write.xlsx(df_list, file = paste0(path,"cor-pvalue-genes-",s,".xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
