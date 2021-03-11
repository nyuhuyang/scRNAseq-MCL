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

cell.types = c("B and MCL cells","T_cells:CD8+", "T_cells:CD4+","Monocytes","NK_cells","Myeloid cells")
file.names =c("B", "CD8T", "CD4T","Monocytes","NK","Myeloid")
Samples <- c("All_samples","Normal","N04","PtU01","PtU02","PtU03","PtU04",
             "Pt2_30Pd","Pt10_LN2Pd","Pt11_LN1","Pt11_1",
             "Pt11_14","Pt11_28","Pt13_BMA1","Pt13_1a",
             "Pt13_1b","Pt16_3Pd","Pt17_LN1","Pt17_2",
             "Pt17_7","Pt17_12","Pt17_31","Pt19_BM2Pd",
             "Pt20_1","Pt25_SB1","Pt25_1","Pt25_1_8",
             "Pt25_24","Pt25_25Pd","Pt25_AMB25Pd","Pt27_1",
             "Pt27_1_8","Pt27_12","PtB13_Ibp","PtB13_Ib1",
             "PtB13_IbR","Pt28_LN1","Pt28_1","Pt28_4",
             "Pt28_28")
opts = data.frame("cell.types" = rep(cell.types, each = length(Samples)),
                  "file.names" = rep(file.names, each = length(Samples)),
                  "samples" = rep(Samples,      time = length(cell.types)),
                  stringsAsFactors = F)
print(cell.type <- opts$cell.types[args])
print(s <- opts$samples[args])
print(file <- opts$file.names[args])
heatmap = FALSE
jpg = FALSE
# need 128 GB memory
if(cell.type == "B and MCL cells"){
        sub_object = readRDS(file = "data/MCL_41_B_20200225.rds")
        sub_object$orig.ident %<>% gsub("N01|N02|N03","Normal",.)
} else {
        (load(file="data/MCL_41_harmony_20200225.Rda"))
        DefaultAssay(object) = "SCT"
        Idents(object) = "Doublets"
        object %<>% subset(idents = "Singlet")
        object@meta.data[object$SCT_snn_res.0.8 %in% c(5,16), "cell.types"] = "Monocytes"
        Idents(object) = "cell.types"
        sub_object <- subset(object, idents = cell.type)
        sub_object$orig.ident %<>% gsub("N01|N02|N03","Normal",.)
}
Idents(sub_object) = "orig.ident"

test_genes <- read.csv("output/20200623/test_genes.csv",header = F,stringsAsFactors = F)
test_genes = test_genes$V1
test_genes = test_genes[test_genes %in% rownames(sub_object)]

if(s != "All_samples") sub_object <- subset(sub_object, idents = s)
## Column clustering (adjust here distance/linkage methods to what you need!)
sub_object <- FindVariableFeatures(object = sub_object, selection.method = "vst",
                               num.bin = 20,nfeatures = 3000,
                               mean.cutoff = c(0.1, 8),
                               dispersion.cutoff = c(1, Inf))
v_genes <- unique(c(VariableFeatures(sub_object),test_genes))
# remove MT-
MT <- grepl("^MT-",v_genes)
if(any(MT)) v_genes = v_genes[!MT]

y = sub_object[["SCT"]]@data[v_genes,]
system.time(cor_res <- Hmisc::rcorr(t(as.matrix(y)), type="spearman"))
cor_res$r[is.na(cor_res$r)] = 0

save.path <- paste0(path,file,"-Correlation-pvalues/", s,"/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

saveRDS(cor_res, file = paste0(path,file,"-Correlation-pvalues/",
                               args,"-",file,"-",s,".rds"))

if(heatmap){
        jpeg(paste0(save.path,"heatmap-cor-",file,"-",s,".jpeg"), units="in", width=10, height=7,res=600)
        heatmap.2(cor_res$r,trace="none")
        dev.off()
}

if(jpg){
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
                jpeg(paste0(save.path,"cor-pvalue-",gene,"-",s,".jpeg"), units="in", width=10, height=7,res=600)
                g <- ggline(df, x = "correlation", y = "log10.p.value",
                            numeric.x.axis = TRUE,
                            ylab = "-log10(p-value)",
                            xlab = "Spearman Correlation",
                            font.label = list(size = 25, face = "plain"),
                            label = "genes",             # Add point labels
                            label.select = c(head(df$genes,n),tail(df$genes,n)),
                            repel = TRUE,
                            title = paste("Genes correlated with",gene,"in",s,cell.type))+
                        theme_minimal() + TitleCenter()
                print(g)
                dev.off()
                colnames(df)[3] = " -log10(p-value)"
                df_list[[i]] = df
                Progress(i,length(test_genes))
        }
        names(df_list) = test_genes
        write.xlsx(df_list, file = paste0(save.path,"cor-pvalue-summary-",s,".xlsx"),
                   colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
}
