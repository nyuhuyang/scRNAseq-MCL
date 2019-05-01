########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("R.utils","Seurat","dplyr","kableExtra","ggplot2"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
        }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("doc")) dir.create("doc")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
args <- commandArgs(trailingOnly = TRUE)
if(length(args)==0) stop("Sample info must be provided in excel!")
df_samples <- readxl::read_excel(args[1])
print(df_samples)
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:12)))
df_samples <- df_samples[sample_n,]
samples = sample
attach(df_samples)
#df_samples[sample_n,] %>% kable() %>% kable_styling()
table(tests);nrow(df_samples)
list_samples <- lapply(colnames(df_samples), function(col) df_samples[,col])
names(list_samples) = colnames(df_samples)
keep = sapply(list_samples, function(n) length(n[!is.na(n)])>1)
list_samples =list_samples[keep]

# check missing data
current <- list.files("data/scRNA-seq")
(current <- current[!grepl(".Rda|RData",current)])
(missing_data <- list_samples$sample.id[!(list_samples$sample.id %in% current)])

# select species
if(unique(list_samples$species) == "Homo_sapiens") species <- "hg19"
if(unique(list_samples$species) == "Mus_musculus") species <- "mm10"
if(unique(list_samples$species) == "Danio_rerio") species <- "danRer10"
if(species == "hg19") suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
if(species == "mm10") suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))

message("Copying the datasets")
if(length(missing_data)>0){
        # Move files from Download to ./data and rename them
        for(missing_dat in missing_data){
                old.pth  <- paste("~/Downloads", missing_dat,"outs",
                                  "filtered_feature_bc_matrix",sep = "/")
                list.of.files <- list.files(old.pth)
                new.folder <- paste("./data/scRNA-seq", missing_dat,"outs",
                                    "filtered_gene_bc_matrices",species,sep = "/")
                if(!dir.exists(new.folder)) dir.create(new.folder, recursive = T)
                # copy the files to the new folder
                file.copy(paste(old.pth, list.of.files, sep = "/"), new.folder)
                print(list_files <- list.files(new.folder))
        }
}

message("Loading the datasets")
## Load the dataset
Seurat_raw <- list()
Seurat_list <- list()
for(i in 1:length(list_samples$sample)){
        Seurat_raw[[i]] <- Read10X(data.dir = paste0("data/",list_samples$sample.id[i],
                                   "/outs/filtered_gene_bc_matrices/",species))
        colnames(Seurat_raw[[i]]) = paste0(list_samples$sample[i],"_",colnames(Seurat_raw[[i]]))
        Seurat_list[[i]] <- CreateSeuratObject(Seurat_raw[[i]],
                                                 min.cells = 0,
                                               min.features = 0)
        Seurat_list[[i]]@meta.data$tests <- list_samples$tests[i]

}
remove(Seurat_raw);GC()


#======1.1.2 QC before merge =========================
# if args 2 is passed
args[2] = as.character(args[2])
if(is.na(args[2])){
        message("Starting QC")
        cell.number <- sapply(Seurat_list, function(x) length(colnames(x)))
        QC_list <- lapply(Seurat_list, function(x) as.matrix(GetAssayData(x, slot = "counts")))
        median.nUMI <- sapply(QC_list, function(x) median(colSums(x)))
        median.nGene <- sapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0])))))
        
        min.nUMI <- sapply(QC_list, function(x) min(colSums(x)))
        min.nGene <- sapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0])))))
        
        QC.list <- cbind(df_samples,cell.number, median.nUMI, median.nGene, 
                         min.nUMI,min.nGene, row.names = samples)
        write.csv(QC.list,paste0(path,"QC_list.csv"))
        #QC.list %>% kable() %>% kable_styling()
        remove(QC_list,median.nUMI,median.nGene,min.nUMI,min.nGene,QC.list);GC()
}
#========1.1.3 merge ===================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), Seurat_list)
remove(Seurat_list);GC()

# read and select mitochondial genes
if(species == "hg19") mito = "^MT-"
if(species == "mm10") mito = "^mt-" # not Mt-
if(species == "danRer10") mito = "^mt-"
message("mito.genes:")

(mito.features <- grep(pattern = mito, x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = mito)
Idents(object) = factor(Idents(object),levels = samples)
g1 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=8),legend.position="none")
})
save(g1,file= paste0(path,"g1","_",length(sample_n),"_",gsub("-","",Sys.Date()),".Rda"))
