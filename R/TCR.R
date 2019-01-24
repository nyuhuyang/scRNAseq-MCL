########################################################################
#
#  0 setup environment, install libraries if nLynchessary, load libraries
# 
# ######################################################################

library(Seurat)
library(magrittr)
library(dplyr)
library(kableExtra)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("./data/")) dir.create("data")
########################################################################
#
#  4 TCR analysis
# 
# ######################################################################
#======4.1 copy the data files and read.csv =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/181227_scRNAseq_info.xlsx",sheet = "TCR")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:7)))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
samples <- df_samples$sample[sample_n]
sample.id <- df_samples$sample.id[sample_n]

# check missing data
(current <- list.files("data")[!grepl(".Rda|RData",list.files("data"))])
(missing_data <- sample.id[!(sample.id %in% current)])

if(length(missing_data)>0){
        # Move files from Download to ./data and rename them
        for(missing_dat in missing_data){
                old.pth  <- paste("~/Downloads", missing_dat,"outs",sep = "/")
                list.of.files <- list.files(old.pth)
                new.folder <- paste("./data", missing_dat,"outs",sep = "/")
                if(!dir.exists(new.folder)) dir.create(new.folder, recursive = T)
                # copy the files to the new folder
                file.copy(paste(old.pth, list.of.files, sep = "/"), new.folder)
                list.files(new.folder) %>% head
        }
}

# read csv
contig_list <- mapply(function(x,y) {
        temp_dat <- read.csv(paste0("./data/", x,"/outs/all_contig_annotations.csv"))
        temp_dat[,"barcode"] = paste0(y,"_", temp_dat[,"barcode"])
        return(as.data.frame(temp_dat))
        }, sample.id, samples, SIMPLIFY = FALSE)
contigs <- bind_rows(contig_list)
table(contigs$high_confidence)
#contigs = contigs[(contigs$high_confidence == "True"),]

# Prepare contigs
colnames(contigs)[colnames(contigs) == "barcode"] = "Barcode"
Barcode_table <- table(contigs$Barcode) %>% as.data.frame
table(Barcode_table$Freq)
Barcode_table[which(Barcode_table$Freq == 40),]
contigs[(contigs$Barcode == "Pt-LM_TGCACCTCAACTGGCC-1"),] %>% 
        kable %>% kable_styling
contigs$cdr3 %in% "None" %>% table
contigs_cdr3 = contigs[!(contigs$cdr3 %in% "None"),]
contigs_cdr3$chain %>% table

contigs_cdr3 = contigs_cdr3[(contigs_cdr3$chain %in% "TRB"),]
contigs_cdr3 = contigs_cdr3[!duplicated(contigs_cdr3$Barcode),]
contigs_cdr3 = rbind.data.frame(contigs_cdr3,contigs)
contigs_cdr3 = contigs_cdr3[!duplicated(contigs_cdr3$Barcode),]

#======4.3 load  SingleCellExperiment =========================
#meta.data
(load(file="data/MCL_Harmony_20_20181231.Rda"))
meta.data <- object@meta.data
meta.data$Barcode = rownames(meta.data)
(remove <- which(colnames(meta.data) %in%c("Sample","is_cell_control",
                                           "pct_counts_in_top_500_features_Mito")))
meta.data = meta.data[,-c(3,5:40)]

# merge contigs
meta.data_new <- left_join(meta.data, contigs_cdr3,by = "Barcode")
rownames(meta.data_new) = meta.data_new$Barcode
dim(meta.data); dim(meta.data_new);dim(contigs_cdr3)
table(meta.data_new$Barcode) %>% table
table(duplicated(meta.data_new$Barcode))
which(duplicated(meta.data_new$Barcode)) %>% head
head(meta.data_new,2)

raw_clonotype_table <- table(meta.data_new$raw_clonotype_id) %>% as.data.frame
table(raw_clonotype_table$Freq)
raw_clonotype_table[which(raw_clonotype_table$Freq == 14),]
meta.data_new[(meta.data_new$raw_clonotype_id %in% "clonotype50"),] %>% 
        kable %>% kable_styling
# count cdr3===================
cdr3_table <- table(meta.data_new$cdr3) %>% as.data.frame
table(cdr3_table$Freq)
cdr3_table = cdr3_table[!(cdr3_table$Var1 %in% "None"),]
cdr3_table[order(cdr3_table$Freq,decreasing = T),] %>% .[1:20,] %>%
        kable %>% kable_styling

cdr3_table[which(cdr3_table$Freq == 5613),]
meta.data_new[(meta.data_new$cdr3 %in% "CASSKGQETQYF"),] %>% .[1:5,] %>%
        kable %>% kable_styling
cdr3_table[which(cdr3_table$Freq == 120),]
meta.data_new[(meta.data_new$cdr3 %in% "CRVSGL*AKNIQYF"),] %>% 
        kable %>% kable_styling
# generate clonotype_id===================
cdr3_table$clonotype_id = NA
cdr3_table = cdr3_table[order(cdr3_table$Freq,decreasing = F),]
cdr3_table[(cdr3_table$Freq %in% 1),"clonotype_id"] = "clonotype1"
(N = which((cdr3_table$Freq >1))) %>% head
cdr3_table[N,"clonotype_id"] = paste0("clonotype",2:(length(N)+1))
tail(cdr3_table)
colnames(cdr3_table)[1] = "cdr3"
meta.data_new2 <- left_join(meta.data_new, cdr3_table,by = "cdr3")
dim(meta.data_new2);dim(meta.data_new);dim(cdr3_table)
rownames(meta.data_new2) = meta.data_new2$Barcode
meta.data_new2 = meta.data_new2[object@cell.names,]

top <-  meta.data_new2 %>% 
        group_by(clonotype_id) %>% 
        top_n(300, Freq) %>% .[,c("Barcode","orig.ident","res.0.6","singler1sub",
                                  "chain","v_gene","d_gene","j_gene","c_gene","productive",
                                  "cdr3","Freq","clonotype_id")]
write.csv(top,file = paste0(path,"TCR.csv"))
#======4.4 TSNEplot =========================
object@meta.data = meta.data_new2

object <- SetAllIdent(object, id = "clonotype_id")
g <- TSNEPlot.1(object = object, do.label = F, group.by = "ident",
                        do.return = TRUE, no.legend = F, 
                        #colors.use = ExtractMetaColor(object),
                        pt.size = 1,label.size = 6 )+
        ggtitle("Tsne plot of all clonotype")+
        theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"/TSNEplot_clonotype_L.jpeg"), units="in", width=10, height=7,res=600)
print(g)
dev.off()





table(object@meta.data$orig.ident)
table(object@ident)
object@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",object@meta.data$orig.ident)

df_samples <- readxl::read_excel("doc/181227_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
object <- SetAllIdent(object, id="singler1sub")
tests <- paste0("test",c(7))
for(test in tests){
        sample_n = which(df_samples$tests %in% c("control",test))
        samples <- unique(df_samples$sample[sample_n])
        
        cell.use <- rownames(object@meta.data)[object@meta.data$orig.ident %in% 
                                                       c("Normal",samples)]
        subset.object <- SubsetData(object, cells.use = cell.use)
        subset.object@meta.data$orig.ident %>% unique %>% sort %>% print
        g <- SplitTSNEPlot(subset.object,group.by = "ident",split.by = "orig.ident",
                           select.plots = c(1,3,2,4),#c(6:8,1:5)
                           no.legend = T,do.label =F,label.size=3,size=20,
                           return.plots =T, label.repel = T,force=2)
        jpeg(paste0(path,test,"_TSNEPlot.jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, c(g, nrow = 2)))
        dev.off()
}