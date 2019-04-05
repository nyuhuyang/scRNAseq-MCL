########################################################################
#
#  0 setup environment, install libraries if nLynchessary, load libraries
# 
# ######################################################################

library(Seurat)
library(magrittr)
library(dplyr)
library(tidyr)
library(kableExtra)
library(RColorBrewer)
library(SingleR)
library(ggplot2)
source("../R/Seurat_functions.R")
source("R/util.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data/")) dir.create("data")
########################################################################
#
#  4 TCR analysis
# 
# ######################################################################
#====== 4.0 Pt =========================
# bulk Pt-17
Pt_17_TCR <- read.table("data/Pt-17-TCR.tsv",sep = "\t",header = T)
head(Pt_17_TCR)
change_name = colnames(Pt_17_TCR) %in% c("sample_name","amino_acid")
colnames(Pt_17_TCR)[change_name] = c("orig.ident","cdr3")
order.by = c("Pt17_C03_TCRB","Pt17_C07_TCRB","Pt17_C31_TCRB")
TCR_Longitud_Plot(Pt_17_TCR,order.by, "Pt-17 bulk",do.log = T,top=300,
                  size =0.2,color = "dodgerblue4", do.print = T)

TCRPairPlot(Pt_17_TCR,c("Pt17_C03_TCRB","Pt17_C07_TCRB"),top=300)
TCRPairPlot(Pt_17_TCR,c("Pt17_C07_TCRB","Pt17_C31_TCRB"),top=1000)

# bulk Pt-25
Pt_25_TCR <- read.table("data/Pt-25-TCR.tsv",sep = "\t",header = T)
head(Pt_25_TCR,1)
change_name = colnames(Pt_25_TCR) %in% c("sample_name","amino_acid")
colnames(Pt_25_TCR)[change_name] = c("orig.ident","cdr3")
order.by = c("Pt25_C04_TCRB","Pt25_C20_TCRB","Pt25_C24_TCRB")
TCR_Longitud_Plot(Pt_25_TCR,order.by, "Pt-25 bulk",top=300,
                  size =0.2,color = "dodgerblue4", do.print = F)

TCRPairPlot(Pt_25_TCR,c("Pt25_C20_TCRB","Pt25_C24_TCRB"),top=1000)


#====== 4.1 copy the data files and read.csv =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/190320_scRNAseq_info.xlsx",sheet = "TCR")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:8)))
df_samples <- df_samples[sample_n,]
attach(df_samples)
(samples = sample)
# check missing data
(current <- list.files("data")[!grepl(".Rda|RData",list.files("data"))])
(missing_data <- sample.id[!(sample.id %in% current)])

if(length(missing_data)>0){
        # Move files from Download to ./data and rename them
        for(missing_dat in missing_data){
                old.pth  <- paste("~/Downloads", missing_dat,"outs",sep = "/")
                list.of.files <- list.files(old.pth)
                contig <- c("filtered_contig_annotations.csv",
                            "all_contig_annotations.csv")
                new.folder <- paste("./data", missing_dat,"outs",sep = "/")
                if(!dir.exists(new.folder) & all(contig %in% list.of.files)) {
                    dir.create(new.folder, recursive = T)
                    file.copy(paste(old.pth, contig, sep = "/"), new.folder)
                list.files(new.folder) %>% head
                }
        }
}

# read csv
contig_list <- mapply(function(x,y) {
        temp_dat <- read.csv(paste0("data/", x,"/outs/filtered_contig_annotations.csv"))
        temp_dat[,"barcode"] = paste0(y,"_", temp_dat[,"barcode"])
        return(as.data.frame(temp_dat))
        }, sample.id, samples, SIMPLIFY = FALSE)
contigs <- dplyr::bind_rows(contig_list)
table(contigs$is_cell)
table(contigs$high_confidence)
contigs = contigs[(contigs$high_confidence == "True"),]

# Prepare contigs
colnames(contigs)[colnames(contigs) == "barcode"] = "Barcode"
Barcode_table <- table(contigs$Barcode) %>% as.data.frame

table(Barcode_table$Freq)
Barcode_table[which(Barcode_table$Freq == 10),]
contigs[(contigs$Barcode == "Pt-11-C31_GATGCTAGTAGTGAAT-1"),] %>% 
        kable %>% kable_styling
contigs$cdr3 %in% "None" %>% table

contigs_cdr3 = contigs[!(contigs$cdr3 %in% "None"),]
table(contigs_cdr3$productive)
contigs_cdr3 = contigs_cdr3[(contigs_cdr3$productive %in% "True"),]
contigs_cdr3$chain %>% table

contigs_cdr3 = contigs_cdr3[(contigs_cdr3$chain %in% "TRB"),]
contigs_cdr3 = contigs_cdr3[order(contigs_cdr3$umis,decreasing = T),]

contigs_cdr3 = contigs_cdr3[!duplicated(contigs_cdr3$Barcode),]
dim(contigs_cdr3)
#======4.3 prepare meta.data =========================
#meta.data
(load(file="data/MCL_Harmony_30_20190320.Rda"))
meta.data <- object@meta.data
meta.data = meta.data[grep("T_cells",meta.data$singler1sub),]
dim(meta.data)
meta.data$Barcode = rownames(meta.data)

# merge contigs
meta.data_new <- left_join(meta.data, contigs_cdr3,by = "Barcode")
rownames(meta.data_new) = meta.data_new$Barcode
dim(meta.data); dim(meta.data_new);dim(contigs_cdr3)
table(duplicated(meta.data_new$Barcode))
head(meta.data_new,1)

# count cdr3===================
cdr3_table <- table(meta.data_new$cdr3) %>% as.data.frame
table(cdr3_table$Freq)
cdr3_table = cdr3_table[!(cdr3_table$Var1 %in% "None"),]

cdr3_table[order(cdr3_table$Freq,decreasing = T),] %>% .[1:20,] %>%
        kable %>% kable_styling
# generate clonotype_id===================
Sum = sum(cdr3_table$Freq)
cdr3_table$frequency = cdr3_table$Freq/Sum
colnames(cdr3_table)[1] = "cdr3"
meta.data_new2 <- left_join(meta.data_new, cdr3_table,by = "cdr3")
dim(meta.data);dim(meta.data_new2);dim(cdr3_table)
top <-  meta.data_new2[,c("Barcode","orig.ident","singler1sub",
                                  "high_confidence","chain","v_gene","d_gene","j_gene","c_gene",
                                  "productive","reads","umis","cdr3","Freq","frequency")]
top <- top[complete.cases(top),]
write.csv(top,file = paste0(path,"TCR.csv"))

# Formenti, Rudqvist, Nature Medicine 2018
Pt_11_meta.data = meta.data_new2[(meta.data_new2$groups %in% "Pt-11"),]
order.by = c("Pt-11-LN-C1", "Pt-11-C1","Pt-11-C14","Pt-11-C28")
TCR_Longitud_Plot(Pt_11_meta.data,order.by, group="Pt-11",top=50,do.log = T,
                  size =1,color = "dodgerblue4", do.print = T)

Pt_17_meta.data = meta.data_new2[(meta.data_new2$groups %in% "Pt-17"),]
order.by = c("Pt-17-LN-C1","Pt-17-C2","Pt-17-C7","Pt-17-C31")
TCR_Longitud_Plot(Pt_17_meta.data,order.by, group="Pt-17",top=50,do.log = T,
                  size =1,color = "dodgerblue4", do.print = T)

Pt_25_meta.data = meta.data_new2[(meta.data_new2$groups %in% "Pt-25"),]
order.by = c("Pt-25-SB-C1","Pt-25-C1","Pt-25-C1D8","Pt-25-C24","Pt-25-AMB-C25","Pt-25-C25")
order.by = c("Pt-25-SB-C1","Pt-25-C1","Pt-25-C1D8","Pt-25-C24","Pt-25-C25")

Pt_25_meta.data = Pt_25_meta.data[(Pt_25_meta.data$orig.ident %in% order.by),]
TCR_Longitud_Plot(Pt_25_meta.data,order.by, "Pt-25",top=30,do.log = T,
                  size =1,color = "dodgerblue4", do.print = T)

Pt_AA13_meta.data = meta.data_new2[(meta.data_new2$groups %in% "Pt-AA13"),]
order.by = c("Pt-AA13-Ib-p","Pt-AA13-Ib-1","Pt-AA13-Ib-R")
TCR_Longitud_Plot(Pt_AA13_meta.data,order.by, group="Pt_AA13",top=20,do.log = T,
                  size =1,color = "dodgerblue4", do.print = T)

AFT_03_meta.data = meta.data_new2[(meta.data_new2$groups %in% c("AFT-03")),]
order.by = c("AFT-03-C1D1","AFT-03-C1D8")
TCR_Longitud_Plot(AFT_03_meta.data,order.by, group="AFT-03",top=30,do.log = T,
                  size =1,color = "dodgerblue4", do.print = T)

AFT_04_meta.data = meta.data_new2[(meta.data_new2$groups %in% c("AFT-04")),]
order.by = c("AFT-04-LN-C1D1","AFT-04-C1D1","AFT-04-C1D8")
TCR_Longitud_Plot(AFT_04_meta.data,order.by, group="AFT-04",top=30,do.log = T,
                  size =1,color = "dodgerblue4", do.print = T)

Untreated_meta.data = meta.data_new2[(meta.data_new2$groups %in% c("Untreated")),]
order.by = c("Pt-RM","Pt-MS","Pt-LM","Pt-1475")
TCR_Longitud_Plot(Untreated_meta.data,order.by, group="Untreated",top=30,do.log = T,
                  size =1,color = "dodgerblue4", do.print = T)

# TSNE========
T_cell <- SubsetData(object , cells.use = meta.data_new2$Barcode)
T_cell %<>% SetAllIdent("res.0.6")
TSNEPlot.1(T_cell,colors.use = ExtractMetaColor(T_cell),do.label = T)
T_cell <- SubsetData(T_cell,ident.remove = c(0,1,4,5,6))
remove_cell <- FeaturePlot(T_cell,features.plot = "CD3E",do.identify = T)
remove_cell <- c("Pt-17-C7_AGTAGTCTCAAGATCC-1","Pt-25-C1D8_TCTGAGATCCATGCTC-1")
all_cell <- meta.data_new2$Barcode
cells.use <- all_cell[!(meta.data_new2$Barcode %in% remove_cell)]
T_cell <- SubsetData(T_cell,cells.use = cells.use)
T_cell %<>% SetAllIdent("singler1sub")
TSNEPlot.1(T_cell,colors.use = ExtractMetaColor(T_cell),do.label = T,do.print = T)

# prepare Persistent frequency and Enriched frequency
meta.data_new3 = meta.data_new2[(meta.data_new2$Barcode %in% T_cell@cell.names),]
T_cell;dim(meta.data_new2);dim(meta.data_new3)
rownames(meta.data_new3) = meta.data_new3$Barcode
meta.data_new3$Freq[is.na(meta.data_new3$Freq)]=0
meta.data_new3$frequency[is.na(meta.data_new3$frequency)]=0
meta.data_new3$TCR_type = NA
meta.data_new3[meta.data_new3$Freq == 1,"TCR_type"] = "Singlets"
meta.data_new3[meta.data_new3$Freq > 1,"TCR_type"] = "Multiplets"
table(meta.data_new3$TCR_type)
colnames(meta.data_new3)[colnames(meta.data_new3) %in% "Freq"] = "total_TCR"
colnames(meta.data_new3)[colnames(meta.data_new3) %in% "frequency"] = "total_frequency"

meta.data_new4 <- split(meta.data_new3, meta.data_new3$groups) %>%
        lapply(function(x) Frequency(x,remove.na =F,remove.dup =F)) %>% bind_rows
dim(meta.data_new3);dim(meta.data_new4)
colnames(meta.data_new4)[colnames(meta.data_new4) %in% "Freq"] = "Patient_specific_TCR"
colnames(meta.data_new4)[colnames(meta.data_new4) %in% "frequency"] = "Patient_specific_frequency"

common <- meta.data_new4$total_TCR >meta.data_new4$Patient_specific_TCR
table(common)
meta.data_new4[common,"total_TCR"] = "Common"
meta.data_new4$Patient_specific = NA
meta.data_new4[meta.data_new4$total_TCR >0,"Patient_specific"] = TRUE
meta.data_new4[common,"Patient_specific"] =FALSE
table(meta.data_new4$Patient_specific)

meta.data_new5 <- split(meta.data_new4, meta.data_new4$orig.ident) %>%
        lapply(function(x) Frequency(x,remove.na =F,remove.dup =F)) %>% bind_rows
dim(T_cell@meta.data);dim(meta.data_new3);dim(meta.data_new4);dim(meta.data_new5)
colnames(meta.data_new5)[colnames(meta.data_new5) %in% "Freq"] = "Enriched_TCR"
colnames(meta.data_new5)[colnames(meta.data_new5) %in% "frequency"] = "Enriched_frequency"
for(i in 40:ncol(meta.data_new5)){
        meta.data_new5[,i] = as.numeric(meta.data_new5[,i])
        meta.data_new5[is.na(meta.data_new5[,i]),i]=0
}
rownames(meta.data_new5) = meta.data_new5$Barcode
meta.data_new5 =meta.data_new5[T_cell@cell.names,]
T_cell@meta.data = meta.data_new5

#====== Generate TCR plots ========
#' @param key chose between "counts" or "frequency"
TCR_Split_TSNEPlot<- function(object, split.by="orig.ident",feature = "persistent_freqency",
                              order.by,threshold =NULL,nrow =2,no.legend=T,
                              do.return=TRUE, do.print= FALSE){
        if(is.null(threshold)) {
                g <-hist(object@meta.data[,feature])
                threshold =g$mids[1]
        }
        if(grepl("frequency",feature)) {
                legend_title <- "freqency"
                } else legend_title <- "clones"
        
        g <- lapply(order.by,function(sample) {
                SubsetData(object, ident.use = sample) %>%
                        SingleFeaturePlot.1(threshold= threshold,
                                            feature = feature,
                                            title = paste(sample),
                                            no.legend = no.legend,
                                            legend.title =legend_title)
        })
        g1 <- do.call(plot_grid, c(g, nrow = nrow))+
                      ggtitle(feature)+
                      theme(text = element_text(size=15),
                            plot.title = element_text(hjust = 0.5),
                            axis.text.x  = element_text(angle=70, vjust=0.5))
        if(do.print){
                path <- paste0("output/",gsub("-","",Sys.Date()),"/")
                if(!dir.exists(path)) dir.create(path, recursive = T)
                jpeg(paste0(path,deparse(substitute(object)),"_",feature,".jpeg"), units="in", 
                     width=10, height=7,res=600)
                print(g1)
                dev.off()
        }
        if(do.return) return(g1)

}


FeaturePlot(T_cell,"total_TCR")
SingleFeaturePlot.1(T_cell,feature = "total_TCR",threshold = 1.1, no.legend = F,
                    legend.title = "total TCR",do.print=T, do.return = T)
SingleFeaturePlot.1(T_cell,feature = "total_frequency",threshold = 0.0001, no.legend = F,
                    legend.title = "total TCR frequency",do.print=T, do.return = T)

T_cell %<>% SetAllIdent(id="orig.ident")

order.by = c("Pt-17-LN-C1","Pt-17-C2","Pt-17-C7","Pt-17-C31")
T_cell_Pt_17 = SubsetData(T_cell,ident.use = order.by)
TCR_Split_TSNEPlot(object=T_cell_Pt_17,split.by="orig.ident",
                   feature = "Enriched_frequency",
                   order.by = order.by,no.legend=F,threshold =0.005,
                   do.return = T,do.print = T)

order.by = c("Pt-25-SB-C1","Pt-25-C1","Pt-25-C1D8","Pt-25-C24","Pt-25-C25")
T_cell_Pt_25 = SubsetData(T_cell,ident.use = order.by)
TCR_Split_TSNEPlot(T_cell_Pt_25,split.by="orig.ident",
                   feature = "Patient_specific_frequency",
                   order.by = order.by,no.legend=T,threshold =NULL,
                   do.return = F,do.print = T)

order.by = c("Pt-RM","Pt-MS","Pt-LM","Pt-1475")
T_cell_Untreated = SubsetData(T_cell,ident.use = order.by)
TCR_Split_TSNEPlot(T_cell_Untreated,split.by="orig.ident",
                   feature = "Patient_specific_frequency",
                   order.by = order.by,no.legend=T,threshold =NULL,
                   do.return = T,do.print = F)









# label Persistent TCR---skip-----
cdr3_samples <- table(meta.data_new4$cdr3, meta.data_new4$orig.ident) %>% 
        as.data.frame %>% spread(Var1,Freq) 
rownames(cdr3_samples) = cdr3_samples$Var2
cdr3_samples = cdr3_samples[,-1] %>% t %>% as.data.frame

(groups <- unique(meta.data_new4$groups))
for(i in 1:length(groups)) {
        (sample_in_group <- meta.data_new4$orig.ident[meta.data_new4$groups %in%
                                                              groups[i]] %>% unique)
        cdr3_samples[,paste0(groups[i]],"_pos")] = apply(cdr3_samples[,sample_in_group], 1,
                                         function(x) sum(length(x[x>0])))
        cdr3_samples[,paste0(groups[i]],"counts")] = rowSums(cdr3_samples[,sample_in_group])
}

n = length(unique(meta.data_new4$orig.ident))
cdr3_samples = cdr3_samples[,(n+1):ncol(cdr3_samples)]
cdr3_samples_Common = rowSums(cdr3_samples)

