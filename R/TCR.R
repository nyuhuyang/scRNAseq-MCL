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
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data/")) dir.create("data")
########################################################################
#
#  4 TCR analysis
# 
# ######################################################################
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
head(meta.data_new,2)

# count cdr3===================
cdr3_table <- table(meta.data_new$cdr3) %>% as.data.frame
table(cdr3_table$Freq)
cdr3_table = cdr3_table[!(cdr3_table$Var1 %in% "None"),]
cdr3_table[order(cdr3_table$Freq,decreasing = T),] %>% .[1:20,] %>%
        kable %>% kable_styling

cdr3_table[which(cdr3_table$Freq == 143),]
meta.data_new[(meta.data_new$cdr3 %in% "CASSFLVKEQYF"),] %>% .[1:5,] %>%
        kable %>% kable_styling
cdr3_table[which(cdr3_table$Freq == 51),]
meta.data_new[(meta.data_new$cdr3 %in% "CASRQGLDTEAFF"),] %>% 
        kable %>% kable_styling
# generate clonotype_id===================
cdr3_table$log_freq = log10(cdr3_table$Freq)
colnames(cdr3_table)[1] = "cdr3"
meta.data_new2 <- left_join(meta.data_new, cdr3_table,by = "cdr3")
dim(meta.data);dim(meta.data_new2);dim(cdr3_table)


top <-  meta.data_new2 %>% 
        group_by(log_freq) %>% 
        top_n(2000, Freq) %>% .[,c("Barcode","orig.ident","singler1sub",
                                  "high_confidence","chain","v_gene","d_gene","j_gene","c_gene",
                                  "productive","reads","umis","cdr3","Freq","log_freq")]
write.csv(top,file = paste0(path,"TCR.csv"))

#====== 4.3 Pt =========================
Frequency <- function(df, col.name = "cdr3",top=100){
        cdr3_table <- table(df[,col.name]) %>% as.data.frame
        table(cdr3_table$Freq) %>% print
        Sum <- sum(cdr3_table$Freq)
        if(colnames(Pt_17_meta.data)[1] == "nGene"){
                cdr3_table$frequency = cdr3_table$Freq/Sum   
        }
        cdr3_table = cdr3_table[order(cdr3_table$Freq,decreasing = T),]
        cdr3_table = cdr3_table[cdr3_table$Var1!="na",]
        colnames(cdr3_table)[1] = "cdr3"
        df_new <- left_join(df, cdr3_table,by = "cdr3")
        df_new = df_new[,-which(colnames(df_new) %in% "Freq")]
        df_new = df_new[!duplicated(df_new$cdr3),]
        df_new = df_new[order(df_new$freq,decreasing = T),]
        top <- min(top,nrow(df_new))
        return(df_new[1:top,])
}


#' Generate Trends of top 20 TCR clones in longitudinal samples
#' @param TCR_data TCR data.frame report. must contain columns like c("orig.ident","cdr3","frequency")
#' @param order.by character vector indicates sample longitudinal order
#' @param group character, the summary of all samples. For output figure titile only
#' @param size ggplot point size
#' @param top integer, include top N TCR clones from each time points
#' @param x.margin integer, blank margin on x axis
#' @example order.by = c("Pt17_C03_TCRB","Pt17_C07_TCRB","Pt17_C31_TCRB")
#' TCR_Longitud_Plot(Pt_17_TCR,order.by, "Pt-17 bulk")
TCR_Longitud_Plot <- function(TCR_data, order.by, group, 
                              size=2,colors = singler.colors,x.margin=0.125,
                              do.return = FALSE, do.print = TRUE,top=20){
        colnames(TCR_data) = tolower(colnames(TCR_data))
        TCR_data$orig.ident = as.character(TCR_data$orig.ident)
        TCR_data <- split(TCR_data, TCR_data$orig.ident) %>% 
                lapply(function(x) Frequency(x,top=top)) %>% bind_rows
        keep_col <- c("orig.ident","cdr3","frequency")
        TCR_data = TCR_data[order(TCR_data$frequency,decreasing = T),]
        TCR_data <- TCR_data[,keep_col]
        # fill up 0
        TCR_data %<>% spread(key="orig.ident", value="frequency",fill = 0)
        TCR_data %<>% gather(key="orig.ident",value="frequency",-cdr3)

        TCR_data$orig.ident = factor(TCR_data$orig.ident,levels = order.by)
        TCR_data$samples <- as.integer(TCR_data$orig.ident)
        g1 <- ggplot(data=TCR_data, aes(x=samples, y=frequency, group=cdr3, color=cdr3)) +
                geom_line() + 
                geom_point(size = size)+
                scale_fill_manual(values = colors)+
                theme_minimal()+
                xlab("longitudinal samples")+
                scale_x_continuous(limits=c(min(TCR_data$samples)-x.margin,
                                            max(TCR_data$samples)+x.margin),
                                   breaks=(min(TCR_data$samples)):(max(TCR_data$samples)), 
                                   labels=order.by)+
                scale_y_continuous(labels = scales::percent)+
                ggtitle(paste("Trends of top",top,"TCR clones in longitudinal",group,
                              "samples"))+
                theme(text = element_text(size=15),			
                      legend.position="none",
                      plot.title = element_text(hjust = 0.5),
                      axis.text.x  = element_text(angle=70, vjust=0.5))
        if(do.print){
                path <- paste0("output/",gsub("-","",Sys.Date()),"/")
                if(!dir.exists(path)) dir.create(path, recursive = T)
                jpeg(paste0(path,group,"_TCR.jpeg"), units="in", width=10, height=7,res=600)
                print(g1)
                dev.off()
        }
        if(do.return) return(g1)
}


#' Generate Trends of top 20 TCR clones in longitudinal samples
#' @param TCR_data TCR data.frame report. must contain columns like c("orig.ident","cdr3","frequency")
#' @param order.by length =2character vector indicates pariwise sample at x-axis and y-axis
#' @param top integer, include top N TCR clones from each time points
#' @param size ggplot point size
#' @param margin integer, blank margin on both x and y axis
#' @example TCRPairPlot(Pt_17_TCR,c("Pt17_C03_TCRB","Pt17_C07_TCRB"))
TCRPairPlot <- function(TCR_data, order.by, size=2,margin=0.1,
                        do.return = FALSE, do.print = TRUE,top=20){
        colnames(TCR_data) = tolower(colnames(TCR_data))
        TCR_data = TCR_data[(TCR_data$orig.ident %in% order.by),]
        TCR_data <- split(TCR_data, TCR_data$orig.ident) %>% 
                lapply(function(x) Frequency(x,top = top)) %>% bind_rows
        keep_col <- c("orig.ident","cdr3","frequency")
        TCR_data = TCR_data[order(TCR_data$frequency,decreasing = T),]
        TCR_data <- TCR_data[,keep_col]
        # fill up 0
        TCR_data %<>% spread(key="orig.ident", value="frequency",fill = 0)
        TCR_data = TCR_data[,1:3]
        colnames(TCR_data)[2:3] =c("sample1","sample2")
        TCR_data$counts = as.numeric(TCR_data$sample1>0)+
                as.numeric(TCR_data$sample2>0)
        
        g1 <- ggplot(data=TCR_data, aes(x=sample1, y=sample2, 
                                             group=counts, color=counts)) +
                geom_point(size = 2)+
                scale_fill_manual(values = singler.colors)+
                theme_minimal()+
                xlab(paste("frequency in",order.by[1]))+
                ylab(paste("frequency in",order.by[2]))+
                scale_x_sqrt(labels = scales::percent,expand=c(0,0),
                             limits=c(0,max(TCR_data$sample1)*(1+margin)))+
                scale_y_sqrt(labels = scales::percent,expand=c(0,0),
                             limits=c(0,max(TCR_data$sample2)*(1+margin)))+
                ggtitle(paste("Top",top,"TCR clones between",
                              paste(order.by,collapse = " and ")))+
                theme(text = element_text(size=15),			
                      legend.position="none",
                      plot.title = element_text(hjust = 0.5),
                      axis.text.x  = element_text(angle=70, vjust=0.5))

        if(do.print){
                path <- paste0("output/",gsub("-","",Sys.Date()),"/")
                if(!dir.exists(path)) dir.create(path, recursive = T)
                jpeg(paste0(path,paste0(order.by,collapse = "_"),"_pairTCR.jpeg"), units="in", width=10, height=7,res=600)
                print(g1)
                dev.off()
        }
        if(do.return) return(g1)
}
# bulk Pt-17
Pt_17_TCR <- read.table("data/Pt-17-TCR.tsv",sep = "\t",header = T)
head(Pt_17_TCR)
change_name = colnames(Pt_17_TCR) %in% c("sample_name","amino_acid")
colnames(Pt_17_TCR)[change_name] = c("orig.ident","cdr3")
order.by = c("Pt17_C03_TCRB","Pt17_C07_TCRB","Pt17_C31_TCRB")
TCR_Longitud_Plot(Pt_17_TCR,order.by, "Pt-17 bulk",top=20)

TCRPairPlot(Pt_17_TCR,c("Pt17_C03_TCRB","Pt17_C07_TCRB"),top=1000)

# bulk Pt-25
Pt_25_TCR <- read.table("data/Pt-25-TCR.tsv",sep = "\t",header = T)
head(Pt_25_TCR,1)
change_name = colnames(Pt_25_TCR) %in% c("sample_name","amino_acid")
colnames(Pt_25_TCR)[change_name] = c("orig.ident","cdr3")
order.by = c("Pt25_C04_TCRB","Pt25_C20_TCRB","Pt25_C24_TCRB")
TCR_Longitud_Plot(Pt_25_TCR,order.by, "Pt-25 bulk",top=20)

TCRPairPlot(Pt_25_TCR,c("Pt25_C20_TCRB","Pt25_C24_TCRB"),top=1000)


# Formenti, Rudqvist, Nature Medicine 2018

Pt_17_meta.data = meta.data_new2[(meta.data_new2$groups %in% "Pt-17"),]
# Estimate frequency of elements in a specified column
order.by = c("Pt-17-LN-C1","Pt-17-C2","Pt-17-C7","Pt-17-C31")
TCR_Longitud_Plot(Pt_17_meta.data,order.by, "Pt-17")

Pt_25_meta.data = meta.data_new2[(meta.data_new2$groups %in% "Pt-25"),]
order.by = c("Pt-25-SB-C1","Pt-25-C1","Pt-25-C1D8","Pt-25-C24","Pt-25-AMB-C25","Pt-25-C25")
TCR_Longitud_Plot(Pt_25_meta.data,order.by, "Pt-25")

order.by = c("Pt-25-SB-C1","Pt-25-C1","Pt-25-C1D8","Pt-25-C24","Pt-25-C25")
Pt_25_meta.data = Pt_25_meta.data[(Pt_25_meta.data$orig.ident %in% order.by),]
TCR_Longitud_Plot(Pt_25_meta.data,order.by, "Pt-25")

# TSNE
T_cell <- SubsetData(object , cells.use = meta.data_new2$Barcode)
T_cell %<>% SetAllIdent("res.0.6")
TSNEPlot.1(T_cell,colors.use = ExtractMetaColor(T_cell),do.label = T)
T_cell <- SubsetData(T_cell,ident.remove = c(0,1,4,5,6))
remove_cell <- FeaturePlot(T_cell,features.plot = "CD3E",do.identify = T)
all_cell <- meta.data_new2$Barcode
cells.use <- all_cell[!(meta.data_new2$Barcode %in% remove_cell)]
T_cell <- SubsetData(T_cell,cells.use = cells.use)
T_cell %<>% SetAllIdent("singler1sub")
TSNEPlot.1(T_cell,colors.use = ExtractMetaColor(T_cell),do.label = T,do.print = T)


meta.data_new3 = meta.data_new2[(meta.data_new2$Barcode %in% T_cell@cell.names),]
T_cell;dim(meta.data_new2);dim(meta.data_new3)
rownames(meta.data_new3) = meta.data_new3$Barcode
meta.data_new3$Freq[is.na(meta.data_new3$Freq)]=0
meta.data_new3$log_freq[is.na(meta.data_new3$log_freq)]=0
T_cell@meta.data = meta.data_new3
SingleFeaturePlot.1(T_cell,feature = "Freq",threshold = 1.5, no.legend = F,
                legend.title = "total TCR clones",do.print=T)
SingleFeaturePlot.1(T_cell,feature = "log_freq",threshold = 0.001, no.legend = F,
                    legend.title = "total log10 TCR clones",do.print=T)
