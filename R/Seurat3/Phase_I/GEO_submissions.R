# Generate raw_counts_txt


library(magrittr)
# read fastq files
links = read.csv("doc/fastq_link.txt",stringsAsFactors = F) %>% .[,1] %>% 
        gsub("/*.fastq.gz","",., fixed=TRUE) %>%
        gsub(".*/Sample_","Sample_",.) 
        
links

# submit fastq to GEO
if(!dir.exists("data/geo_submission_fastq")) dir.create("data/geo_submission_fastq")
samples <- list.dirs("data/fastq")
samples <- samples[3:length(samples)]
for(s in samples[2:length(samples)]) {
        fastq_fils <- list.files(s)
        file.rename(from = paste0(s,"/",fastq_fils),
                    to = paste0("data/geo_submission_fastq/",fastq_fils))
}
for(s in samples[1:length(samples)]) unlink(s,recursive = T)

# write counts to GEO
library(Matrix)
library(coop)
library(data.table)
source("../R/Seurat3_functions.R")
(load(file="data/MCL_41_harmony_20191230.Rda"))

object_list <- SplitObject(object, split.by = "orig.ident")
counts_list <- lapply(object_list, function(x) x@assays$RNA@counts)
dt_list <- lapply(counts_list, function(x) as.data.table(x,keep.rownames=TRUE))
(samples <- sapply(counts_list, function(x) sub("_.*","",colnames(x)[1])))
#
for(i in seq_along(samples)){
        fwrite(dt_list[[i]],file = paste0("data/geo_submission_fastq/MCL_10X_raw_counts_",
               samples[i],".txt"),
               sep = "\t",row.names = FALSE, col.names = TRUE)
        Progress(i, length(samples))
}

#format(object.size(object_list),unit = "GB")
#format(object.size(counts_list),unit = "GB")
#format(object.size(dt_list),unit = "GB")
#format(object.size(DT),unit = "GB")

#' calculate sparsity of dgCMatrix. be awarid to convert integer to numeric for large dgCMatrix.
sparsity <- function(SMT){
        return(1 - nnzero(SMT)/(as.numeric(nrow(SMT)) * ncol(SMT)))
}

# linux
files=$(ls)
echo $files
md5sum $files > MCL_fastq_counts_checksums.md5
md5sum -c BladderCancer_bam_checksums.md5
"https://www.howtoforge.com/linux-md5sum-command/"

ncftp
set passive on
set so-bufsize 33554432
open ftp://geoftp:rebUzyi1@ftp-private.ncbi.nlm.nih.gov
cd uploads/yah2014_Ixm16U6w
mkdir MCL_data
cd MCL_data
put -R *
