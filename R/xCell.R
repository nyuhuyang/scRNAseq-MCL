# xCell
source("http://bioconductor.org/biocLite.R")
biocLite((c("GSVA","GSEABase")))
library(devtools)
install_github("dviraran/xCell")

library(xCell)
exprMatrix = read.table(expr_file,header=TRUE,row.names=1, as.is=TRUE)

xCellAnalysis(xCell.data)

# SingleR
# http://comphealth.ucsf.edu/sample-apps/SingleR/SingleR_specifications.html
devtools::install_github('dviraran/SingleR')

library(xCell)
source("https://bioconductor.org/biocLite.R")
biocLite("GSVAdata")
library(GSVAdata)
MCL.expression <- read.delim2("./output/MCL_expression.txt",row.names = 1)
LM22 <- read.delim2("./data/LM22.txt",row.names = 1)
LM22 <- data.matrix(LM22) # convert dataframe to matrix
head(sapply(LM22,class))

# convert dataframe column to list r
genes <- rownames(LM22)
LM22.list <- lapply(1:ncol(LM22), function(x) LM22[,x])
LM22.list <- setNames(LM22.list, colnames(LM22))
lapply(LM22.list,head)
LM22.list <- lapply(LM22.list, function(x) x[order(x,decreasing = T)])
lapply(LM22.list,head)
LM22.list <- lapply(LM22.list, function(x) names(x[1:300]))
lapply(LM22.list,length)
genes <-unique(unlist(LM22.list))
length(genes)
results_xCell <- xCellAnalysis.1(expr=MCL.expression, signatures = LM22.list)
