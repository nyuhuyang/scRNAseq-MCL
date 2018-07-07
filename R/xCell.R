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
