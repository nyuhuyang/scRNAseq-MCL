# 1) 	Install R packages from within R:
install.packages('Rserve') 
install.packages('e1071')  
install.packages('parallel')
install.packages('colorRamps')
source("https://bioconductor.org/biocLite.R")
biocLite("preprocessCore")

# 2) 	Start Rserve with option i or ii:
# i) Within R console:
library(Rserve)
Rserve(args="--no-save")
get()
##> SOCK_ERROR: bind error #48(address already in use)
# solve: restart the computer

# ii) From Unix shell:
path=/Library/Frameworks/R.framework/Versions/3.5/Resources/library/Rserve/libs/
R CMD $path/Rserve --no-save
cd R/CIBERSORT_package
java -Xmx3g -Xms3g -jar CIBERSORT.jar -M MCL_expression.txt -B LM22.txt
#>R version 3.5.1 (2018-07-02)
#Error: long vectors not supported yet: qap_encode.c:36
# Fatal error: unable to initialize the JIT
" didn't solve"
