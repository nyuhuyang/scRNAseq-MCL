p <- read.delim2("./output/MCL_expression.txt",row.names = 1)
len <- ncol(p)
for(i in 1:ncol(p)) 
        {
        d <-p[,i]
        l <- p[,len] #accessing the length column
        cS <- sum(as.numeric(p[,i])) #Total mapped reads per sample 
        rpkm[[i]] <- (10^9)*(as.numeric(p[,i]))/(as.numeric(l)*cS)
        rpkm[[1]] <- p[[1]]
        }

write.table(rpkm,"output_table_rpkm.txt",sep="\t",quote=F,row.names=F)