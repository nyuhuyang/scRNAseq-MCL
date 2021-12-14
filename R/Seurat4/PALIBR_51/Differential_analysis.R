invisible(lapply(c("magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

csv_list <- list.files(pattern = "FC0.01",path = "output/20210511",full.names = T)
deg_list <- pbapply::pblapply(csv_list, function(x){
    tmp = read.csv(x,row.names = 1)
    tmp = tmp[tmp$cluster == "KO", ]
    tmp$cluster %<>% paste0("/WT")
    tmp = tmp[abs(tmp$avg_log2FC) >0.1,]
    tmp = tmp[order(tmp$avg_log2FC,decreasing = T), ]
    #tmp = tmp[tmp$avg_log2FC > 0, ]
    tmp

})
names(deg_list) = gsub(".*FC0.01_","",csv_list) %>% sub("\\.csv","",.)
openxlsx::write.xlsx(deg_list, file =  paste0(path,"differential_expressed.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))
