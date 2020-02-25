
(csv_list = list.files("Yang/B_pairwise_heatmaps/DE_analysis_files"))
for(i in seq_along(csv_list)){
        gde.markers = read.csv(paste0(path,"DE_analysis_files/", csv_list[i]),
                               row.names = 1, stringsAsFactors = F)
        print(table(gde.markers$cluster1.vs.cluster2))
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        Top_n = 40
        top <-  gde.markers %>% group_by(cluster1.vs.cluster2) %>% 
                top_n(Top_n, avg_logFC) %>% as.data.frame()
        write.csv(top,file = paste0(path,"DE_analysis_files/top40_", csv_list[i]))
}
