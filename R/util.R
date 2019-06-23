#remove NA columns and NA rows, remove duplicate Gene_name
CleanUp <- function(df){
    
        rm_NA_col <- df[which(df[,1] == "Gene_name"),] %>% is.na %>% as.vector
        df = df[,!rm_NA_col]
        rm_NA_row <- apply(df,1, function(x) all(is.na(x)))
        df = df[!rm_NA_row,]
        colnames(df) = df[which(df[,1] == "Gene_name"),] %>% as.character
        df = df[-which(df[,1] == "Gene_name"),]
        
        rm_col <- colnames(df) %in% c("Gene_id","biotype")
        df = df[,!rm_col]
        df = RemoveDup(df)
        
        return(df)
}

# re-calculate the cdr3 frequency in data.frame
Frequency <- function(df, key = "cdr3",top=NULL,remove.na =T,
                      remove.dup =T, verbose = FALSE){
        
        #colnames(df) = tolower(colnames(df))
        if(any(colnames(df) %in% "frequency")){
                df = df[,-which(colnames(df) %in% "frequency")]
        }
        cdr3_table <- table(df[,key]) %>% as.data.frame
        if(verbose) table(cdr3_table$Freq) %>% print
        Sum <- sum(cdr3_table$Freq)
        if(colnames(df)[1] == "ngene"|colnames(df)[1] == "nGene"){
        cdr3_table$frequency = cdr3_table$Freq/Sum
        }
        cdr3_table = cdr3_table[order(cdr3_table$Freq,decreasing = T),]
        if(remove.na) cdr3_table = cdr3_table[cdr3_table$Var1!="na",]
        colnames(cdr3_table)[1] = "cdr3"
        df_new <- left_join(df, cdr3_table,by = "cdr3")
        if(remove.dup) df_new = df_new[!duplicated(df_new$cdr3),]
        df_new$Freq[is.na(df_new$Freq)]=0
        df_new$frequency[is.na(df_new$frequency)]=0
        df_new = df_new[order(df_new$freq,decreasing = T),]
        if(is.null(top)) top = nrow(df_new)
        top = min(top,nrow(df_new))
        return(df_new[1:top,])
}


#remove duplicate rownames with lower rowsumns
#' @param mat input as data.frame with gene name
#' @export mat matrix with gene as rownames, no duplicated genes
RemoveDup <- function(mat){
        gene_id <- as.matrix(mat[,1])
        mat <- mat[,-1]
        if(!is.matrix(mat)) mat <- sapply(mat,as.numeric)
        rownames(mat) <- 1:nrow(mat)
        mat[is.na(mat)] = 0
        mat <- cbind(mat, "rowSums" = rowSums(mat))
        mat <- mat[order(mat[,"rowSums"],decreasing = T),]
        gene_id <- gene_id[as.numeric(rownames(mat))]
        remove_index <- duplicated(gene_id)
        mat <- mat[!remove_index,]
        rownames(mat) <- gene_id[!remove_index]
        return(mat[,-ncol(mat)])
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
TCR_Longitud_Plot <- function(TCR_data, order.by, group,color=NULL,
                              size=2,colors = singler.colors,x.margin=0.125,
                              do.return = TRUE, do.print = FALSE,do.log =TRUE, top=20){
        
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
        g1 <- ggplot(data=TCR_data, aes(x=samples, y=frequency, group=cdr3,colour=cdr3)) +
                geom_line(size = size) +
                geom_point(size = size)+
                #scale_fill_manual(values = colors)+
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
        if(do.log) g1 = g1 + scale_y_log10()
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
                        do.return = TRUE, do.print = FALSE,top=20){
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


#' select 1/4 of cell from control
ScaleDown <- function(object, control=c("BH","DJ","MD","NZ")){
    
        normal_cells = lapply(control, function(x){
        rownames(object@meta.data)[(object@meta.data$orig.ident %in% x)]
        })
        set.seed(101)
        remove_normal_cells = lapply(normal_cells, function(x) sample(x, size = length(x)*3/4)) %>% unlist
        table(object@cell.names %in% remove_normal_cells)
        cell.use <- object@cell.names[!(object@cell.names %in% remove_normal_cells)]
        object <- SubsetData(object, cells.use = cell.use)
        object@meta.data$orig.ident = gsub(paste(control,collapse = "|"),"Normal",object@meta.data$orig.ident)
        
        return(object)
}

#' select 10% of cell from control
#' Seurat 3
ScaleDown.1 <- function(object, control=c("BH","DJ","MD","NZ"), scale.down=0.1){
    
    normal_cells = lapply(control, function(x){
        rownames(object@meta.data)[(object@meta.data$orig.ident %in% x)]
    })
    set.seed(101)
    remove_normal_cells = lapply(normal_cells, function(x) {
        sample(x, size = length(x)*(1-scale.down))
        }) %>% unlist
    table(colnames(object) %in% remove_normal_cells)
    cell.use <- colnames(object)[!(colnames(object) %in% remove_normal_cells)]
    object <- object[, cell.use]
    object@meta.data$orig.ident = gsub(paste(control,collapse = "|"),"Normal",object@meta.data$orig.ident)
    
    return(object)
}
