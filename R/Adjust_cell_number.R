#######################################
# Normalize cell number
######################################
library(gplots)
library(tidyr)
library(Seurat)
library(kableExtra)
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("./data/")) dir.create("data")

(load(file="data/MCL_Harmony_24_20190128.Rda"))

CBC <- readxl::read_excel("doc/190120 MCL CBC summary.xlsx",sheet = "Sheet2")
CBC %>% kable %>% kable_styling()
CBC <- readxl::read_excel("doc/190120 MCL CBC summary.xlsx") %>% as.data.frame()
CBC %>% kable %>% kable_styling()
sapply(CBC,class)
range(CBC$Monocytes)
median(CBC$Monocytes)
###########################
#==== CBC barchat=======
###########################
y =  CBC[,c("sample","Lymphocytes","Eosinophils","Neutrophils","Basophils","Monocytes")]

df <- gather(y, key = "Cell.Type", value = "counts", -sample)
#df = melt(y,variable.name = "Cell.Type", value.name="counts")
p <-ggplot(df, aes(sample, counts))+
        geom_bar(stat = "identity", aes(fill = Cell.Type)) +
        ylab("million counts/ml")+
        ggtitle("CBC with differential summary")+
        theme(text = element_text(size=15),     							
              plot.title = element_text(size=20,hjust = 0.5),
              axis.text.x = element_text(angle = 90, hjust = 1))+
        scale_fill_manual(values = c("#00B0F6","#A3A500","#00BF7D","#F8766D","#E76BF3"))
jpeg(paste0(path,"/barchart_CBC.jpeg"), units="in", width=10, height=7,res=600)
p
dev.off()

#-------CBC barchat no Lymphocytes----------
y =  CBC[,c("sample","Eosinophils","Neutrophils","Basophils","Monocytes")]

df <- gather(y, key = "Cell.Type", value = "counts", -sample)
#df = melt(y,variable.name = "Cell.Type", value.name="counts")
p <-ggplot(df, aes(sample, counts))+
        geom_bar(stat = "identity", aes(fill = Cell.Type)) +
        ylab("million counts/ml")+
        ggtitle("CBC with differential summary without Lymphocytes")+
        theme(text = element_text(size=15),     							
              plot.title = element_text(size=20,hjust = 0.5),
              axis.text.x = element_text(angle = 90, hjust = 1))+
        scale_fill_manual(values = c("#00B0F6","#A3A500","#F8766D","#E76BF3"))
jpeg(paste0(path,"/barchart_CBC_noL.jpeg"), units="in", width=10, height=7,res=600)
p
dev.off()

#------- CBC percentage barchat Lymphocytes/monoctes---------
# first Adjust cell number in later section
x =  x_y[,c("Lymphocytes.cbc","Monocytes.cbc")]
colnames(x) = gsub("\\.cbc","",colnames(x))
total = rowSums(x)
pcts = sapply(x, function(s) s /total )
x[,1:ncol(x)] = pcts
x$sample = rownames(x)
df1 <- gather(x, key = "Cell.Type", value = "percentage", -sample)
#df = melt(y,variable.name = "Cell.Type", value.name="counts")
p1 <-ggplot(df1, aes(sample, percentage))+
        geom_bar(stat = "identity", aes(fill = Cell.Type)) +
        ylab("percentage")+
        #ggtitle("Eosinophils + Neutrophils + Basophils+ Monocytes")+
        theme(text = element_text(size=15),    
              legend.position="none",
              plot.title = element_text(size=20,hjust = 0.5),
              axis.text.x = element_text(angle = 90, hjust = 1))+
        scale_fill_manual(values = c("#00BF7D","#E76BF3"))
jpeg(paste0(path,"/barchart_CBC.p.Lympho_mono.jpeg"), units="in", width=10, height=7,res=600)
p1
dev.off()

y = y[-which(y$sample == "Pt-17-C7"),]
sapply(y,mean) *mean(total_col)
###########################
# scRNA-seq 
###########################
#==== scRNA-seq cell number barchat Lymphocytes/monoctes =======
# before merge
y <- table(object@meta.data$singler1main,object@meta.data$orig.ident) %>%
        #y <- table(object@meta.data$singler1sub,object@meta.data$orig.ident) %>%
        as.data.frame %>% spread(Var2, Freq)
rownames(y) = y$Var1
y = y[,-1]
y = as.data.frame(t(y))
# after merge
y =  x_y[,-grep("\\.cbc",colnames(x_y))]
colnames(y) = gsub("\\.sc","",colnames(y))
#-----------------------
y[,"B_cells"] = y[,"B_cells"] +y[,"MCL"]
y = y[,-(which(colnames(y) %in% c("DC","Erythrocytes","Fibroblasts",
                                    "MCL","HSC","total")))]
#df = y
total = rowSums(y)
df = sapply(y, function(s) s/total) %>% as.data.frame()
colnames(df) = paste0(c("Lymphocytes:","","","Lymphocytes:","Lymphocytes:"),colnames(df))
df$sample = rownames(y)
df <- gather(df, key = "Cell.Type", value = "percentage", -sample)

p2 <-ggplot(df, aes(sample, percentage))+
        geom_bar(stat = "identity", aes(fill = Cell.Type)) +
        #ylab("cell number")+
        ylab("percentage")+
        #ggtitle("cell type summary of all samples in scRNA-seq")+
        theme(text = element_text(size=15),  
              legend.position = "none",
              plot.title = element_text(size=20,hjust = 0.5),
              axis.text.x = element_text(angle = 90, hjust = 1))+
        scale_fill_manual(values = c("#4DAF4A","#A65628","#377EB8",
                                     "#e94749","#FB9A99"))
jpeg(paste0(path,"/barchart_scRNA.LM.jpeg"), units="in", width=10, height=7,res=600)
p2
dev.off()

jpeg(paste0(path,"/barchart_CBC_vs_scRNA.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1,p2)+
        ggtitle("Compare Lymphocytes/Monocytes ratio in CBC vs scRNA-seq")+
        theme(text = element_text(size=15),  
              legend.position = "none",
              plot.title = element_text(size=20,hjust = 0.5))
dev.off()

#==== DJ only barchat=======
#df_DJ = df[df$sample =="DJ",] # with Lymphocytes
df_DJ <- gather(y[y$sample =="DJ",], key = "Cell.Type", value = "percentage", -sample)# without Lymphocytes
p1 <-ggplot(df_DJ, aes(sample, percentage))+
        geom_bar(stat = "identity", aes(fill = Cell.Type)) +
        ylab("percentage")+
        ggtitle("CBC of normal patient")+
        theme(text = element_text(size=15),     							
              plot.title = element_text(size=20,hjust = 0.5))+
        #scale_fill_manual(values = c("#00B0F6","#A3A500","#00BF7D","#F8766D","#E76BF3"))
        scale_fill_manual(values = c("#00B0F6","#A3A500","#F8766D","#E76BF3"))
jpeg(paste0(path,"/barchart_CBC_DJ~.jpeg"), units="in", width=10, height=7,res=600)
p1
dev.off()

(Mean = apply(y[,-c(1)],2,mean))
#Mean = c("mean",Mean)
#names(Mean)[1] ="sample"
#y = rbind(y, t(as.data.frame(Mean)))

#==== scRNA cell type barchat=======
object <- SetAllIdent(object, id = "orig.ident")
Normal_object <- SubsetData(object, ident.use = c("BH","DJ","MD","NZ"))
df_Normal <- table(Normal_object@meta.data$singler1main,Normal_object@meta.data$orig.ident) %>%
        #prop.table(margin = 2) %>% 
        as.data.frame %>% spread(Var2, Freq)

rownames(df_Normal) = df_Normal$Var1
df_Normal = df_Normal[,-1]
df_Normal["B_cells",] = df_Normal["B_cells",] +df_Normal["MCL",]
df_Normal = df_Normal[-which(rownames(df_Normal) =="MCL"),]

df_Normal$Cell.Type = rownames(df_Normal)
df <- gather(df_Normal, key = "sample", value = "counts", -Cell.Type)

p <-ggplot(df, aes(sample, counts))+
        geom_bar(stat = "identity", aes(fill = Cell.Type)) +
        #ylab("cell number in single RNA-seq")+
        ylab("cell percentage in single RNA-seq")+
        ggtitle("cell type summary of normal samples")+
        theme(text = element_text(size=15),     							
              plot.title = element_text(size=20,hjust = 0.5),
              axis.text.x = element_text(angle = 90, hjust = 1))+
        scale_fill_manual(values = c("#B3DE69","#FED9A6","#E41A1C","#FFFFCC",
                                       "#FF7F00","#e94749","#FB9A99","#A65628","#377EB8"))
jpeg(paste0(path,"/barchart_normal.p.jpeg"), units="in", width=10, height=7,res=600)
p
dev.off()


#==== T cells : NK cells : Monocytes percentage =======
df <- table(object@meta.data$singler1main,object@meta.data$orig.ident) %>%
        as.data.frame %>% spread(Var2, Freq)
df <- df[df$Var1 %in% c("Macrophages","Monocytes","NK_cells","T_cells"),]
(total_col <- colSums(df[,-1]))
pcts = apply(df[,-1],1, function(x) x / total_col)
df[,2:ncol(df)] = t(pcts)

colnames(df)[1] = "Cell.Type"
df <- gather(df, key = "sample", value = "percentage", -Cell.Type)
p <-ggplot(df, aes(sample, percentage))+
        geom_bar(stat = "identity", aes(fill = Cell.Type)) +
        #ylab("cell number in single RNA-seq")+
        ylab("percentage")+
        ggtitle("Macrophages","Monocytes, NK cells,T cells in all scRNA-seq")+
        theme(text = element_text(size=15),     							
              plot.title = element_text(size=20,hjust = 0.5),
              axis.text.x = element_text(angle = 90, hjust = 1))+
        scale_fill_manual(values = c("#e94749","#FB9A99","#A65628","#377EB8"))
jpeg(paste0(path,"/barchart_all.p.jpeg"), units="in", width=10, height=7,res=600)
p
dev.off()


############################
# Adjust cell number
############################

#==== scRNA-seq percentage barchat Lymphocytes/monoctes =======
x =  CBC
rownames(x) =x$sample
x = x[,-1]

colnames(x) = paste0(colnames(x),".cbc")

#y <- table(object@meta.data$singler1main,object@meta.data$orig.ident) %>%
y <- table(object@meta.data$singler1sub,object@meta.data$orig.ident) %>%
        as.data.frame %>% spread(Var2, Freq)
rownames(y) = y$Var1
y = y[,-1]
y = as.data.frame(t(y))
y$total = rowSums(y)
colnames(y) = paste0(colnames(y),".sc")

x_y = merge(x, y, by ="row.names")#,all.y = TRUE)
head(x_y)
rownames(x_y) = (x_y$Row.names)
x_y = x_y[,-1]

# Find L+M cbc value to adjust
x_y$L_M.cbc = x_y$Lymphocytes.cbc + x_y$Monocytes.cbc
mean(x_y$L_M.cbc);median(x_y$L_M.cbc)
#use MD as control
x_y$L_M.cbc = x_y$L_M.cbc / x_y["MD","L_M.cbc"]

#==== scRNA-seq cell number barchat Before ajusting =======
y =  x_y[,-grep("\\.cbc",colnames(x_y))]
colnames(y) = gsub("\\.sc","",colnames(y))
y = y[,-(which(colnames(y) %in% c("total")))] #
y = y[,-(which(colnames(y) %in% c("total","MCL")))] # remove MCL
dim(y)
y$sample = rownames(y)

y <- gather(y, key = "Cell.Type", value = "cell.number", -sample)

singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
singler_colors1[duplicated(singler_colors1)];singler_colors2[duplicated(singler_colors2)]
length(singler_colors1);length(singler_colors2)


g1 <-ggplot(y, aes(sample, cell.number))+
        geom_bar(stat = "identity", aes(fill = Cell.Type)) +
        ylab("cell number")+
        ggtitle("cell type summary of in scRNA-seq before adjusting cell number")+
        theme(text = element_text(size=15),  
              legend.position = "none",
              plot.title = element_text(size=20,hjust = 0.5),
              axis.text.x = element_text(angle = 90, hjust = 1))+
        #scale_fill_manual(values = c(singler_colors1)) # include MCL
        #scale_fill_manual(values = c(singler_colors1[-7])) # no MCL
        #scale_fill_manual(values = c(singler_colors2))
        scale_fill_manual(values = c(singler_colors2[-15])) # no MCL
jpeg(paste0(path,"/barchart_scRNA.jpeg"), units="in", width=10, height=7,res=600)
g1
dev.off()

#==== scRNA-seq cell number barchat After ajusting =======
# after adjust cell number=============================
y =  x_y[,-grep("\\.cbc",colnames(x_y))]
colnames(y) = gsub("\\.sc","",colnames(y))
#y = y[,-(which(colnames(y) %in% c("total")))] #
y = y[,-(which(colnames(y) %in% c("total","MCL")))] #
y.adj = apply(y,2, function(s) s*x_y$L_M.cbc) %>% as.data.frame
y.adj$sample = rownames(y)
y.adj <- gather(y.adj, key = "Cell.Type", value = "cell.number", -sample)

g2 <-ggplot(y.adj, aes(sample, cell.number))+
        geom_bar(stat = "identity", aes(fill = Cell.Type)) +
        ylab("cell number")+
        ggtitle("After")+
        theme(text = element_text(size=15),  
              legend.position = "none",
              plot.title = element_text(size=20,hjust = 0.5),
              axis.text.x = element_text(angle = 90, hjust = 1))+
        #scale_fill_manual(values = c(singler_colors1)) # include MCL
        #scale_fill_manual(values = c(singler_colors1[-7])) # no MCL
        #scale_fill_manual(values = c(singler_colors2))
        scale_fill_manual(values = c(singler_colors2[-15])) # no MCL

jpeg(paste0(path,"/barchart_scRNA~.jpeg"), units="in", width=10, height=7,res=600)
g2
dev.off()
        
jpeg(paste0(path,"/barchart_scRNA_adjusting_sub.noMCL.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g1,g2)+
        ggtitle("Before and after adjusting scRNA-seq cell number")+
        theme(text = element_text(size=15),  
              legend.position = "none",
              plot.title = element_text(size=20,hjust = 0.5))
dev.off()


# fill up missing value--------------------
x_y$m.cbc_L = x_y$Monocytes.sc / x_y$WBC_L.cbc/1000000*100
(m.cbc_L <- median(x_y[which(!is.na(x_y$WBC.cbc) & (x_y$Monocytes.sc >50)),"m.cbc_L"]))
hist(x_y[which(!is.na(x_y$WBC.cbc) & (x_y$Monocytes.sc >50)),"m.cbc_L"])
x_y[which(is.na(x_y$WBC_L.cbc)),]
x_y[which(is.na(x_y$WBC_L.cbc)),"m.cbc_L"] = m.cbc_L
x_y[which(is.na(x_y$WBC_L.cbc)),"WBC_L.cbc"] = 
        x_y$Monocytes.sc[which(is.na(x_y$WBC_L.cbc))] / m.cbc_L/1000000*100

x_y
x_y$L_T_ratio.x = x_y$Lymphocytes.x/(x_y$Monocytes.x+x_y$Lymphocytes.x)
x_y$L_T_ratio.y = x_y$Lymphocytes.y/(x_y$total)

x_y[,c("sample","L_T_ratio.x","L_T_ratio.y")]%>%
        kable %>% kable_styling()

x_y  %>% kable %>% kable_styling()

x_y$Lymphocytes.x + x_y$Monocytes.x




z <- table(object@meta.data$singler1sub,object@meta.data$orig.ident) %>%
        as.data.frame %>% spread(Var2, Freq)
rownames(z) = z$Var1
z = z[,-1]
z = as.data.frame(t(z))
z$total <- rowSums(z)

x_y_z = merge(x_y[,c("sample","L_T_ratio.x")], z, by ="row.names")
rownames(x_y_z) = x_y_z$Row.names
x_y_z = x_y_z[,-c(1:2)]
x_y_z$adj_total = x_y_z$L_T_ratio.x *x_y_z$total 
x_y_z[,c("total","adj_total")] %>% kable %>% kable_styling()

