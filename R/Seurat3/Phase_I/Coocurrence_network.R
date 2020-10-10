# https://zhuanlan.zhihu.com/p/170466920
#载入所需R包；
library(igraph)
library(Hmisc)
library(networkD3)
orig.idents <- "Pt25_SB1"
source("Rshiny/B_MCL_41_Pt25_network/util.R")
#工作目录设置；

#读入OTU绝对丰度表；
otu=read.table("otu_data.xls" ,header=T,row.names = 1,sep = "\t")
#转成矩阵，之后的相关性计算需要矩阵对象；
otu <-as.matrix(otu)
dim(otu)
head(otu)

#将丰度值大于1的值替换为1，便于计算不同otu的检测率；
dt <-otu
dt[dt>1] <-1
#将样本发现率低于20%的otu过滤掉；
no <-which(rowSums(dt)/ncol(dt)>0.2)
length(no)
otu <-otu[no,]

#相关性系数计算 =================
#主要利用Hmisc包的rcorr()函数计算spearman相关系数，当然，也可计算pearson相关系数。rcorr()函数默认是对列进行两两计算，故须转置，计算otu间的相关性。
#计算相关性系数；
sp.cor <- rcorr(t(otu),type="spearman")

#提取r、p值矩阵；
load("Rshiny/B_MCL_41_Pt25_network/data/B_MCL_41_Pt25_network~.Rda")

r.cor <- cor_list$Pt25_SB1
p.cor <- pvalue_list$Pt25_SB1

#使用Benjamini-Hochberg("FDR-BH")法进行多重检验校正；
p.adj  <-  p.adjust(p.cor, method="BH")

#指定阈值；
r.cutoff=0.5
p.cutoff=0.001

#copy files
r.matrix <- as.matrix(r.cor)
p <- p.adj

#将矩阵中不符合条件的r值替换为0；
r.matrix[which(abs(r.matrix) <= r.cutoff)]=0
r.matrix[which(p>p.cutoff)]=0

M = lower.tri(r.matrix, diag = T)
r.matrix[M] = t(r.matrix)[M]# fill up Triangular matrix
dim(r.matrix)

#删掉相关系数矩阵数据全都为0的行和列；
r.matrix <- r.matrix[which(rowSums(abs(r.matrix))!=0),
                     which(colSums(abs(r.matrix))!=0)]
#重新添加对角线 的 1
diag(r.matrix) = 1

#查看过滤后的矩阵；
dim(r.matrix)
table(rownames(r.matrix) == colnames(r.matrix))
r.matrix[1:7,1:7]

# 生成网络图 ===================
# 从结构来看，对称矩阵是展示网络图非常合适的数据结构，
# 比如相关系数矩阵是两个OTU的丰度相关系数，
# 而网络图恰是Source和Target两个OTU结点的连线，相关系数可作为权重。
# 如果是n维矩阵（即n个OTU），除去对角线，共有n(n-1)/2个关系对，当然用组合公式也可算得。

# Convert to igraph
#使用邻接矩阵（即相关系数矩阵）创建网络；
g1 <- graph.adjacency(r.matrix,weight=T,mode="undirected")

# Remove duplicate edges
#去掉冗余的边（multiple edges、loop edges）；
g1 <- simplify(g1)

# Find group membership
#计算群体结构（short random walks）；
wt <- cluster_walktrap(g1, steps = 6)
members <- membership(wt)
#使用默认颜色列表；
V(g1)$color  <- c$membership+1

# Convert igraph to list for networkD3
sj_list <- igraph_to_networkD3(g1, group = members)

# Plot as a forceDirected Network
plot <- forceNetwork(Links = sj_list$links, Nodes = sj_list$nodes,
             colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"),
             Source = 'source',
             Target = 'target', NodeID = 'name', Group = 'group',
             fontSize = 20,
             zoom = TRUE, linkDistance = 50)

# 导出graph对象 =================
#将网络图导出为"graphml"、"gml"格式，方便导入Gephi中使用；
write_graph(plot, "g1.graphml", format="graphml")
write_graph(g1, "g1.gml", format="gml")
#支持的格式有"edgelist", "pajek", "ncol", "lgl","graphml", "dimacs", "gml", "dot", "leda"), ...
