#https://bioconductor.org/packages/devel/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf
setwd('D:/')
vst = read.table('./data/89_F0405_for_heatmap_0.3.txt',sep='\t')
#col = read.table(file='./data2/my_data_colname_90.txt')
#col=col[-which(col$V1=='U3ZG'),]
#write.table(col,file='./data/my_data_colname_89.txt')
col = read.table('./data/my_data_colname_89.txt',sep = '\t')

gene=vst[,1]
rownames(vst)=gene
vst = vst[,-1]
colnames(vst)=col$V1

#crpc만 뽑기(사실은hspc섞여있음)
vst = vst[,1:36]
colnames(vst) = col[1:36,]

#all ccrg gene
our_ccrg_gene = as.data.frame(readxl::read_excel('./data2/Data_Sheet_1_Molecular Classification Based on Prognostic and Cell Cycle-Associated Genes in Patients With Colon Cancer.xlsx ', sheet=4))
our_ccrg_gene

#
dim(vst[rownames(vst) %in% our_ccrg_gene$Gene,])#264 36
d=as.matrix(vst[rownames(vst) %in% our_ccrg_gene$Gene,])
#BiocManager::install('ConsensusClusterPlus')
library(ConsensusClusterPlus)

colnames(vst)
title = 'CRPC subtyping'

results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
#results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=123456789,plot="png")
#results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=542789435,plot="png")



icl = calcICL(results,title=title,plot="png")
icl[['clusterConsensus']]

results

library(pheatmap)
results[[3]]
annot = as.data.frame(sapply(results[[3]]$consensusClass,function(x) ifelse(x==1,'cluster1',
                                                                            ifelse(x==2,'cluster2','cluster3'))))
colnames(annot)='Cluster'
table(annot$Cluster)

result_k3=results[[3]]$consensusMatrix
colnames(result_k3) = rownames(annot)

pheatmap(result_k3, annotation_col = annot)
clusters = as.data.frame(results[[3]]$consensusClass)
colnames(clusters)='cluster'
write.table(clusters, './data/clustering_89.txt',sep = '\t')



#annot = as.data.frame(sapply(results[[5]]$consensusClass,function(x) ifelse(x==1,'cluster1',
                                                                            ifelse(x==2,'cluster2',
                                                                                   ifelse(x==3,'cluster3',
                                                                                          ifelse(x==4,'cluster4',
                                                                                                 ifelse(x==5,'cluster5','cluster6')))))))
annot=as.data.frame(sapply(results[[2]]$consensusClass, function(x) ifelse(x==1,'cluster1','cluster2')))
colnames(annot)='Cluster'
z=results[[2]]$consensusMatrix
colnames(z)=rownames(annot)
pheatmap(z,annotation_col = annot)


