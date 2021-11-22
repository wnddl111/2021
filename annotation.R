library(SeuratObject)
library(Seurat)
library(harmony)
setwd("/Users/juyoung/Desktop/immunotherapy_jy/singlecell_jy/")

combined=readRDS("./data/harmony_combined_neighbor_0.05.rds")
cluster_dist=as.data.frame(table(combined$seurat_clusters)) #총 23개 각 cluster별로 1000개씩 랜덤하게 추출

##############
#data size down
##############
#data크기가 1만개가 넘어가는 cell에 대해서만 down시켜줌,, 통계적으로 맞을 거란 생각은 안들지만 

set.seed(1)
num=0
for (i in cluster_dist$Freq){
  if(i>10000){
    random=sample(x=1:i,size=i*0.1,replace = F)
    assign(paste0('cluster',num),combined[,which(combined$seurat_clusters==num)][,random])
  }
  
  num=num+1
}
#

#1개이상인 cluster0-4번이랑 나머지 merge
left_combined=combined[,which(combined$seurat_clusters %in% c(5:23))]
obj_list=c(cluster1,cluster2,cluster3,cluster4,left_combined)

down_combined=merge(cluster0,y=obj_list, project= 'down_combined')
down_combined#44484개 cell로 줄었땅
##############
#basic
##############

down_combined_test=NormalizeData(down_combined)
down_combined_test=FindVariableFeatures(down_combined_test, selection.method = 'vst', nfeatures=2000)
all.genes <- rownames(down_combined_test)
down_combined_test <- ScaleData(down_combined_test, features =all.genes)
down_combined_test <- RunPCA(down_combined_test, features = VariableFeatures(object=down_combined_test))

ElbowPlot(down_combined_test)
down_combined_test <- FindNeighbors(down_combined_test, dims=1:20)
down_combined_test <- FindClusters(down_combined_test, resolution=0.5)
down_combined_test <- RunUMAP(down_combined_test, dims = 1:20)
DimPlot(down_combined_test, reduction = "umap")

down_combined_test$seurat_clusters

table(colnames(down_combined)==colnames(down_combined_test))#all true
down_combined_test$orig_cluster=down_combined$seurat_clusters
#test에 ORIG CLUSTER추가해서 그림그릴 때 보여줄 수 있으려나?
##############
#singleR
##############

#refernece data가 필요함
#BiocManager::install('celldex')
#BiocManager::install('scRNAseq')
library(celldex)
library(SingleR)
library(scRNAseq)
hpca.se <- HumanPrimaryCellAtlasData()

#seurat obj가 아니라 singlecellexperiment가 input임
down_combined_test_se=as.SingleCellExperiment(down_combined_test)
colLabels(down_combined_test_se)=down_combined_test_se$seurat_clusters #이걸통해서 cluster별로 지정할 수 있음 
#assay.type.ref = An integer scalar or string specifying the assay of ref containing the relevant expression matrix
pred.down.combined_test <-SingleR(test=down_combined_test_se, ref=hpca.se, assay.type.ref = 1,
                             labels= hpca.se$label.main,
                             clusters=colLabels(down_combined_test_se))
pred.down.combined_test.table<-as.data.frame(pred.down.combined_test) #각 row는 single cell에 대한 예측결과, first.label이 미세 조정 전 결과, label이  조정 후 결과 
table(pred.down.combined_test$labels)
a=as.data.frame(down_combined_test$seurat_clusters,down_combined_test$orig_cluster)
View(a)
new.cluster.ids <- pred.down.combined_test$labels
#names(new.cluster.ids) <-levels(down_combined_test)
#down_combined_test <- RenameIdents(down_combined_test, new.cluster.ids)
Idents(down_combined_test)=new.cluster.ids
DimPlot(down_combined_test, reduction = 'umap', label=T, pt.size = 0.5)

pred.down.combined

