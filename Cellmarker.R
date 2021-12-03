library(Seurat)
library(SeuratObject)
library(dplyr)

setwd('./singlecell_jy/')

data(final_singleR_annotation2) #singleR 돌려서 annotation2로 한 seurat obj (4만개)
data("pred.final.down2") #annotation2 singleR결과 

as.data.frame(pred.final.down2$labels) 
Idents(final)=final$seurat_clusters #이렇게 하는 거 말고 좋은 방안이있나?/ findallmarker로 돌리려면 cluster가 ident가 되야해서 바꿈,,,

final.markers <- FindAllMarkers(final, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
table(findMarker$avg_log2FC > 0)#all true == only.pos=T

findMarker=as.data.frame(final.markers %>%
  group_by(cluster) %>%
  slice_max(n=30, order_by= avg_log2FC)) #5개는 너무 적어서 30개로 바꿈

table(findMarker$cluster)#cluster23의 경우 avg_log2FC가 positive한 경우가 별로 없나봄 3개뽑힘 (ㄴㄴ cluster에 해당되는 cell이 21개 밖에 없어서 그럼)



findMarker_subset=findMarker[,c('cluster','gene')]
write.csv(findMarker_subset,file='./data/findMarker.csv')

findMarker_subset = read.csv('./data/findMarker.csv')
all_findMarker=read.csv('./data/server_all_findMarker.csv')
table(all_findMarker$gene==findMarker_subset$gene)#f 64개나 있음 
table(all_findMarker$cluster)



#server_all_findMarker -> n=5개로 해서 server에서 돌린 결과 
