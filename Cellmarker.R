library(Seurat)
library(SeuratObject)
library(dplyr)

setwd('./singlecell_jy/')

data(final_singleR_annotation2) #singleR 돌려서 annotation2로 한 seurat obj
data("pred.final.down2") #annotation2 singleR결과 

as.data.frame(pred.final.down2$labels)
Idents(final)=final$seurat_clusters #이렇게 하는 거 말고 좋은 방안이있나? 
final.markers <- FindAllMarkers(final, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

findMarker=as.data.frame(final.markers %>%
  group_by(cluster) %>%
  slice_max(n=5, order_by= avg_log2FC))

table(findMarker$cluster)#cluster23의 경우 avg_log2FC가 positive한 경우가 별로 없나봄 ㄷ3개뽑힘
(findMarker$gene)
table(findMarker$avg_log2FC > 0)#all true

findMarker_subset=findMarker[,c('cluster','gene')]
write.csv(findMarker_subset,file='./data/findMarker.csv')

findMarker_subset = read.csv('./data/findMarker.csv')
all_findMarker=read.csv('./data/server_all_findMarker.csv')
table(all_findMarker$gene==findMarker_subset$gene)#f 64개나 있음 
table(all_findMarker$cluster)



#server_all_findMarker -> n=5개로 해서 server에서 돌린 결과 