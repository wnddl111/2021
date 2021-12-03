data(16_cluster_final.rda) #최종파일 n_final로 저장됨

final.markers <- FindAllMarkers(n_final, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25) #전체 마커 찾음 -> python돌

###############
#cd8+ cell로 명명된 것들 끼리 비교해서 marker찾음 -> python sub_annotation~
###############

cluster_markers_name <- c('cluster0_include.markers','cluster16_include.markers','cluster14_include.markers', 'cluster7_include.markers',
                          'cluster6_include.markers','cluster5_include.markers','cluster4_include.markers',
                          'cluster1_include.markers','cluster8_include.markers')

cluster_num <- as.numeric(gsub('\\D','',cluster_markers_name))

for (markers in cluster_markers_name){
  ident1=as.numeric(gsub('\\D','',markers))
  ident2=cluster_num[-which(cluster_num==ident1)]
  
  assign(markers, FindMarkers(n_final,ident.1 = ident1, ident.2 = ident2, min.pct = 0.25))
  
  write.csv(get(markers), file=paste0('./data/all_find_list/all_list_',markers,'.csv'))
  assign(markers, as.data.frame(get(markers) %>%
                                  slice_max(n=30, order_by = avg_log2FC)))
  assign(paste0(markers,'_genes'), rownames(get(markers)))
  write.csv(get(paste0(markers,'_genes')), file=paste0('./data/all_find_list/',markers,'.csv'))
}

###############
#NK cell
###############

cluster_markers_name <- c('cluster13_include.markers','cluster15_include.markers')

cluster_num <- as.numeric(gsub('\\D','',cluster_markers_name))

for (markers in cluster_markers_name){
  ident1=as.numeric(gsub('\\D','',markers))
  ident2=cluster_num[-which(cluster_num==ident1)]
  
  assign(markers, FindMarkers(n_final,ident.1 = ident1, ident.2 = ident2, min.pct = 0.25))
  
  write.csv(get(markers), file=paste0('./data/all_find_list/all_list_',markers,'.csv'))
  assign(markers, as.data.frame(get(markers) %>%
                                  slice_max(n=30, order_by = avg_log2FC)))
  assign(paste0(markers,'_genes'), rownames(get(markers)))
  write.csv(get(paste0(markers,'_genes')), file=paste0('./data/all_find_list/',markers,'.csv'))
}

###############
#B cell
###############

cluster_markers_name <- c('cluster3_include.markers','cluster11_include.markers')

cluster_num <- as.numeric(gsub('\\D','',cluster_markers_name))

for (markers in cluster_markers_name){
  ident1=as.numeric(gsub('\\D','',markers))
  ident2=cluster_num[-which(cluster_num==ident1)]
  
  assign(markers, FindMarkers(n_final,ident.1 = ident1, ident.2 = ident2, min.pct = 0.25))
  
  write.csv(get(markers), file=paste0('./data/all_find_list/all_list_',markers,'.csv'))
  assign(markers, as.data.frame(get(markers) %>%
                                  slice_max(n=30, order_by = avg_log2FC)))
  assign(paste0(markers,'_genes'), rownames(get(markers)))
  write.csv(get(paste0(markers,'_genes')), file=paste0('./data/all_find_list/',markers,'.csv'))
}

#####################

