#R studio server에서 진행
#최종 20만개 파일로 진행
#16-23 cluster는 cluster 16으로 합침

#.libPaths('/data/juyoung/test_library/x86_64-apple-darwin17.0/')
.libPaths('/opt/R/4.1.1/lib/R/library')

setwd('/data/juyoung/scRNAseq/for_local/immunotherapy_jy/singlecell_jy/')

library(scRNAseq)
library(Seurat)
library(SeuratObject)
library(SingleR)
#BiocManager::install('celldex')
library(celldex)

hpca.se <- HumanPrimaryCellAtlasData()
bpe <- BlueprintEncodeData()

#20만개 
final=readRDS('/data/juyoung/scRNAseq/for_local/immunotherapy_jy/singlecell_jy/data/harmony_combined_neighbor_0.05.rds')
final$seurat_clusters
final_se=as.SingleCellExperiment(final)

colLabels(final_se) = final_se$seurat_clusters

pred.final.down.all <- SingleR(test=final_se, ref=list(BPE=bpe, HPCA=hpca.se), assay.type.ref = 1,
                                                   labels= list(bpe$label.main,hpca.se$label.main),
                                                   clusters=colLabels(final_se))

pred.final.down.all$labels
getwd()
save(pred.final.down.all, file='./data/pred.final.down.all.rda')

#cluster 23개 -> 16개로 줄임
library(dplyr)
n_final=final
new.cluster.ids <-c(seq(0,16),rep(16,7))
length(new.cluster.ids) #24

names(new.cluster.ids) <- levels(n_final)
n_final <- RenameIdents(n_final, new.cluster.ids)
DimPlot(n_final, reduction = 'umap', label = T, pt.size = 0.5, raster=F)+ NoLegend()

save(n_final,file='./data/16_cluster_final.rda') #최종파일

n_final.markers <- FindAllMarkers(n_final, only.pos=T, min.pct = 0.25, logfc.threshold = 0.25)

n_final_se = as.SingleCellExperiment(n_final)
colLabels(n_final_se) =n_final_se$ident

pred.final.down.all.16 <- SingleR(test=n_final_se, ref=list(BPE=bpe, HPCA=hpca.se), assay.type.ref = 1,
                            labels= list(bpe$label.main,hpca.se$label.main),
                            clusters=colLabels(n_final_se))

pred.final.down.all.16.data= as.data.frame(pred.final.down.all.16$labels)
save(pred.final.down.all.16.data, file='./data/pred.final.down.all.16.rda')


new.cluster.ids <- pred.final.down.all.16$labels
names(new.cluster.ids) <-levels(n_final)
n_final <- RenameIdents(n_final, new.cluster.ids)

DimPlot(n_final, reduction='umap',label=T, pt.size = 0.5, raster=F)
