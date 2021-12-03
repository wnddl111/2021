##############
#singleR
#합치는 방법이 총 3갠데 2개를 시도해봤고 최종적으로 2
##############

#refernece data가 필요함
#BiocManager::install('celldex')
#BiocManager::install('scRNAseq')
setwd('./singlecell_jy/')
data("harmony_combined_40000cell")#4만개로 줄인 data사용/ 나중에 server에는 전체 데이터 사용함 

library(celldex)
library(SingleR)
library(scRNAseq)

#####################
#ref2
hpca.se <- HumanPrimaryCellAtlasData() #두번째 annotation용 방법
bpe <- BlueprintEncodeData()#두번째 annotation용 방법

#singleR돌려서 예측한 label확인
pred.final.down2 <- SingleR(test=final_se, ref=list(BPE=bpe, HPCA=hpca.se), assay.type.ref = 1,
                            labels= list(bpe$label.main,hpca.se$label.main),
                            clusters=colLabels(final_se))

table(pred.final.down2$labels)
table(pred.final.down2$reference) #1 -> 21 / 2 -> 3 1이 이김(bpe) 어떤 ref를 기준으로 annotation이 많이 됐는지 보는 것 

head(pred.final.down2$orig.results$BPE$labels)
head(pred.final.down2$orig.results$HPCA$labels)


plotScoreHeatmap(pred.final.down2)
plotDeltaDistribution(pred.final.down2)
save(pred.final.down2, file='./data/pred.final.down2.rda')

new.cluster.ids <- pred.final.down2$labels #예측한 파일 기준으로 ident 바꿈 
names(new.cluster.ids) <-levels(final)
final <- RenameIdents(final, new.cluster.ids)

save(final, file='./data/final_singleR_annotation2.rda')


#####################
#ref1
hpca2 <- hpca.se
hpca2$label.main <- paste0('HPCA.', hpca2$label.main) #어느 ref에서 온건지 헷갈리지 않으려고 설정하는 것 

bpe2 <- bpe
bpe2$label.main <- paste0('BPE.', bpe2$label.main)

shared <- intersect (rownames(hpca2), rownames(bpe2))
length(shared) #16267
dim(hpca2)#19363
dim(bpe2)#19859

union <- cbind(hpca2[shared,], bpe2[shared,])

#seurat obj가 아니라 singlecellexperiment가 input임

final_se=as.SingleCellExperiment(final)
colLabels(final_se)=final_se$seurat_clusters #이걸통해서 cluster별로 지정할 수 있음 


#assay.type.ref = An integer scalar or string specifying the assay of ref containing the relevant expression matrix
pred.final.down <-SingleR(test=final_se, ref=union, assay.type.ref = 1,
                                  labels= union$label.main,
                                  clusters=colLabels(final_se))
pred.final.down.table<-as.data.frame(pred.final.down) #각 row는 single cell에 대한 예측결과, first.label이 미세 조정 전 결과, label이  조정 후 결과 
table(pred.final.down$labels)
table(pred.final.down$pruned.labels)#이렇게 통합하면 각 ref에 대한 batch effect가 존재
#result -> 이렇게 하면 ref별로 따로 naming을 해서 annotation이 어렵고
#the marker set is likely to contain genes responsible for uninteresting batch effects
#cell type이 아니라 batch effect에 의해서 차이나는 유전자가 있을 수 있다 그말이래 


##############
##############

getwd()
png("./figure/singleR_annotation2.png",width=1000)
p1 <- DimPlot(final, reduction = 'umap', label=T, pt.size = 0.5)
p1#그림이 안보이는 거지 그려짐 ! 저장됨 
dev.off()

png("./figure/singleR_annotation_cluster2.png", width=1000)
p1 <- DimPlot(final,group.by = 'seurat_clusters', reduction = 'umap', label=T, pt.size = 0.5)
p1#그림이 안보이는 거지 그려짐 ! 저장됨 
dev.off()

#singleR_annotation.pdf -> singleR돌려서 cluster name 대신에 cell type으로 annotation한 결과
#singleR_annotation_cluster -> 비교를 위해서 cluster name으로 annotation한 거, labeling함
#singleR_annotation_cluster_F-> labeling 안함, 겹치는 게 있어서 
#끝에 2 붙은 것들은 이제 HPCA+BPE 2번째 방법으로 합쳤을 때 결과임

