#https://swaruplab.bio.uci.edu/tutorial/integration/integration_tutorial.html#harmony

library(Seurat)
library(harmony)
setwd('/data/juyoung/scRNAseq/for_local/immunotherapy_jy/singlecell_jy')

l1=readRDS('./l1_last.rds')
g1=readRDS('./g1_last.rds')
bcc=readRDS('./bcc_last.rds')
scc=readRDS('./scc_last.rds')
mel=readRDS('./mel_last.rds')

#min.cell, min.feature -> create seurat obj
#mel,bcc,scc 는 count data 통으로 진행
#bcc,scc sample별로 하니까 유전자 소실 큼, mel도 sample이 적어서 따로 하는 건 불가
#g1,l1은 sample 별로 진행=> project별로 변경.(0103)

#미토콘드리아도 각각 따로 진행 -> mel은 유전자가 m인게 없었음
obj_list=c(l1,g1,bcc,scc,mel)
#202112228 미래에서 왔다 공통 유전자만 뽑아서 integration해야하는데 합집합 유전자로 뽑아버렸다 -- 
gene_list=c(rownames(l1))
for ( i in obj_list){
  gene_list = intersect(gene_list, rownames(i))
}
length(gene_list)#12395->12748(G1,L1도 PROJECT별로 했더니 올라감)

#fififififififinal obj 
final_l1=l1[gene_list,]
final_g1=g1[gene_list,]
final_bcc=bcc[gene_list,]
final_scc=scc[gene_list,]
final_mel=mel[gene_list,]

#반응성 데이터 합치기 
final_g1$response=ifelse(final_g1$sample=='40784154','responder','non-responder')
table(final_g1$sample) #40784154 - 7409

final_l1$response = ifelse(final_l1$sample=='39745926','responder','non-responder')
table(final_l1$sample) #'39745926' 1180

final_bcc$response = ifelse(final_bcc$sample %in% c('su005','su06','su007','su008','su010'),'non-responder','responder')
table(final_bcc$sample %in% c('su005','su06','su007','su008','su010')) #11747

final_scc$response = ifelse(final_scc$sample %in% c('su013','su014'),'non-responder','responder')
table(final_scc$sample %in% c('su013','su014')) #5400

final_mel$response = ifelse(final_mel$sample=='Mel04.3','responder','non-responder')
table(final_mel$sample =='Mel04.3') #78

table(final_g1$response) #responder 7409
table(final_l1$response) #responder 1180
table(final_bcc$response) #non-responder 11747
table(final_scc$response) #non-responder 5400
table(final_mel$response) #responder 78

#새롭게 OBJ_LIST 만들어 주기 
obj_list=c(final_l1,final_g1,final_bcc,final_scc,final_mel)

combined=merge(obj_list[[1]],y=obj_list[2:5], project= 'combined')
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = 'vst', nfeatures = 2000)

combined <-ScaleData(combined, verbose=F) #check 안되서 SERVER로
combined <-RunPCA(combined, npcs=30, verbose=F)
#noF_combined<-RunPCA(noF_combined,npcs-30, verbose=F)

combined <- RunUMAP(combined, reduction='pca', dims = 1:20)
combined <- FindNeighbors(combined, reduction='pca')
combined <- FindClusters(combined, resolution = 0.05)

save(combined,file='./data/harmony_combined_from_server_0102.rda')

#data("harmony_combined_from_server_noNF")
#noF_combined <- combined
#data("harmony_combined_from_server") #비교해보니까 차이가난다!! 

#plotting settings
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ggplot2)
theme_set(theme_cowplot())

umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)

fig_dir='./figure/'

pdf(paste0(fig_dir, "qc_violin_plot.pdf"), width=30, height=10)
VlnPlot(combined,group.by='cancer_type',features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#harmony
combined#12748 221753
meta=cbind(combined$orig.ident,combined$platform,combined$type,combined$treatment.group,combined$cancer_type,combined$sample, combined$response)
meta=as.data.frame(meta)

colnames(meta)=c('orig.ident','platform','type','treatment.group','cancer_type','sample','response')
#combined <- AddMetaData(combined,metadata=as.character(meta[colnames(combined),]$platform),col.name = 'platform')
combined <- AddMetaData(combined,metadata = meta)

#dim 정하기
pdf(paste0(fig_dir, "noNF_pca_harmony.pdf"), width=10, height=10)
p1 <- DimHeatmap(combined,dims = 1:20, cells=500,balanced=T)
p1
dev.off()

#combined <- JackStraw(combined, num.replicate = 100) 
#combined <- ScoreJackStraw(combined, dims=1:20) #최ㅐ대 20
#pdf(paste0(fig_dir, "pca_jack_harmony.pdf"), width=10, height=10)
#p1 <- JackStrawPlot(combined,dims=1:20)
#p1
#dev.off()

#pdf(paste0(fig_dir,'elbow_harmony.pdf'),width=10,height=10)
#p1 <- ElbowPlot(combined)
#p1
#dev.off()#15/10쯤 하면 될 것 같은 느낌 => 교수님은 20까지 해도 될 것 같대
#그래야 rare type 볼 수 있다고 그럼 이것도 따로 해놓아야겠다

combined <- RunHarmony(combined,'cancer_type')#0102
#dim will change
combined <- RunUMAP(combined,reduction='harmony',dims=1:20)
#Harmony converged after 9 iterations
#Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.RNA.harmony; see ?make.names for more details on syntax validity
#There were 16 warnings (use warnings() to see them)

combined <- FindNeighbors(combined,reduction='harmony',dims = 1:20) #각 세포마다 distance를 구해서 비슷한 cell type인지 아닌지 찾음
combined <- FindClusters(combined, resolution=0.05)
combined#12748 221753
#combined@meta.data
pdf(paste0(fig_dir, "pc20_run_platrom_harmony_PLATFORM_neighbor_reol0.05_0102.pdf"), width=20, height=10)
p1 <- DimPlot(combined, split.by = 'platform', reduction = 'umap',raster=F)
p1#그림이 안보이는 거지 그려짐 ! 저장됨 
dev.off()

pdf(paste0(fig_dir, "pc20_run_platrom_harmony_CANCER_neighbor_0102.pdf"), width=30, height=10)
p1 <- DimPlot(combined, split.by = 'cancer_type', reduction = 'umap',raster=F)
p1#그림이 안보이는 거지 그려짐 ! 저장됨 
dev.off()

pdf(paste0(fig_dir, "run_platrom_harmony_treatment_neighbor_0102.pdf"), width=20, height=10)
p1 <- DimPlot(combined, split.by = 'treatment.group', reduction = 'umap',raster=F)
p1#그림이 안보이는 거지 그려짐 ! 저장됨 
dev.off()

save(combined,file='./data/harmony_combined_neighbor_0.05_0102.rda')
#saveRDS(combined, './data/harmony_combined_neighbor_0.05_0102.rds')
