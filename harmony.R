#https://swaruplab.bio.uci.edu/tutorial/integration/integration_tutorial.html#harmony

library(Seurat)
library(harmony)
setwd('/Users/juyoung/Desktop/immunotherapy_jy/singlecell_jy')

l1=readRDS('./l1_last.rds')
g1=readRDS('./g1_last.rds')
bcc=readRDS('./bcc_last.rds')
scc=readRDS('./scc_last.rds')
mel=readRDS('./mel_last.rds')

#min.cell, min.feature -> create seurat obj
#mel,bcc,scc 는 count data 통으로 진행
#bcc,scc sample별로 하니까 유전자 소실 큼, mel도 sample이 적어서 따로 하는 건 불가
#g1,l1은 sample 별로 진행

#미토콘드리아도 각각 따로 진행 -> mel은 유전자가 m인게 없었음
obj_list=c(l1,g1,bcc,scc,mel)
#check
#이부분은 없는데 해줘야할 것 같아서 함...
obj_list <- lapply(X=obj_list, FUN=function(x){
  x <-NormalizeData(x) #이거는 반영이 됐나? 어떻게 확인하지? 이것도 결국 merge하고 다시함
  x <-FindVariableFeatures(x,selection.method='vst',nfeatures=2000)
  #각각 찾은 건 merge에 반영이 안되네
})

combined=merge(obj_list[[1]],y=obj_list[2:5], project= 'combined')

#combined <-ScaleData(combined, verbose=F) #check 안되서 SERVER로
save(combined,file='./data/harmony_combined_to_server.rda')
data()
combined <-RunPCA(combined, npcs=30, verbose=F)

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

pdf(paste0(fig_dir, "qc_violin_plot.pdf"), width=10, height=10)
VlnPlot(combined,group.by='cancer_type',features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0,)
dev.off()

#harmony
data("harmony_combined_from_server")

meta=cbind(combined$orig.ident,combined$platform,combined$type,combined$treatment.group,combined$cancer_type,combined$sample)
meta=as.data.frame(meta)

colnames(meta)=c('orig.ident','platform','type','treatment.group','cancer_type','sample')
#combined <- AddMetaData(combined,metadata=as.character(meta[colnames(combined),]$platform),col.name = 'platform')
combined <- AddMetaData(combined,metadata = meta)

#dim 정하기
pdf(paste0(fig_dir, "pca_harmony.pdf"), width=10, height=10)
p1 <- DimHeatmap(combined,dims = 1:20, cells=500,balanced=T)
p1
dev.off()

combined <- JackStraw(combined, num.replicate = 100) 
combined <- ScoreJackStraw(combined, dims=1:20) #최ㅐ대 20
pdf(paste0(fig_dir, "pca_jack_harmony.pdf"), width=10, height=10)
p1 <- JackStrawPlot(combined,dims=1:20)
p1
dev.off()

pdf(paste0(fig_dir,'elbow_harmony.pdf'),width=10,height=10)
p1 <- ElbowPlot(combined)
p1
dev.off()#15/10쯤 하면 될 것 같은 느낌

combined <- RunHarmony(combined,'platform')#
#dim will change
combined <- RunUMAP(combined,reduction='harmony',dims=1:10)
#경고: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.RNA.harmony; see ?make.names for more details on syntax validity
#경고메시지(들): 
# 1: Quick-TRANSfer stage steps exceeded maximum (= 11087050) 
#2: Quick-TRANSfer stage steps exceeded maximum (= 11087050) 

combined <- FindNeighbors(combined,reduction='harmony',dims = 1:10) #각 세포마다 distance를 구해서 비슷한 cell type인지 아닌지 찾음
combined <- FindClusters(combined, resolution=0.05)
#combined@meta.data
pdf(paste0(fig_dir, "run_platrom_harmony_PLATFORM_neighbor_reol0.05.pdf"), width=20, height=10)
p1 <- DimPlot(combined, split.by = 'platform', reduction = 'umap',raster=F)
p1#그림이 안보이는 거지 그려짐 ! 저장됨 
dev.off()

pdf(paste0(fig_dir, "run_platrom_harmony_CANCER_neighbor.pdf"), width=30, height=10)
p1 <- DimPlot(combined, split.by = 'cancer_type', reduction = 'umap',raster=F)
p1#그림이 안보이는 거지 그려짐 ! 저장됨 
dev.off()

pdf(paste0(fig_dir, "run_platrom_harmony_treatment_neighbor.pdf"), width=20, height=10)
p1 <- DimPlot(combined, split.by = 'treatment.group', reduction = 'umap',raster=F)
p1#그림이 안보이는 거지 그려짐 ! 저장됨 
dev.off()

saveRDS(combined, './data/harmony_combined_neighbor_0.05.rds')

#dim을 그냥 30으로 하고 했더니 cluster너무 많이 생기니까 이걸 조정할 필요있을 듯
#그리고 cancer 기준으로 halmony돌리고 platform봤는데 괜찮더라 

#platform 기준으로 하고 해도 된다!! find neighbor 함 -> 55개로 분리됨

