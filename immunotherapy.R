library(Seurat)
library(SeuratObject)
setwd('/Users/juyoung/Desktop/singlecell_jy')
#################
####gse115978####
#melanoma - post data 만 사용함 
#################
gse115978_counts=read.csv(file='/Users/juyoung/Downloads/GSE115978/GSE115978_counts.csv', header = T,row.names=1, as.is=T)
gse115978_cell_annotations=read.csv(file='/Users/juyoung/Downloads/GSE115978/GSE115978_cell.annotations.csv',header = T,row.names=1, as.is=T)

dim(gse115978_counts)#23686 7186
gse115978<- CreateSeuratObject(counts= gse115978_counts, project = 'gse115978',min.cells = 3, min.features = 200)
#grep("APOBEC3A-B",rownames(gse115978),value=T) #feature name이 '_'는 가진다고 warning 뜸 -> "APOBEC3A_B" "C4B_2"  -> 알아서 바꿔줌  
gse115978 #22454 7186

length(rownames(gse115978_cell_annotations))#7186
length(colnames(gse115978))#7186
length(colnames(gse115978_counts))#7186

gse115978#22454 7186
table(colnames(gse115978)==rownames(gse115978_cell_annotations))
gse115978$treatment.group=gse115978_cell_annotations$treatment.group

#add####################
#sample
gse115978$sample = gse115978_cell_annotations$sample
#cancer_type
gse115978$cancer_type = rep('melanoma',ncol(gse115978))
########################

gse115978_fu=gse115978[,which(gse115978$treatment.group=='post.treatment')]
table(gse115978_fu$treatment.group)#all post.treatment 3556

#b/fu 
gse115978_fu$treatment.group=rep('FU',3556)
#p/t
gse115978_fu$type= rep('T',3556)

#final check
table(gse115978_fu$type)
table(gse115978_fu$treatment.group)
table(gse115978_fu$sample)
table(gse115978_fu$cancer_type)
#final save
#save(gse115978_fu, file='/Users/juyoung/Desktop/singlecell_jy/data/gse115978_fu.rda')
save(gse115978_fu, file='/Users/juyoung/Desktop/singlecell_jy/data/gse115978_fu_add.rda')

data(gse115978_fu)

#add###############
for ( i in unique(gse115978_fu$sample)){
  assign(i, gse115978_fu[,which(gse115978_fu$sample==i)])
}
object_list=list()
for (i in unique(gse115978_fu$sample)){
  object_list=append(object_list,get(i))
}

assign('melanoma', merge(object_list[[1]],y=object_list[2:16], add.cell.ids =unique(gse115978_fu$sample), project= 'melanoma'))

#cehck
VlnPlot(melanoma,features = 'nFeature_RNA') #여전히 얘네만의 orig.ident로 나옴
Idents(melanoma)=melanoma$sample #이코드로 바꿀 수 있음

##################


#################
####gse123813####
#################

gse123813_bcc_counts=read.csv(file='/Users/juyoung/Downloads/GSE123813_bcc_scRNA_counts.csv',header = T,row.names = 1,sep='\t')
gse123813_bcc_metadata=read.csv(file='/Users/juyoung/Downloads/GSE123813_bcc_all_metadata.txt',header = T,row.names = 1,sep='\t')
gse123813_bcc <-CreateSeuratObject(counts= gse123813_bcc_counts, project = 'gse123813_bcc',min.cells = 3, min.features = 200)

#

gse123813_scc_counts=read.csv(file='/Users/juyoung/Downloads/GSE123813_scc_scRNA_counts.txt',header = T,row.names = 1,sep='\t')
gse123813_scc_metadata=read.csv(file='/Users/juyoung/Downloads/GSE123813_scc_metadata.txt',header = T,row.names = 1,sep='\t')
gse123813_scc <-CreateSeuratObject(counts= gse123813_scc_counts, project = 'gse123813_scc',min.cells = 3, min.features = 200)

#p/t
gse123813_bcc$type=rep('T',ncol(gse123813_bcc))
gse123813_scc$type=rep('T',ncol(gse123813_scc))
#b/fu
length(colnames(gse123813_bcc))#53029 #뭐가 사라진거지?
length(colnames(gse123813_bcc_counts))#53030
length(rownames(gse123813_bcc_metadata))#53030


################################
#없어진 친구를 찾기위해 
a=colnames(gse123813_bcc_counts)
b=colnames(gse123813_bcc)

a=as.data.frame(a)
colnames(a)='col'
dim(a)#53030 1

b=as.data.frame(b)
colnames(b)='col'
dim(b)#53029

setdiff(a$col,b$col)#"bcc.su010.post.tcell_ATTACTCCATCGACGC" #min.cell이거에 의해서 걸러진 샘플같기도함 ->282라 아닌데.. 
setdiff(b$col,a$col)# 없음 

sum(gse123813_bcc_counts[,"bcc.su010.post.tcell_ATTACTCCATCGACGC"])#282
###############################
#임시로 일단 저거 없이 진행
###############################
tempor_gse123813_bcc_metadata=gse123813_bcc_metadata[-which(rownames(gse123813_bcc_metadata)=='bcc.su010.post.tcell_ATTACTCCATCGACGC'),]
gse123813_bcc$treatment.group=ifelse(tempor_gse123813_bcc_metadata$treatment=='post','FU','B')

gse123813_scc$treatment.group=ifelse(gse123813_scc_metadata$treatment=='post','FU','B')

#add####################
#sample
gse123813_bcc$sample = tempor_gse123813_bcc_metadata$patient
gse123813_scc$sample = gse123813_scc_metadata$patient
#cancer_type
gse123813_bcc$cancer_type = rep('bcc',ncol(gse123813_bcc))
gse123813_scc$cancer_type = rep('scc',ncol(gse123813_scc))

for ( i in unique(gse123813_bcc$sample)){
  assign(i, gse123813_bcc[,which(gse123813_bcc$sample==i)])
}

for ( i in unique(gse123813_scc$sample)){
  assign(i, gse123813_scc[,which(gse123813_scc$sample==i)])
}

object_list=list()
for (i in unique(gse123813_bcc$sample)){
  object_list=append(object_list,get(i))
}

assign('bcc', merge(object_list[[1]],y=object_list[2:11], add.cell.ids =unique(gse123813_bcc$sample), project= 'bcc'))

for (i in unique(gse123813_scc$sample)){
  object_list=append(object_list,get(i))
}

assign('scc', merge(object_list[[1]],y=object_list[2:4], add.cell.ids =unique(gse123813_scc$sample), project= 'scc'))

Idents(bcc)=bcc$sample
Idents(scc)=scc$sample

###############
###위암,간암###
###############

#server - /homw/juyoung/scRNAseq/gse11~fol_local 이걸로 seuratobject다 만듬 -> 다시 local로 옮김 
#seurat object는 각각 rda 파일로 만들어서 /dat/juyoung/scRNAseq에 있음 (그냥 ./for_local만들고 안으로 집어넣음)
#rds파일로 저장 
#하나씩 치기 귀찮아서 command line 서버에서 만든거 가져옴(r_command.txt)
gastric_1058BP=readRDS('./data/gastric_1058BP.rds')
gastric_1058BT=readRDS('./data/gastric_1058BT.rds')
gastric_1058FUP=readRDS('./data/gastric_1058FUP.rds')
gastric_1058FUT=readRDS('./data/gastric_1058FUT.rds')
gastric_3622BP=readRDS('./data/gastric_3622BP.rds')
gastric_3622BT=readRDS('./data/gastric_3622BT.rds')
gastric_3622FUP=readRDS('./data/gastric_3622FUP.rds')
gastric_3622FUT=readRDS('./data/gastric_3622FUT.rds')
gastric_39054637BP=readRDS('./data/gastric_39054637BP.rds')
gastric_39054637BT=readRDS('./data/gastric_39054637BT.rds')
gastric_39054637FUP=readRDS('./data/gastric_39054637FUP.rds')
gastric_39054637FUT=readRDS('./data/gastric_39054637FUT.rds')
gastric_40784154BP=readRDS('./data/gastric_40784154BP.rds')
gastric_40784154BT=readRDS('./data/gastric_40784154BT.rds')
gastric_40784154FUP=readRDS('./data/gastric_40784154FUP.rds')
gastric_40784154FUT=readRDS('./data/gastric_40784154FUT.rds')
gastric_7727BP=readRDS('./data/gastric_7727BP.rds')
gastric_7727BT=readRDS('./data/gastric_7727BT.rds')
gastric_7727FUP=readRDS('./data/gastric_7727FUP.rds')
gastric_7727FUT=readRDS('./data/gastric_7727FUT.rds')
liver_10665612BP=readRDS('./data/liver_10665612BP.rds')
liver_10665612FUP=readRDS('./data/liver_10665612FUP.rds')
liver_12375821BP=readRDS('./data/liver_12375821BP.rds')
liver_12375821FUP=readRDS('./data/liver_12375821FUP.rds')
liver_22448359BP=readRDS('./data/liver_22448359BP.rds')
liver_22448359FUP=readRDS('./data/liver_22448359FUP.rds')
liver_32982982BP=readRDS('./data/liver_32982982BP.rds')
liver_32982982FUP=readRDS('./data/liver_32982982FUP.rds')
liver_36335616BP=readRDS('./data/liver_36335616BP.rds')
liver_36335616FUP=readRDS('./data/liver_36335616FUP.rds')
liver_38949491BP=readRDS('./data/liver_38949491BP.rds')
liver_38949491FUP=readRDS('./data/liver_38949491FUP.rds')
liver_39745926BP=readRDS('./data/liver_39745926BP.rds')
liver_39745926FUP=readRDS('./data/liver_39745926FUP.rds')
liver_40183243BP=readRDS('./data/liver_40183243BP.rds')
liver_40183243FUP=readRDS('./data/liver_40183243FUP.rds')
liver_40218172BP=readRDS('./data/liver_40218172BP.rds')
liver_40218172FUP=readRDS('./data/liver_40218172FUP.rds')
liver_41226138BP=readRDS('./data/liver_41226138BP.rds')
liver_41226138FUP=readRDS('./data/liver_41226138FUP.rds')

########################
gastric_name=list('gastric_1058BP','gastric_1058BT','gastric_1058FUP','gastric_1058FUT','gastric_3622BP','gastric_3622BT','gastric_3622FUP','gastric_3622FUT','gastric_39054637BP','gastric_39054637BT',
                  'gastric_39054637FUP','gastric_39054637FUT','gastric_40784154BP','gastric_40784154BT','gastric_40784154FUP','gastric_40784154FUT','gastric_7727BP','gastric_7727BT','gastric_7727FUP','gastric_7727FUT')
object_list=list()
for (i in gastric_name){
  object_list=append(object_list,get(i))
}

sample_list=list('1058BP','1058BT','1058FUP','1058FUT','3622BP','3622BT','3622FUP','3622FUT', '39054637BP','39054637BT','39054637FUP','39054637FUT','40784154BP','40784154BT', '40784154FUP','40784154FUT', '7727BP','7727BT','7727FUP','7727FUT')
#length(sample_list) 
assign('gastric', merge(object_list[[1]],y=object_list[2:20], add.cell.ids =sample_list, project= 'gastric'))
########################
liver_name=list('liver_10665612BP',
                'liver_10665612FUP',
                'liver_12375821BP',
                'liver_12375821FUP',
                'liver_22448359BP',
                'liver_22448359FUP',
                'liver_32982982BP',
                'liver_32982982FUP',
                'liver_36335616BP',
                'liver_36335616FUP',
                'liver_38949491BP',
                'liver_38949491FUP',
                'liver_39745926BP',
                'liver_39745926FUP',
                'liver_40183243BP',
                'liver_40183243FUP',
                'liver_40218172BP',
                'liver_40218172FUP',
                'liver_41226138BP',
                'liver_41226138FUP')
object_list=list()
for (i in liver_name){
  object_list=append(object_list,get(i))
}

sample_list=list('10665612BP',
                 '10665612FUP',
                 '12375821BP',
                 '12375821FUP',
                 '22448359BP',
                 '22448359FUP',
                 '32982982BP',
                 '32982982FUP',
                 '36335616BP',
                 '36335616FUP',
                 '38949491BP',
                 '38949491FUP',
                 '39745926BP',
                 '39745926FUP',
                 '40183243BP',
                 '40183243FUP',
                 '40218172BP',
                 '40218172FUP',
                 '41226138BP',
                 '41226138FUP')
assign('liver', merge(object_list[[1]],y=object_list[2:20], add.cell.ids =sample_list, project= 'liver'))


###########
###batch###
###########
#1.firm
#install.packages("devtools")
#install.packages('pbkrtest')
library(pbkrtest)
library(devtools)
#library(FIRM)

renites::install_version()
devtools::has_devel()
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
install_github("mingjingsi/FIRM")

#안됨 
#window는 됨
###########

#2.seurat integration
#10x와 smartseq 우선은 platform끼리 merge하지말고 구분할 컬럼을 만들어서 해보자
gastric$platform=rep('10x',ncol(gastric))
liver$platform=rep('10x',ncol(liver))
melanoma$platform=rep('Smartseq',ncol(melanoma)) #ㅎgse115978_fu를 했어야 함 
bcc$platform=rep('10x',ncol(bcc))
scc$platform=rep('10x',ncol(scc))

#test용이라 생각하고 복사하자
g1=gastric
l1=liver
#bcc=bcc
#scc=scc
mel=melanoma

#save(mel,file='./data/mel.rda')
#save(bcc,file='./data/bcc.rda')
#save(scc,file='./data/scc.rda')
#####again########
setwd('/Users/juyoung/Desktop/singlecell_jy')
data(mel)
data(g1)
data(bcc)
data(scc)
data("l1")

mel#3556
g1#200781
scc#14504
bcc#92869
l1#37679


library(stringr)

g1$type=ifelse(str_detect(g1$orig.ident,'P'),'P','T')
l1$type=ifelse(str_detect(l1$orig.ident,'P'),'P','T')

g1$treatment.group=(ifelse(str_detect(g1$orig.ident,'B'),'B','FU'))
l1$treatment.group=(ifelse(str_detect(l1$orig.ident,'B'),'B','FU'))

g1$cancer_type='gastric'
l1$cancer_type='liver'

g1_sample=c('1058','3622','39054637','40784154','7727')

sample=list()
for ( i in 1:ncol(g1)){
  sample=append(sample,print(gsub('\\D','',g1$orig.ident[i])))
}

g1[,which(grep(g1$orig.ident %in% i)]

g1$sample=g1$orig.ident
l1$sample=l1$orig.ident

g1$orig.ident

#MT 제거
#mt not remove

PercentageFeatureSet(mel, pattern = "^MT")
g1[['percent.mt']] <- PercentageFeatureSet(g1,pattern='^MT-')
l1[['percent.mt']] <- PercentageFeatureSet(l1,pattern='^MT-')
bcc[['percent.mt']] <- PercentageFeatureSet(bcc,pattern='^MT-')
scc[['percent.mt']] <- PercentageFeatureSet(scc,pattern='^MT-')
mel[['percent.mt']] <- PercentageFeatureSet(mel,pattern='^MT-')

mean(g1$percent.mt) #18.63595
max(g1$percent.mt) #96.9913
min(g1$percent.mt) #0

g1[,which(g1$percent.mt<5)] #200781 -> 16125
g1[,which(g1$percent.mt<10)] #200781 -> 98212
g1[,which(g1$percent.mt<15)] #200781 -> 129705 #10이나 15나 비슷하니까 sample살려서 15로 감 

mean(l1$percent.mt) #10.13572
max(l1$percent.mt) #95.67324
min(l1$percent.mt) #0

l1[,which(l1$percent.mt<5)] #37679 -> 12646
l1[,which(l1$percent.mt<10)] #37679 -> 27783 #이걸루 
l1[,which(l1$percent.mt<15)] #37679 -> 32163

mean(mel$percent.mt) #0
max(mel$percent.mt) #0
min(mel$percent.mt) #0

table(mel$percent.mt) #모두 0 

mean(bcc$percent.mt) #3.831172
max(bcc$percent.mt) #9.991197
min(bcc$percent.mt) #0

bcc[,which(bcc$percent.mt<5)] #53029 ->40970

mean(scc$percent.mt) #3.479254
max(scc$percent.mt) #9.991776
min(scc$percent.mt) #0

scc[,which(scc$percent.mt<5)] #26016 ->21841

#
l1 <- subset(l1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
g1 <- subset(g1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)
bcc <- subset(bcc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
scc <- subset(scc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
mel <-subset(mel,subset = nFeature_RNA > 200 & nFeature_RNA <2500) #mt는 다 0이니까 이것만 진행함 

# Visualize QC metrics as a violin plot
VlnPlot(g1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

l1=readRDS('./l1mt.rds')
g1=readRDS('./g1mt.rds')
bcc=readRDS('./bccmt.rds')
scc=readRDS('./sccmt.rds')
mel=readRDS('./melmt.rds')

getwd()
#normalize, find varable features
seurat_intergrate_list=list(g1,l1,bcc,scc,mel)
seurat_intergrate_list <- lapply(X=seurat_intergrate_list, FUN=function(x){
  x <-NormalizeData(x)
  x <-FindVariableFeatures(x,selection.method='vst',nfeatures=2000)
})

features <- SelectIntegrationFeatures(object.list = seurat_intergrate_list)

#anchors <- FindIntegrationAnchors(object.list=seurat_intergrate_list,anchor.features =features)
#combined <- IntegrateData(anchorset = anchors)
#2시간넘게 도는 중 -> 서버로옮김
setwd("/Users/juyoung/Desktop/singlecell_jy/data")
combined=readRDS('./combined.rds')#207657

DefaultAssay(combined) <-'integrated'
#defaultassay를 integrated로 설정하거나 rna로 할 수 있음
#여기서 integrated는 합쳐진 value를 가지고 분석을 한다는 거소
#rna는 original value에서 분석을 한다는 거다 
#https://www.biostars.org/p/399789/
combined <-ScaleData(combined, verbose=F)
combined <-RunPCA(combined, npcs=30, verbose=F)
combined <-RunUMAP(combined, reduction = 'pca',dims=1:30)
combined <-FindNeighbors(combined, reduction = 'pca',dims=1:30)
combined <-FindClusters(combined, resolution = 0.1)

combined$seurat_clusters
p1 <- DimPlot(combined,reduction='umap',group.by ='seurat_clusters',raster=F,size=0.05)
p2 <- DimPlot(combined, reduction = 'umap',label=T,repel=T,raster=F)
p2
p1

table(combined$orig.ident)
DimPlot(combined, reduction = 'type', split.by='')

#identify conserved cell type markers
DefaultAssay(combined) <-'RNA'
#BiocManager::install('multtest')
#install.packages('metap')
library(multtest)
library(metap)
options(future.globals.maxSize=8000*1024^2)
marker_1 <-FindConservedMarkers(combined, ident.1=1, grouping.var = 'platform',verbose = F)
marker_1 #믿을만한 값이 안나옴 



##########



###########

universe <-rownames(gse115978)
for (gastric in sample_list_gastric){
  universe=intersect(universe,rownames(get(gastric)))
}
for (liver in sample_list_liver){
  universe=intersect(universe,rownames(get(liver)))
}
universe=intersect(universe, rownames(gse123813))
length(universe) #5747

##########
##########
all_list = append(sample_list_gastric,sample_list_liver)
all_list[[41]]='gse115978'
all_list[[42]]='gse123813'

#subset
for ( i in all_list){
  assign(i,get(i)[universe,])
}

#누가 유전자 개수 줄이는 데 크게 기여했는지 보려고 
for (i in all_list){
  print(length(rownames(get(i))))
}

library(devtools)
devtools::install_github('satijalab/seurat-data',force=T)
devtools::install_github('satijalab/seurat-wrappers')
library(ggplot2)
library(patchwork)
library(SeuratWrappers)







