library(Seurat)
library(SeuratObject)
setwd('/Users/juyoung/Desktop/immunotherapy_jy/singlecell_jy')
################
#local
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
#ident
Idents(gse115978_fu)= gse115978_fu$sample
#final check
table(gse115978_fu$type)
table(gse115978_fu$treatment.group)
table(gse115978_fu$sample)
table(gse115978_fu$cancer_type)

#final save
#save(gse115978_fu, file='/Users/juyoung/Desktop/singlecell_jy/data/gse115978_fu.rda')
save(gse115978_fu, file='/Users/juyoung/Desktop/singlecell_jy/data/gse115978_fu_add.rda')
data("gse115978_fu_add")
#cehck
VlnPlot(gse115978_fu,features = 'nFeature_RNA') 

#################
#testgse123813###
#################
gse123813_bcc_counts=read.csv(file='/Users/juyoung/Downloads/GSE123813_bcc_scRNA_counts.csv',header = T,row.names = 1,sep='\t')
gse123813_bcc_metadata=read.csv(file='/Users/juyoung/Downloads/GSE123813_bcc_all_metadata.txt',header = T,row.names = 1,sep='\t')

dim(gse123813_bcc_metadata)#53030 7
dim(gse123813_bcc_counts)#23309 53030

gse123813_bcc_metadata$treatment.group=ifelse(gse123813_bcc_metadata$treatment=='post','FU','B')

bcc_seurat_name_list=list()
for (i in unique(gse123813_bcc_metadata$patient)){
  a=gse123813_bcc_metadata[which(gse123813_bcc_metadata==i),]
  for ( j in unique(gse123813_bcc_metadata$treatment.group)){
    assign(paste0(i,j,'T_'),rownames(a[which(a$treatment.group==j),]))
    bcc_seurat_name_list=append(bcc_seurat_name_list,paste0(i,j,'T_'))    
  }
}

#create seurat obj하기전 확인
#유전자개수는 변화 없으니까 cell개수만 봄
for (i in bcc_seurat_name_list){
  print(length(get(i)))
}

#sample별로 구분
for (i in bcc_seurat_name_list){
  assign(i,CreateSeuratObject(counts= gse123813_bcc_counts[,get(i)], project = i,min.cells = 3, min.features = 200))
}

#나눴을 때 loss정도를 확인해봄
#loss가 큰 것도 존재함 (약 35%)
for (i in bcc_seurat_name_list){
  print(ncol(GetAssay(get(i))))
}

#나눠진 걸로 merge를 해보자 어떻게 되려나? 
object_list=list()
for (i in bcc_seurat_name_list){
  object_list=append(object_list,get(i))
}

assign('bcc_use_split_count', merge(object_list[[1]],y=object_list[2:22], add.cell.ids =bcc_seurat_name_list, project= 'bcc_use_split_count'))
bcc_use_split_count#22928 53029 왜 이렇게 숫자가 커질 수 있는걸까? outerjoin

#따로 있는거랑 merge랑의 비교
#따로 있을 때 유전자는?
split_su004FUT_gene=rownames(su004FUT@assays$RNA@counts)
#merge한거에서 su004FUT에 해당하는 거 한번에 비교해볼랬는데 걍 한 cell로만 해봄 
test_su004fut_for_merge=list()
for (i in get('su004FUT_')){
  test_su004fut_for_merge=append(test_su004fut_for_merge,paste0('su004FUT_',i))  
}
#"su004FUT_bcc.su004.post.tumor.cd45_TTTATGCAGAACAACT" 이걸로해봄

test_cell_merge_count=as.data.frame(GetAssayData(object = bcc_use_split_count, slot = "counts")[,"su004FUT_bcc.su004.post.tumor.cd45_TTTATGCAGAACAACT"])
#여기서 split때는 없었던 유전자에 대해서는 무슨값을 가지나?
table(test_cell_merge_count[-which(rownames(test_cell_merge_count) %in% split_su004FUT_gene),]) #모두 0으로 표시하는 구나 7292
##그럼 통째로 하면 loss가 없나?
gse123813_bcc_use_rawcount <-CreateSeuratObject(counts= gse123813_bcc_counts, project = 'gse123813_bcc_raw',min.cells = 3, min.features = 200)
gse123813_bcc_use_rawcount #22961 53029 엄청 안줄긴하구낭.-. 

test_use_rawcount=as.data.frame(gse123813_bcc_use_rawcount@assays$RNA@data[,"bcc.su004.post.tumor.cd45_TTTATGCAGAACAACT"])
table(test_use_rawcount[-which(rownames(test_use_rawcount) %in% split_su004FUT_gene),])#0 7319 1 6 얘는 모듀 0이 아니굼 
#################
#last gse123813##
#################
#project별로 통째로 진행 
gse123813_bcc <-CreateSeuratObject(counts= gse123813_bcc_counts, project = 'gse123813_bcc',min.cells = 3, min.features = 200)
#
gse123813_scc_counts=read.csv(file='/Users/juyoung/Downloads/GSE123813_scc_scRNA_counts.txt',header = T,row.names = 1,sep='\t')
gse123813_scc_metadata=read.csv(file='/Users/juyoung/Downloads/GSE123813_scc_metadata.txt',header = T,row.names = 1,sep='\t')

dim(gse123813_scc_counts)#18347 26016
gse123813_scc <-CreateSeuratObject(counts= gse123813_scc_counts, project = 'gse123813_scc',min.cells = 3, min.features = 200)
#gse123813_scc 17559 26016
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

setdiff(a$col,b$col)#"bcc.su010.post.tcell_ATTACTCCATCGACGC" #min.cell이거에 의해서 걸러진 샘플같기도함 ->맞음
setdiff(b$col,a$col)# 없음 

nrow(gse123813_bcc_counts[which(gse123813_bcc_counts$"bcc.su010.post.tcell_ATTACTCCATCGACGC">0),])#195 
#여기 

###############################
#저거 없이 진행
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

Idents(gse123813_bcc)=gse123813_bcc$sample
Idents(gse123813_scc)=gse123813_scc$sample

###############
###위암,간암###
#server
###############

setwd('/data/juyoung/scRNAseq/for_local/') #얘는 위치가 여기임 
gastric_1058BP=readRDS('./gastric_1058BP.rds')
gastric_1058BT=readRDS('./gastric_1058BT.rds')
gastric_1058FUP=readRDS('./gastric_1058FUP.rds')
gastric_1058FUT=readRDS('./gastric_1058FUT.rds')
gastric_3622BP=readRDS('./gastric_3622BP.rds')
gastric_3622BT=readRDS('./gastric_3622BT.rds')
gastric_3622FUP=readRDS('./gastric_3622FUP.rds')
gastric_3622FUT=readRDS('./gastric_3622FUT.rds')
gastric_39054637BP=readRDS('./gastric_39054637BP.rds')
gastric_39054637BT=readRDS('./gastric_39054637BT.rds')
gastric_39054637FUP=readRDS('./gastric_39054637FUP.rds')
gastric_39054637FUT=readRDS('./gastric_39054637FUT.rds')
gastric_40784154BP=readRDS('./gastric_40784154BP.rds')
gastric_40784154BT=readRDS('./gastric_40784154BT.rds')
gastric_40784154FUP=readRDS('./gastric_40784154FUP.rds')
gastric_40784154FUT=readRDS('./gastric_40784154FUT.rds')
gastric_7727BP=readRDS('./gastric_7727BP.rds')
gastric_7727BT=readRDS('./gastric_7727BT.rds')
gastric_7727FUP=readRDS('./gastric_7727FUP.rds')
gastric_7727FUT=readRDS('./gastric_7727FUT.rds')


liver_10665612BP=readRDS('./liver_10665612BP.rds')
liver_10665612FUP=readRDS('./liver_10665612FUP.rds')
liver_12375821BP=readRDS('./liver_12375821BP.rds')
liver_12375821FUP=readRDS('./liver_12375821FUP.rds')
liver_22448359BP=readRDS('./liver_22448359BP.rds')
liver_22448359FUP=readRDS('./liver_22448359FUP.rds')
liver_32982982BP=readRDS('./liver_32982982BP.rds')
liver_32982982FUP=readRDS('./liver_32982982FUP.rds')
liver_36335616BP=readRDS('./liver_36335616BP.rds')
liver_36335616FUP=readRDS('./liver_36335616FUP.rds')
liver_38949491BP=readRDS('./liver_38949491BP.rds')
liver_38949491FUP=readRDS('./liver_38949491FUP.rds')
liver_39745926BP=readRDS('./liver_39745926BP.rds')
liver_39745926FUP=readRDS('./liver_39745926FUP.rds')
liver_40183243BP=readRDS('./liver_40183243BP.rds')
liver_40183243FUP=readRDS('./liver_40183243FUP.rds')
liver_40218172BP=readRDS('./liver_40218172BP.rds')
liver_40218172FUP=readRDS('./liver_40218172FUP.rds')
liver_41226138BP=readRDS('./liver_41226138BP.rds')
liver_41226138FUP=readRDS('./liver_41226138FUP.rds')
########################
#20220102 위암이랑 간암도 sample별로 할 때랑 project로 할 때 소실이 컸음... 그래서 다시 돌리고 진행함 
########################
gastric_name=list('gastric_1058BP','gastric_1058BT','gastric_1058FUP','gastric_1058FUT','gastric_3622BP','gastric_3622BT','gastric_3622FUP','gastric_3622FUT','gastric_39054637BP','gastric_39054637BT',
                  'gastric_39054637FUP','gastric_39054637FUT','gastric_40784154BP','gastric_40784154BT','gastric_40784154FUP','gastric_40784154FUT','gastric_7727BP','gastric_7727BT','gastric_7727FUP','gastric_7727FUT')
object_list=list()
for (i in gastric_name){
  object_list=append(object_list,get(i))
}

gene=rownames(gastric_1058BP)
for( i in object_list){
  print(table(gene==rownames(i)))
} #all true 예상대로 feature 모두 같음 

sample_list=list('1058BP','1058BT','1058FUP','1058FUT','3622BP','3622BT','3622FUP','3622FUT', '39054637BP','39054637BT','39054637FUP','39054637FUT','40784154BP','40784154BT', '40784154FUP','40784154FUT', '7727BP','7727BT','7727FUP','7727FUT')
#length(sample_list)) 
assign('gastric', merge(object_list[[1]],y=object_list[2:20], add.cell.ids =sample_list, project= 'gastric'))
gastric #36601 214946 
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

gene=rownames(liver_41226138FUP)
for( i in object_list){
  print(table(gene==rownames(i)))
} #all true 예상대로 feature 모두 같음 

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

g1<-CreateSeuratObject(counts=gastric@assays$RNA@counts, project='gastric', min.cells = 3, min.features = 200)
g1 #28050 200781 

l1<-CreateSeuratObject(counts=liver@assays$RNA@counts, project='liver', min.cells = 3, min.features = 200)
l1 #20496 37679

g1$platform=rep('10x',ncol(g1)) #g1
l1$platform=rep('10x',ncol(l1))#l1 
mel$platform=rep('Smartseq2',ncol(melanoma)) #ㅎgse115978_fu를 했어야 함 
bcc$platform=rep('10x',ncol(bcc))
scc$platform=rep('10x',ncol(scc))

#test용이라 생각하고 복사하자
bcc=gse123813_bcc
scc=gse123813_scc
mel= gse115978_fu

setwd('/data/juyoung/scRNAseq/for_local/immunotherapy_jy/singlecell_jy')
#save(mel,file='./data/mel.rda')
#save(bcc,file='./data/bcc.rda')
#save(scc,file='./data/scc.rda')
#save(g1,file='./data/g1.rda')
#save(l1,file='./data/l1.rda')

#####again########
data(mel)
data(g1)
data(bcc)
data(scc)
data("l1")

mel#3556
g1#200781
scc#26016
bcc#53029
l1#37679

library(stringr)

g1$type=ifelse(str_detect(g1$orig.ident,'P'),'P','T')
l1$type=ifelse(str_detect(l1$orig.ident,'P'),'P','T')

g1$treatment.group=(ifelse(str_detect(g1$orig.ident,'B'),'B','FU'))
l1$treatment.group=(ifelse(str_detect(l1$orig.ident,'B'),'B','FU'))

g1$cancer_type='gastric'
l1$cancer_type='liver'

Idents(g1)=gsub('\\D','',g1$orig.ident)
Idents(l1)=gsub('\\D','',l1$orig.ident)

g1$sample=gsub('\\D','',g1$orig.ident)
l1$sample=gsub('\\D','',l1$orig.ident)

#MT 제거
#mt not remove

grep('MT-',rownames(mel@assays$RNA@data))#3376 INMT-FAM188B 
mel@assays$RNA@data[3376:3377,]#MT gene이 애초에 없음 #미리 없애줬나봄 

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
l1[,which(l1$percent.mt<10)] #37679 -> 27783 
l1[,which(l1$percent.mt<15)] #37679 -> 32163#이걸루 

mean(mel$percent.mt) #0
max(mel$percent.mt) #0
min(mel$percent.mt) #0

table(mel$percent.mt) #모두 0 

mean(bcc$percent.mt) #3.831172
max(bcc$percent.mt) #9.991197
min(bcc$percent.mt) #0

bcc[,which(bcc$percent.mt<7)] #53029 -> 49519

mean(scc$percent.mt) #3.479254
max(scc$percent.mt) #9.991776
min(scc$percent.mt) #0

scc[,which(scc$percent.mt<5)] #26016 ->21841

#
l1 <- subset(l1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)
g1 <- subset(g1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)
bcc <- subset(bcc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 7)
scc <- subset(scc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#mel <-subset(mel,subset = nFeature_RNA > 200 & nFeature_RNA <2500) #mt는 다 0이니까 이것만 진행함 
#ss2는 doublet이 생기지 않는다? ㅇㅇ

# Visualize QC metrics as a violin plot
VlnPlot(g1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#saveRDS(l1,file='./l1_last.rds')
#saveRDS(g1,file='./g1_last.rds')
#saveRDS(bcc,file='./bcc_last.rds')
#saveRDS(scc,file='./scc_last.rds')
#saveRDS(mel,file='./mel_last.rds')
