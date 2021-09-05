#https://blog.naver.com/study_very_hard/222061641671
#UMI BARCODE PCR ERROR 설명

library(dplyr)
#install.packages('Seurat')
library(Matrix) #이거 버전이 안맞아서 처음에 안됐음
#Matrix도 다시 설치해줘야했다 
library(Seurat)
library(patchwork) #그래프 그리는데 필요함 
setwd('/Users/juyoung/Desktop/singlecell_jy')
pbmc.data <-Read10X(data.dir='./filtered_gene_bc_matrices/hg19/')
#read10x 는  umi count matrix반환 
#각 cell 에서 detection 되는 feaure수 (gene 수)
# . 으로 표시되는 건 molecule이 발견되지 x 

#seurat object load(non-normalized data)
pbmc <- CreateSeuratObject(counts= pbmc.data, project = 'pbmc3k',min.cells = 3, min.features = 200)
#counts: raw data 지정
#projects: 내가 지정
#min cells : 적어도 특정 유전자를 발현한 세포가 이정도는 있어야한다.
#min feaures:적어도 이 정도 갯수의 유전자는 발현하고 있어야한다. 

# quality 안좋은 세포 거르기
pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc, patter= '^MT-') #^는 시작을 뜻함 $는 끝을 나타냄, MT로 시작하는 걸 걸러라 

#세포가 사멸하거나 세포막이 터지면 세포가 원래 발현하던건 밖으로 나가버림, 반면 미토콘드리아는 그 안에 있는 요소니까 상대적으로 발현이 안정됨/ 그래서 미토콘드리아 발현 비율이 높아 -> 세포상태 별로인것

#visualize qc metrics as a violin plot
VlnPlot(pbmc, features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol=3, cols='dodgerblue')

#nFeature_RNA - 각 Cell에서 감지된 유전자수(그림에서 보여지는 하나의 점이 하나의 셀임!) - mapping된 unique한 유전자 수 
#높다 -> double / multiplet 됐을 수 있음 
#낮다 -> dead/ dying/ empty droplet

#nCount_RNA - 각 CELL에서 몇개의 유전자가 mapping 됐는지 -> depth를 보는 느낌 ! 그래서 feature가 높으면 count가 높을테니까  
##molecule의 정확한 의미가 뭐지? rna말고 뭘 말하는 걸까 
#높다 -> double or multiplet


#각각의 violin plot으로 표현된 그래프를 합쳐서 2차원 그래프로 서로의 연관성을 나타낼 수 있음 

plot1 <- FeatureScatter(pbmc,feature1 = 'nCount_RNA',
                        feature2 = 'percent.mt')
plot2 <- FeatureScatter(pbmc, feature1 = 'nCount_RNA',
                        feature2 = 'nFeature_RNA') #얘는 높을 수 밖에 없지 ㅇ상관관계까 0.95인가 그정도임 
plot1+plot2

#그래프를 보면 nFeature 수가 지나치게 높거나 낮은 애들이 있음 + 미토콘드리아 도 -> 그래서 이거 제거 해야함

pbmc <- subset(pbmc,
               subset= nFeature_RNA >200 & nFeature_RNA < 2500 & percent.mt <5)

#절대적인 숫자는 아니고 분석에 따라 바꿀 수 있음
#교수님음 10%이하로 한 적도 있음 
#상태가 안좋으면 여기서 많이 걸러짐! 

#feature count는 그대로 쓰고 있음 
# percent.mt는 5%이하인 애들만 pbmc로 다시 지정

#normalizing할 차례임(두번 할거임)
#먼저, 각 샘플 마다 total read 기반으로 해주는 것 
#sample a mapping된 수 - 50/ sample b - 30//  total read가 앞이 1000이면 뒤는 100이야 그럼 b가 더 높은 거지! 

#pbmc  <- NormalizeData(pbmc, normalization.method = 'LogNormalize', scale.factor = 10000)
#평균적으로 total_read가 10000 정도라서 이걸로 나누는 것

pbmc <- NormalizeData(pbmc) #위에가 default값이라 그냥 이것만 해줘도됨

#cell간 variation이 높은 데이터를 찾음
#어떤 cell type인지 알기 위해 marker들을 찾으려고 하는 것 
pbmc <- FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)

#모든 세포를 기준으로 2000개를 뽑는 것 -> 만약 여기서 variable~가 없으면 전체적으로 똑같은 cell일테니까 전체를 기준으로 확인해보는 것
#이것만 봐도 어떤 세포type들이 있는지 확인이 가능함 ! 

#selection.method -> 분산이 큰 걸 뽑는 방법
#vst: 첫째, 로컬 다항식 회귀 분석(느림)을 사용하여 로그(분산)와 로그(평균)의 관계에 선을 적합시킵니다. 그런 다음 관측된 평균과 기대 분산을 사용하여 형상 값을 표준화합니다(적합선에 의해 제공됨). 그런 다음 최대치로 클리핑한 후 표준화된 값에 대해 형상 분산을 계산합니다(clip.max 매개 변수 참조 

#nfeatures : 몇개뽑을 건지 (selection 방법이 dispersion / vst일때만 사용함)
#그 cell의 특징을 나타내는 유전자들로 이 cell이 뭔지 찾아 낼 거니까 분산이 큰 유전자들이 뭔지 확인해야 알 수 있겠지

#
top10 <-head(VariableFeatureså(pbmc),10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot=plot1 , points=top10, repel=F)
plot1+ plot2

#정규화 작업, 발현평균은 0으로 분산은 1로 맞추는 과정
#scaling! 다음 작업을 위한 작업! 

# c1 c2 c3
#a 100 200 300
#b 1   7   10
#c 
#b도 b로만 보면의미가 있는데 a가 너무 크다보니 b가 묻힘 
#그래서 유전자별로 분포가 흡사하게 만들어줌 !!! -> 그래서 정규화 ㄱ작업이라 하는듯 

head(pbmc)
all.genes <- rownames(pbmc)
all.genes

#꼭 전체 데이터를 쓸 필요 없고 뽑은 2000개로 해도됨 
pbmc <- ScaleData(pbmc, features = all.genes)

#scale이 끝난 데이터로 pca 시작

pbmc <-RunPCA(pbmc, features = VariableFeatures(object=pbmc))
#고유값 고유벡터의 의미와 계산법
#https://angeloyeo.github.io/2019/07/17/eigen_vector.html
#pca결과 볼 수 있는 여러방법
#https://angeloyeo.github.io/2019/07/27/PCA.html 
# pca가 뭔지 이해하기 
#https://stats.stackexchange.com/questions/222/what-are-principal-component-scores
#pc score가 뭔지 이해하기
#https://www.graphpad.com/guides/prism/latest/statistics/stat_pca_example_pc_scores.htm
#1. 프린트하기
#dims 는 그 pca 그래프가 총 50개 있는데 그중에 1-5까지 본다는 말
print(pbmc[['pca']],dims=1:5, nfeatures= 5)
#pc1, pc2 이게 축인데! 이게 결국 어떤 cell type의 특성을 잘 보여주고 있을지를 확인해보는 거지 
#positive : pc1 축과 상관관계가 + 인 애들은 macrophage특징이 보이네 
#negative : 반대로 얘네는 macrophage랑은 ㄹㅇ 관련이 없구나 ~ 

#2. vizdimloadings - visualize dimensional reduction genes
VizDimLoadings(pbmc, dims=1:2, reduction = 'pca')
#x축은 고유값(공분산 행렬을 svd 했을때 나오는)-이게 결국 상관계수를 말함, 이때 이 고유값이 클수록 분산이 크다는 의미가 되기 때문에 그래프에서 오른쪽으로 갈 수록 유의미한 유전자가 됨  
#고유값은 고유벡터의 방향으로 얼만큼 퍼져있는지를 나타냄
#pc1이라는 벡터와 특정 유전자와의 관계성에 대한 특정값을 본 거다라고 생각해도된데
#숫자는 정확하게 뭔지 모르겠대
#절대값이 큰 순서로 보여준 것 같음 

#3.dimplot - dimentional reuction plot
DimPlot(pbmc, reduction='pca') 
#각각은 CELL을 나타내고 SCATTERPLOT임
#어떤 세포인지 모르는 상태에서 어떤 세포끼리 묶이는 지! 대략적으로 볼 수 있지 

#
DimHeatmap(pbmc, dims=1, cells = 500, balanced = T)
#TOP 500 CELL, +,- SCORE다 뽑겠다 -> BALANCED t
#pc1으로 했을 때 +관계 -관계에 있는 top 500개 cell -> 잘 나눠지더라 ~

DimHeatmap(pbmc, dims=1:15, cells=500, balanced = T)
##################
#몇개의 pc축을 가지고 두번째 dimension을 할거냐! 
#jackstrqw를 통해 우리가 가지고 있는 데이터들을 resampling해준다. 총 데이터에서 랜덤으로 1%씩 묶어 pca를 돌려 차원을 축소해줌 , 귀무가설 이런걸 함

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <-ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims=1:15)
#교수님 해석 -> 점선 위로 갈 수록 관계를 잘 표현한대! 근데 이것만 가지고 보기 어려워서 elbowplot을 봄 
#pvalue값이 기울기임, pc 10-12가 넘어가는 순간 pcalue가 높아짐 
#jackstrawplot은 qq plot을 그린것 
#qq plot은 이 데이터분포가 정규분포를 따르는지 안따르는지를 판단하는 것을 시각적으로 나타낸것/ x축은 데이터의 사분위수, y축은 이론상의 사분위수 이때 직선이 귀무가설이 맞다는 그런 선인가봐 그래서 이 선이랑 데이터 분포가 비슷하면 정규성을 만족한는 것이다. 
#https://bioinformaticsandme.tistory.com/37
ElbowPlot(pbmc) 
#pc에 대한 표편을 보여줌 
#pc 9-10에서 팔꿈치 같이 푹 꺼지는 경사를 볼 수 있음
#majority of true signal is captured in first 10 pc
#한pc9-10정도 까지 잡고 분석을 하면 양질의 분석을 얻을 수 있다는 뜻으로 이해하면 된다. 

#위에서 정한 dimension을 가지고 cell을 성질에 따라 분류함
pbmc <- FindNeighbors(pbmc, dims=1:10) #각 세포마다 distance를 구해서 비슷한 cell type인지 아닌지 
pbmc <- FindClusters(pbmc, resolution = 0.5)#그걸 기준으로 clustering을 해줌 

#resolution -> cluster의 개수를 얼마로 쪼개겠냐(5개/10개)/ cell type이 많을 수록 많이 ㅈ쪼개야하고 그래서 resolution값도 ㄷ올려야한다
#0.4-1.2가 괜찮은 range라 언급돼ㅣ있음 
#이 두함수는 성질이 비슷한 세포들을 한데 묶는 역할을 함

#2번째 reduction
#3차원에서 보면 2d에서 멀어보였던게 가까울 수 있더라! 
#non lineal 기반의 관계성을 보여줌
#근데 왜 pca를 볼까? 가장 정확하고 효율적임 , 왜 순서를 이렇게 하느냐 -> non-linear를 먼저하면 너무 많은 feature들의 loss가 일어날 수 있음 (수학적 공식에 의함 )
#umap 처리 후 dimplot으로 시각화하면 됨
#
pbmc <-RunUMAP(pbmc, dims=1:10)
DimPlot(pbmc, reduction='umap')

#tsne 세포내에서 비슷하더라 , cluster끼리 가까우니까 서로 연관이 있다는 말을 못함 
#umap cluster와 cluster가 ㄱ가까울수록 ㅇ비슷하다
#그래서 umap을 더 많이씀 , 그렇다고 tsne 쓰면 안되는 것도 아님! 

cluster1.markers <-FindMarkers(pbmc, ident.1=1, min.pct = 0.25)
#각 cluster가 무슨 cell인지 모르니까 enrich돼 있는 marker가 뭔지 찾는 것 
#cluster1 번의 25% cell에서 이러한 유전자들이 높더라
head(cluster1.markers, n=5)

#marker 찾는 기법 : scpred? -> annotation이 정확해(but 내가 가진 sample과 동일한 ref가 있어야함)
#sccatch -> ref sample이 없더라도 marker database로 찾아줌
#교수님은 sccatch 찾고 없으면 직접 찾는데..
 
############################
#find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25) #ident.1은 어떤 cluster를 볼건지
##min.pct는 두 그룹에서 나온 cell의 최소비

cluster2.markers

#cluster0과 3이랑 다른 5의 특징을 나타낼 수 있는 마커
cluster5.markers <- FindMarkers(pbmc, ident.1=5, ident.2=c(0,3), min.pct=0.25)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

library(dplyr)
pbmc.markers %>% group_by(cluster) %>% top_n(n=2,wt=avg_log2FC)

#top_n : top/bottom n개의 row를 뽑는 것 , wt= ordering에 사용되는 기준 

cluster0.markers <- FindMarkers(pbmc, ident.1=0,
                                logfc.threshold = 0.25,
                                test.use = 'roc',
                                only.pos = T)
#여기서 사용된 roc는 개별마커에 대한 분류력을 반환함 
#0은 무작위, 1은 perfect
VlnPlot(pbmc, features = c('MS4A1','CD79A'))

VlnPlot(pbmc, features = c('NKG7','PF4'),slot='counts', log=T) # raw count data를 사용할 수 있음 

FeaturePlot(pbmc, features = c('MS4A1','GNLY','CD3E','CD14','FCER1A'))

#상위 20개 MARKER/ 만일 마커가 20개이하면 그거 전체를 PLOTTING함
top10 <-pbmc.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) +NoLegend()
#NOLEGEND 이거는 범주 범위 표시하는 거 없앰

#이 데이터 셋은 표준 마커를 사용하여 편향되지 않은 클러스터링을 알려진 셀 유형에 쉽게 일치 시킬 수 있다. 
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <-levels(pbmc)
levels(pbmc) #0부터 8까지 나옴 cluster 개수

new.cluster.ids #0-navie cd4 t,이런식으로 됨

pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction='umap',label=T, pt.size=0.5) +NoLegend()



