setwd('D:/')
library(corrplot)
vst = read.table(file='./data/F0405_for_heatmap_0.3.txt', sep='\t')
col = read.table(file='./data2/my_data_colname_90.txt')

head(vst)
gene = vst$V1
rownames(vst) = gene
#vst = vst[,-1]
colnames(vst)= col$V1

#튀는 normal 제거
library(pheatmap)
data("tcga_cosmic_F_vst_o")
pheatmap(F_vst_o[,37:90], cluster_rows = F, cluster_cols = F) #u3gz 환자만 튐 61번
#vst = vst[,-61]
colnames(vst)

#crpc vs n 
r_sample=colnames(vst)[1:36] 
n_sample=colnames(vst)[37:89]

r_vst = vst[,r_sample]
n_vst = vst[,n_sample]

#ar gene에 대한 분포 확인
library(ggplot2)

ar_vst = as.data.frame(t(vst['AR',]))
ar_vst$sample = c(rep('T',36),rep('N',53))
#z=as.data.frame(colnames(AR_up_vst))
#write.csv(z, file = './data2/삭제2.csv')
head(ar_vst)
ar_vst %>%
  ggplot(aes(x=sample, y=AR, fill= sample))+
  geom_violin(trim=F)+
  geom_boxplot(width=0.1)+
  ggtitle('AR gene boxplot in CRPC, GTEx Normal')

M= mean(apply(vst['AR',],2,as.numeric))#0.008067849

AR_up_vst = vst[,which(vst['AR',]>=M)]
colnames(AR_up_vst)
AR_down_vst = vst[,which(vst['AR',] < M)] 
colnames(AR_down_vst)

#AR correlation 보기
AR_up_vst=t(AR_up_vst)
AR_down_vst = t(AR_down_vst)

write.table(AR_up_vst,'./data/AR_up_vst.txt',sep='\t')
write.table(AR_down_vst,'./data/AR_down_vst.txt',sep='\t')

#sample별로 AR 분포 보기 

z=as.data.frame(t(vst['AR',r_sample]))
z
z$sample = rownames(z)
p <- ggplot(z, aes(y=AR, x=sample))+
  geom_point()
p <- ggMarginal(p, margins='y', color='purple',size=4)
p
###################시도1
#crpc vs n 
r_sample=colnames(vst)[1:36] 
n_sample=colnames(vst)[37:90]

#r sample에서 ar이 높은 그룹 낮은 그룹으로 표현
r_vst = vst[,r_sample]
colnames(r_vst)

r_vst = apply(r_vst,2, as.numeric)
rownames(r_vst)= rownames(vst)
head(r_vst)
hist(r_vst['AR',])
M= mean(r_vst['AR',])#1.053007

AR_up_r_vst = r_vst[,which(r_vst['AR',]>=M)]
AR_down_r_vst = r_vst[,which(r_vst['AR',] < M)] 

#AR correlation 보기
AR_up_r_vst=t(AR_up_r_vst)
AR_down_r_vst = t(AR_down_r_vst)

write.table(AR_up_r_vst,'./data/AR_up_r_vst.txt',sep='\t')
write.table(AR_down_r_vst,'./data/AR_down_r_vst.txt',sep='\t')


