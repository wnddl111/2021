data("f_merge_gene")
#install.packages('BiocManager')
#BiocManager::install('sva')
#BiocManager::install('gridExtra')
#BiocManager::install('UpSetR')
#BiocManager::install('edgeR')
#install.packages('ggplot2')
library(gridExtra)
library(edgeR)
library(ggplot2)
library(sva)
library(UpSetR)
'''
condition =c(rep('T',36),rep('N',948))
#dim(final_merge)
f <- data.frame(sapply(final_merge[2:985], function(x) as.numeric(as.character(x)))) 
dim(f)#54398 984
f1 <- as.data.frame(final_merge[,1])
f_merge=cbind(f1,f)
head(f_merge)
colnames(f_merge)=colnames(final_merge)
'''
pca_uncorrected_obj = prcomp(f_merge_gene[,37:984])#pca한것
head(pca_uncorrected_obj)
pca_uncorrected = as.data.frame(pca_uncorrected_obj[2]$rotation) #pca value 뽑음
pca_uncorrected[,'condition']=condition

cols <- c("N" = "#1F968BFF")
p1 = ggplot(data=pca_uncorrected, aes(x=PC1, y=PC2))
p1 = p1 + geom_point(size=3)
p1 = p1 + stat_ellipse(type="norm", linetype=2)
p1 = p1 + labs(title="PCA, RNA-seq N samples (uncorrected data)")
p1 = p1 + scale_colour_manual(values = cols)
p1

c=pca_uncorrected[pca_uncorrected$PC1<0.06 & pca_uncorrected$PC2<0.1,]
p2 = ggplot(data=c, aes(x=PC1, y=PC2))
p2 = p2 + geom_point(size=3)
p2 = p2 + stat_ellipse(type="norm", linetype=2)
p2 = p2 + labs(title="PCA, RNA-seq N samples (1st corrected data)")
p2 = p2 + scale_colour_manual(values = cols)
p2

c2=c[c$PC1<0.04 & c$PC2<0,]
p3 = ggplot(data=c2, aes(x=PC1, y=PC2))
p3 = p3 + geom_point(size=3)
p3 = p3 + stat_ellipse(type="norm", linetype=2)
p3 = p3 + labs(title="PCA, RNA-seq N samples (2nd corrected data)")
p3 = p3 + scale_colour_manual(values = cols)
p3

rownames(c2)#577

c3=c2[c2$PC1<0.035 & c2$PC2< -0.01,]
p4 = ggplot(data=c3, aes(x=PC1, y=PC2))
p4 = p4 + geom_point(size=3)
p4 = p4 + stat_ellipse(type="norm", linetype=2)
p4 = p4 + labs(title="PCA, RNA-seq N samples (3rd corrected data)")
p4 = p4 + scale_colour_manual(values = cols)
p4

rownames(c3) #454

set.seed(777)
sample=sample(rownames(c3),54)#

#c_f로 하니까 또 너무 모여있는 것 같아서 1번째 correction한 곳에서 random하게 뽑는 걸로 수정 
#아냐 교수님은 별 상관ㅇ ㅇ없을 것 같대 그리고 얘네만 뽑아서 PCA 그려봤는데 나름 모이면서 퍼져있달까? 괜찮은듯!

c_f=c3[sample,]#최종
p5= ggplot(data=c_f, aes(x=PC1,y=PC2))
p5 = p5 + geom_point(size=3)
p5 = p5 + stat_ellipse(type='norm', linetype=2)
p5 = p5 + labs(title='PCA, RNA-seq N samples (final selected data)')
p5 
'''
set.seed(777)
sample=sample(rownames(c),56)#
c_f=c[sample,]
'''
type=ifelse(rownames(pca_uncorrected) %in% rownames(c_f),'selected','not selected')
pca_uncorrected[,'type']=type
table(type)#54 894

cols <- c("not selected" = "#481567FF", "selected" = "#1F968BFF")
p6 = ggplot(data=pca_uncorrected, aes(x=PC1, y=PC2, color=type))
p6 = p6 + geom_point(size=3)
p6= p6 + stat_ellipse(type="norm", linetype=2)
p6 = p6 + labs(title="PCA, RNA-seq counts for T/N(GTEX) samples", color="type")
p6 = p6 + scale_colour_manual(values = cols)
p6

pca_corrected_obj = prcomp(f_merge_gene[,sample])#
pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation) #pca value 뽑음
cols <- c("N" = "#1F968BFF")
p7 = ggplot(data=pca_corrected, aes(x=PC1, y=PC2))
p7 = p7 + geom_point(size=3)
p7 = p7 + stat_ellipse(type="norm", linetype=2)
p7 = p7 + labs(title="PCA, RNA-seq N samples (corrected data)")
p7 = p7 + scale_colour_manual(values = cols)
p7


save(sample,file='./data/GTEX_54_selected.rda')

