library(Biobase)
data("tcga_prostate_gene_0222_eset")
count = exprs(tcga_prostate_gene_eset)

pdata= pData(tcga_prostate_gene_eset)
rownames(pdata)
dim(pdata)#551 147
head(pdata)

clinical=read.csv('./data/clinical.csv')
intersect(clinical$X_PATIENT, pdata$patient)#495
length(unique(clinical$X_PATIENT))#498
length(unique(pdata$patient))#495

pdata$sample#vial까지 표시
clinical$sampleID #vial표시 안함

a=pdata[pdata$sample=='TCGA-HC-8258-01A',]
b=pdata[pdata$sample=='TCGA-HC-8258-01B',]
all.equal(a,b)
#어차피 clinical에 있는 정보는 vial로 구분 할 수 있는게 아니니까 그대로 써야지 
library(stringr)
pdata_sampleID=sapply(pdata$sample,function(x)
  str_sub(x,1,-2))
pdata$clinial_sampleID=pdata_sampleID 

#clinical에서 pdata랑 일치하는 것만 가져오기
#clinical에만 이 정보가 온전히 있어서 그러는 것
rownames(clinical)=clinical$sampleID
dim(clinical[pdata_sampleID,])#551 105

length(clinical$sampleID)#566
table(rownames(clinical)==clinical$sampleID)#all true
length(unique(clinical$sampleID))#566
length(unique(pdata_sampleID))#548
length(pdata_sampleID)#551

pdata_sampleID=unique(pdata_sampleID)
final_pdata=clinical[pdata_sampleID,]
rownames(final_pdata)
dim(final_pdata)#548 105

exprs_clinical_col=lapply(colnames(count), function(x)
  str_sub(x,1,15)
  )
colnames(count)=exprs_clinical_col
unique(colnames(count))#548

final_count=count[,pdata_sampleID]#겹치는게 있으면 하나면 골라지구나..ㅎ.
###
count[,"TCGA-HC-8258-01"]
cbind(count[,256],count[,255])#vial의차이..
dim(final_count)#38662 548
##
#######eset#####

featureData=new('AnnotatedDataFrame',data=fData(tcga_prostate_gene_eset))
phenoData=new('AnnotatedDataFrame',data=final_pdata)
tcga_prostate_gene_clinical_eset=
  new('ExpressionSet', exprs=final_count, phenoData=phenoData,featureData=featureData)
save(tcga_prostate_gene_clinical_eset,file='./data/tcga_prostate_gene_clinical_0222_eset.rda')

#n,t,r구분
data("tcga_prostate_gene_clinical_0222_eset")
pdata= pData(tcga_prostate_gene_clinical_eset)
exp=exprs(tcga_prostate_gene_clinical_eset)
fdata=fData(tcga_prostate_gene_clinical_eset)

r_pdata <- pdata[pdata$biochemical_recurrence=='YES'& substr(pdata$sampleID,14,15)=='01',]

dim(r_pdata)#58 105
r_sample = rownames(r_pdata)

p_pdata <- pdata[pdata$biochemical_recurrence=='NO'& substr(pdata$sampleID,14,15)=='01',]

dim(p_pdata)#369 105
p_sample = rownames(p_pdata)

n_pdata <- pdata[substr(pdata$sampleID,14,15)=='11',]

dim(n_pdata) #52 105
n_sample = rownames(n_pdata)

##########################
f_pdata = pdata[c(n_sample,p_sample,r_sample),]
f_exp =exp[,c(n_sample,p_sample,r_sample)] 

#01로 끝나면서 recurrence정보 없는게 있어서 줄어듬 
dim(f_pdata)#479 105
dim(f_exp)#38662 479

f_pdata$t_vs_n_group=c(rep('n',52),rep('t',427))
f_pdata$n_vs_p_vs_r_group = c(rep('n',52),rep('p',369),rep('r',58))

##########################
library(DESeq2)
g=rownames(f_exp)
s=colnames(f_exp)

f=as.matrix(f_exp)
storage.mode(f)='integer'

f_pdata$t_vs_n_group
table(is.na(f))#모두 없음 

479*0.1
dim(f[rowSums(f==0) <= 48,])#38662 479->19574 479
dim(f)#38662 479


#안줄었네 모두 0인 유전자는 없다
#여기서 많이 줄이고 들어가는 사람도 있던데
#그래서 전체 sample의 반이상으로 0가지고 있으면 없앰
f=f[rowSums(f==0) <= 48,]
table(apply(f,1, is.infinite))#false 

dds <- DESeqDataSetFromMatrix(
  countData = f,
  colData=f_pdata,
  design=~t_vs_n_group, tidy=F
)

#tidy는 countdata의 첫번째 column이 rownames인지
dds
dds <- DESeq(dds)
n_vs_t_res <- results(dds)
head(results(dds, tidy = F))
n_vs_t_res = as.data.frame(n_vs_t_res)
write.csv(n_vs_t_res,file='./data/prad_zerocnt_0.1_res.csv')

library(dplyr)
#all_equal(n_vs_t_res, n_vs_t_res_1)#true
#dds에 sizefactor이런거 하고나서 하면 뭐가 바뀌나 했는데 ㄴㄴ 
n_vs_t_res <- n_vs_t_res[n_vs_t_res$pvalue<=0.05& abs(n_vs_t_res$log2FoldChange)>=1.5,]
dim(n_vs_t_res)#1082 6 10%
n_vs_t_res <- n_vs_t_res[order(n_vs_t_res$log2FoldChange),]
down_nt_res <- n_vs_t_res[n_vs_t_res$log2FoldChange<=-1.5,]
dim(down_nt_res)#847 6 -> 607 6

nt_bottom=rownames(down_nt_res)[1:25]
##############
n_vs_t_res <- n_vs_t_res[order(-n_vs_t_res$log2FoldChange),]
up_nt_res <- n_vs_t_res[n_vs_t_res$log2FoldChange>=1.5,]
head(up_nt_res)
dim(up_nt_res)#797 6 -> 475 6 
nt_top=rownames(up_nt_res)[1:25]

###################
#log2변환
save(dds, file = './data/nt_dds_0.1.rda') #0.5는 _0.1만 없이 저장 
save(n_vs_t_res,file = './data/n_vs_t_res_0.1.rda')


#vsd <-vst(dds, blind=F) #rna seq에 안쓴대 
#dds <- estimateSizeFactors(dds)
#norm.counts <- counts(dds, normalized=T)
#range(log.norm.counts)
#log.norm.counts <- log2(norm.counts+1)
#hist(log.norm.counts)

range(log.norm.counts)
range(assay(vsd))
mean(assay(vsd))
#
hist(counts(dds, normalized=T))
log2_count <- log2(counts(dds, normalized=T)+1)
range(log2_count)#0 23.97126 -> 0 22.04083
hist(log2_count)
scale_log2_count <- t(scale(t(log2_count)))
range(scale_log2_count)#-12.78 12.25 -> -12.78 12.25


hist(scale_log2_count)

library(pheatmap)
df <-as.data.frame(colData(dds)[,c('t_vs_n_group')])
rownames(df)=colnames(scale_log2_count)
colnames(df)='Type'

n_vs_t_res #top은 5432, down은 -10,-9,-8~
dim(n_vs_t_res)
pheatmap(scale_log2_count[c(nt_top,nt_bottom),], cluster_rows = F,show_rownames =T, show_colnames = F, cluster_cols = F, annotation_col = df,breaks=seq(-1,1, length.out = 100))

####################
#f가 0.1로 0뺀 exp
dim(f)
dim(f[,c(n_sample,p_sample)])#19574 421, 23837 - 0.5기준 
np_f=f[,c(n_sample,p_sample)]
nr_f=f[,c(n_sample,r_sample)]
pr_f=f[,c(p_sample,r_sample)]

np_p =f_pdata[c(n_sample,p_sample),]
nr_p =f_pdata[c(n_sample,r_sample),]
pr_p =f_pdata[c(p_sample,r_sample),]

table(np_p$n_vs_p_vs_r_group)#n 52 p 369
table(nr_p$n_vs_p_vs_r_group)#n 52 r 58
table(pr_p$n_vs_p_vs_r_group)#p 369 r 58

dds_p <- DESeqDataSetFromMatrix(
  countData = np_f,
  colData=np_p,
  design=~n_vs_p_vs_r_group, tidy=F
)

dds_r <- DESeqDataSetFromMatrix(
  countData = nr_f,
  colData=nr_p,
  design=~n_vs_p_vs_r_group, tidy=F
)

dds_pr <- DESeqDataSetFromMatrix(
  countData = pr_f,
  colData=pr_p,
  design=~n_vs_p_vs_r_group, tidy=F
)

dds_p <- DESeq(dds_p)
dds_r <- DESeq(dds_r)
dds_pr <- DESeq(dds_pr)
#요기
res_p <- results(dds_p)
res_r <- results(dds_r)
res_pr <- results(dds_pr)

#
res_p <- as.data.frame(res_p)
res_r <- as.data.frame(res_r)
res_pr <- as.data.frame(res_pr)

write.csv(n_vs_t_res,file='./data/prad_zerocnt_0.1_res.csv')
write.csv(res_p,file='./data/prad_zerocnt_0.1_res_p.csv')
write.csv(res_r,file='./data/prad_zerocnt_0.1_res_r.csv')
write.csv(res_pr,file='./data/prad_zerocnt_0.1_res_pr.csv')


library(dplyr)
res_p <- res_p[res_p$pvalue<=0.05& abs(res_p$log2FoldChange)>=1.5,]
dim(res_p)#1565 6
res_p <- as.data.frame(res_p[order(res_p$log2FoldChange),])
head(res_p)
res_p[1:25,]

p_bottom=rownames(res_p)[1:25]
res_p <- res_p[order(-res_p$log2FoldChange),]
p_top=rownames(res_p)[1:25]
#
res_r <- res_r[res_r$pvalue<=0.05& abs(res_r$log2FoldChange)>=1.5,]
res_r<- as.data.frame(res_r[order(res_r$log2FoldChange),])
dim(res_r)#2273 6
res_r[1:25,]

r_bottom=rownames(res_r)[1:25]
res_r <- res_r[order(-res_r$log2FoldChange),]
r_top=rownames(res_r)[1:25]
#
res_pr <- res_pr[res_pr$pvalue<=0.05& abs(res_pr$log2FoldChange)>=1.5,]
res_pr<- as.data.frame(res_pr[order(res_pr$log2FoldChange),])
dim(res_pr)#78 6
head(res_pr)
res_pr[1:25,]

o_pr_bottom=rownames(res_pr)[1:25]
res_pr <- res_pr[order(-res_pr$log2FoldChange),]
o_pr_top=rownames(res_pr)[1:7]

###################
#
hist(counts(dds_p, normalized=T))
log2_count <- log2(counts(dds_p, normalized=T)+1)
range(log2_count)#0 23.97126
hist(log2_count)
scale_log2_count <- t(scale(t(log2_count)))
range(scale_log2_count)#-12.62 12.19


hist(scale_log2_count)

library(pheatmap)
df <-as.data.frame(colData(dds_p)[,c('n_vs_p_vs_r_group')])
rownames(df)=colnames(scale_log2_count)
colnames(df)='Type'

pheatmap(scale_log2_count[c(p_top,p_bottom),], cluster_rows = F,show_rownames =T, show_colnames = F, cluster_cols = F, annotation_col = df,breaks=seq(-1,1, length.out = 100))
res_p[p_top,]
###################
hist(counts(dds_r, normalized=T))
log2_count <- log2(counts(dds_r, normalized=T)+1)
range(log2_count)#0 24.05558
hist(log2_count)
scale_log2_count <- t(scale(t(log2_count)))
range(scale_log2_count)#-9.3 7.9


hist(scale_log2_count)

library(pheatmap)
df <-as.data.frame(colData(dds_r)[,c('n_vs_p_vs_r_group')])
rownames(df)=colnames(scale_log2_count)
colnames(df)='Type'

pheatmap(scale_log2_count[c(p_top,p_bottom),], cluster_rows = F,show_rownames =T, show_colnames = F, cluster_cols = F, annotation_col = df,breaks=seq(-1,1, length.out = 100))
res_p[p_top,]
###################
hist(counts(dds_pr, normalized=T))
log2_count <- log2(counts(dds_pr, normalized=T)+1)
range(log2_count)#0 22.0171
hist(log2_count)
scale_log2_count <- t(scale(t(log2_count)))
range(scale_log2_count)#-12 12


hist(scale_log2_count)

library(pheatmap)
df <-as.data.frame(colData(dds_r)[,c('n_vs_p_vs_r_group')])
rownames(df)=colnames(scale_log2_count)
colnames(df)='Type'

pheatmap(scale_log2_count[c(o_pr_top,o_pr_bottom),], cluster_rows = F,show_rownames =T, show_colnames = F, cluster_cols = F, annotation_col = df,breaks=seq(-1,1, length.out = 100))
res_pr[o_pr_top,]

heatmap(scale_log2_count[c(o_pr_top,o_pr_bottom),])

