fpkm_only_T = read.table('./data/crpc_T_fpkm_31.txt',sep='\t')
head(fpkm_only_T)
dim(fpkm_only_T)
cnt.zero = rowSums(fpkm_only_T == 0)
table(cnt.zero)
hist(cnt.zero)

colnames(fpkm_only_T)=lapply(colnames(fpkm_only_T), function(x) substr(x,2,6))

fbxo_septin14_inc = fpkm_only_T[,as.numeric(substr(colnames(fpkm_only_T),1,2)) %in% c(08,14,25,29,34,15,03,12,39,11,06)]
fbxo_septin14_noinc = fpkm_only_T[,!as.numeric(substr(colnames(fpkm_only_T),1,2)) %in% c(08,14,25,29,34,15,03,12,39,11,06)]

fpkm_only_T = cbind(fbxo_septin14_inc,fbxo_septin14_noinc)
dim(fpkm_only_T)
dim(fbxo_septin14_inc)

#zero
zero_0.3_fpkm_only_T=fpkm_only_T[-which(cnt.zero > round(ncol(fpkm_only_T)*0.3)),]
#log2
range(zero_0.3_fpkm_only_T) # 0 - 587530.6
log2_zero_0.3_fpkm_only_T = log2(zero_0.3_fpkm_only_T+1)
range(log2_zero_0.3_fpkm_only_T)# 0 -19.16431
boxplot(log2_zero_0.3_fpkm_only_T)
#quantile
library(preprocessCore)
norm_log2_zero_0.3_fpkm_only_T=normalize.quantiles(as.matrix(log2_zero_0.3_fpkm_only_T))
norm_log2_zero_0.3_fpkm_only_T=t(scale(t(norm_log2_zero_0.3_fpkm_only_T)))
boxplot(norm_log2_zero_0.3_fpkm_only_T)
hist(apply(norm_log2_zero_0.3_fpkm_only_T,2,as.numeric))

#immunecellai 넣을려고 만드는 파일
table(is.na(norm_log2_zero_0.3_fpkm_only_T))# all false
rownames(norm_log2_zero_0.3_fpkm_only_T) = rownames(log2_zero_0.3_fpkm_only_T)
data("fbxo_septin_feature")
norm_log2_zero_0.3_fpkm_only_T = as.data.frame(norm_log2_zero_0.3_fpkm_only_T)
range(norm_log2_zero_0.3_fpkm_only_T)
colnames(norm_log2_zero_0.3_fpkm_only_T) = rownames(feature)
Fbox25_septin14_0.3_immunecellai = rbind(t(feature), norm_log2_zero_0.3_fpkm_only_T)
head(Fbox25_septin14_0.3_immunecellai)

write.table(Fbox25_septin14_0.3_immunecellai,file='/data2/fusion/fbxo25_septin14_0.3_fpkm.txt',
          sep='\t')

#실제로 돌릴 때 첫번째 열
#symbol sam1 sam2 sam3 이렇게 되도록 수정


