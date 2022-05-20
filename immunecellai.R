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

write.table(zero_0.3_fpkm_only_T,file='./test_remove.txt',sep='\t')

#log2 이거 하지 말기 ! => 왜냐면 immune cell ai에 넣으면 log2를 해주거든! 
range(zero_0.3_fpkm_only_T) # 0 - 587530.6
log2_zero_0.3_fpkm_only_T = zero_0.3_fpkm_only_T
#log2_zero_0.3_fpkm_only_T = log2(zero_0.3_fpkm_only_T+1)
#range(log2_zero_0.3_fpkm_only_T)# 0 -19.16431
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

write.table(Fbox25_septin14_0.3_immunecellai,file='/data2/fusion/fbxo25_septin14_0.3_fpkm_nolog2.txt',
          sep='\t')

#실제로 돌릴 때 첫번째 열
#symbol sam1 sam2 sam3 이렇게 되도록 수정

##################
#결과분석

fbxo25_septin14_res=read.table(file='./data2/fusion/ImmuCellAI_abundance_result_0.3_fpkm.txt')
fbxo25_septin14_res$fusion = c(rep('Fbxo25_septin14_Inc',11),
                               rep('Fbxo25_septin14_noInc',20))

head(fbxo25_septin14_res)

for_plot=as.data.frame(fbxo25_septin14_res[,c('DC','fusion')],row.names = 1:nrow(fbxo25_septin14_res))
for_plot$immune=c(rep('DC',nrow(for_plot)))
colnames(for_plot)[1]='value'

num = ncol(fbxo25_septin14_res)-1

for (i in 2:num){
  col = colnames(fbxo25_septin14_res)[i]
  ins=cbind(fbxo25_septin14_res[,c(col,'fusion')],c(rep(col,
                                           length(fbxo25_septin14_res[[i]])))) %>%
    as.data.frame(row.names = 1:nrow(.))
  colnames(ins)[3]='immune'
  colnames(ins)[1]='value'
  for_plot=rbind(for_plot,ins)
  
}
head(for_plot)
for_plot =for_plot[-which(for_plot$immune == 'InfiltrationScore'),]
for_plot %>%ggplot(aes(x=immune,y=value,fill=fusion))+ geom_boxplot()+
  theme(axis.text.x=element_text(angle = 90, hjust=1))








