#immunecellai 넣을려고 만드는 파일
vst = read.table('./data/Fbox25_septin14_for_heatmap_0.3.txt', sep='\t')
vst = as.data.frame(vst)

data("31_feature")

colnames(vst) = rownames(feature)
Fbox25_septin14_0.3_immunecellai = rbind(t(feature), vst)
rownames(Fbox25_septin14_0.3_immunecellai)
head(Fbox25_septin14_0.3_immunecellai)
write.table(Fbox25_septin14_0.3_immunecellai,file='/data2/fusion/fbox25_septin14_immunecellai_0.3.txt',
          sep='\t')

#실제로 돌릴 때 첫번째 열
#symbol sam1 sam2 sam3 이렇게 되도록 수정



library(ggplot2)
fbox25_septin24_immunecellai_result = read.table('./data2/fusion/ImmuCellAI_group_result.txt',sep='\t', header=1)
head(fbox25_septin24_immunecellai_result)

ggplot(fbox25_septin24_immunecellai_result,
       aes(fill=))


test = read.table('C:/Users/User/Downloads/ImmuCellAI_example.txt',sep='\t', header=1)
head(test)
rownames(test) = test$Symbol
head(test)
#test = test[-1,-1]
table(is.na(test))

apply(test,2,range)#log2 한거맞음? 겁나 커


data("f_merge_gene_31")
#f_merge_gene_31
data("31_feature")

raw_Fbox25_septin14_0.3_immunecellai = rbind(t(feature), f_merge_gene_31)
rownames(raw_Fbox25_septin14_0.3_immunecellai)
head(raw_Fbox25_septin14_0.3_immunecellai)
write.table(raw_Fbox25_septin14_0.3_immunecellai,file='/data2/fusion/fbox25_septin14_immunecellai_raw.txt',
            sep='\t')
