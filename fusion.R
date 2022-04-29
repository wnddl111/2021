setwd('D:/data2/fusion')

#data 
folders = read.table('./foldername.txt')
for ( i in folders$V1 ){
  assign(i,readr::read_tsv(paste0(i,'/star-fusion.fusion_predictions.abridged.coding_effect.tsv')))
  assign(i,cbind(sample_name = rep(i, nrow(get(i))),get(i)))
  assign(paste0(i,'_inf'), get(i)[which(get(i)[['PROT_FUSION_TYPE']]=='INFRAME'),])
}

#inframe 인것만해서 파일 하나로 합치기
all_inf = TURP01_inf

for (i in folders$V1[-1]){
  assign('all_inf', rbind(all_inf, get(paste0(i,'_inf'))))
}

dim(all_inf)#190 28
library(dplyr)


all_inf = as.data.frame(all_inf)
dim(all_inf)
rownames(all_inf) = NULL
library(writexl)
write_xlsx(all_inf, './all_inf.xlsx')

#sample별로 몇개의 fusion이 있는지 
#inf인것만 뽑으면 4명의 환자가 빠지게 됨
num_of_fusion_per_sample <- all_inf[which(all_inf$SpliceType=='ONLY_REF_SPLICE'),] %>%
  group_by(sample_name) %>%
  summarise(count=n()) %>%
  as.data.frame()

library(ggplot2)

num_of_fusion_per_sample$type = ifelse(num_of_fusion_per_sample$sample_name %in% c('TURP01',
                                                   'TURP04',
                                                   'TURP05',
                                                   'TURP13',
                                                   'TURP27'), 'HSPC','CRPC')

ggplot(num_of_fusion_per_sample,
       aes(x=sample_name, y= count, fill=type))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle=90))+
  geom_text(aes(label=count, vjust = -1))

################
#install.packages('BioCircos')

library(BioCircos)
hg38 = list('1'=248956422,
'2'=242193529,
'3'=198295559,
'4'=190214555,
'5'=181538259,
'6'=170805979,
'7'=159345973,
'8'=145138636,
'9'=138394717,
'10'=133797422,
'11'=135086622,
'12'=133275309,
'13'=114364328,
'14'=107043718,
'15'=101991189,
'16'=90338345,
'17'=83257441,
'18'=80373285,
'19'=58617616,
'20'=64444167,
'21'=46709983,
'22'=50818468,
'X'=156040895,
'Y'=57227415
)

BioCircos(genome = hg38)

rank1_read1 = c('8')
rank2_read1 = '6'
rank3_read1=c('21',
              '14',
              '14'
)      
rank4_read1 = c('21',
                '14',
                '14',
                '21',
                '22',
                '5',
                '4',
                '20',
                '14',
                '15'
)
rank1_read2 = c('7')
rank2_read2 = c('6')
rank3_read2 = c('21','14','14')
rank4_read2 = c(               '21',
                               '14',
                               '14',
                               '21',
                               '22',
                               '5',
                               '4',
                               '20',
                               '14',
                               '10'
)

rank1_read_pos1=c(
  435707)
  
rank2_read_pos1=c(147509387)
rank3_read_pos1=c(39793613,
37623003,
102195120)
rank4_read_pos1 = c(41498119,
                    42601970,
                    33982236,
                    113609152,
                    2975241,
                    69658713,
                    70057474
)

rank1_read_pos2 = c(
  55796092)
rank2_read_pos2=148390134
rank3_read_pos2=c(
  39898476,
37422855,
101762905)

rank4_read_pos2=c(38445621,
                  45372292,
                  34005899,
                  113906745,
                  2751700,
                  68594485,
                  132910068
                  
)

rank1_labels='FBXO25--SEPTIN14'
rank2_labels='SAMD5--SASH1'
rank3_labels=c('IGSF5--PCP4',
'TTC6--MIPOL1',
'WDR20--PPP2R5C')
rank4_labels=c(
'TMPRSS2--ERG',
'POLDIP3--SMC1B',
'SLC45A2--AMACR',
'CAMK2D--ARSJ',
'PTPRA--EBF4',
'SUSD6--RAD51B',
'TLE3--CFAP46')

 
tracklist1 = BioCircosLinkTrack('myLinkTrack', rank1_read1, rank1_read_pos1,
                                           rank1_read_pos1+50000000 , rank1_read2, rank1_read_pos2, rank1_read_pos2 + 750000,
                                           maxRadius = 0.55, labels = rank1_labels, width =5.5)
tracklist2 = tracklist1+BioCircosLinkTrack('myLinkTrack', rank2_read1, rank2_read_pos1,
                                rank2_read_pos1+50000000 , rank2_read2, rank2_read_pos2, rank2_read_pos2 + 750000,
                                maxRadius = 0.55, labels = rank2_labels, width =3.5)
tracklist3 = tracklist2+BioCircosLinkTrack('myLinkTrack', rank3_read1, rank3_read_pos1,
                                rank3_read_pos1+50000000 , rank3_read2, rank3_read_pos2, rank3_read_pos2 + 750000,
                                maxRadius = 0.55, labels = rank3_labels, width =1.5)
tracklist4 = tracklist3+BioCircosLinkTrack('myLinkTrack', rank4_read1, rank4_read_pos1,
                                rank4_read_pos1+50000000 , rank4_read2, rank4_read_pos2, rank4_read_pos2 + 750000,
                                maxRadius = 0.55, labels = rank4_labels, width =1)
#, width = '0.3em')
BioCircos(tracklist4 ,genome = hg38)

dev.off()

?BioCircosLineTrack


########################
#Top1인 fbx025-septin14로 그룹을 나눠서 deg 확인해볼까
setwd('D:/')
data(f_merge_gene)
dim(f_merge_gene)#36506 984 biomart로 gene symbol로 바꿔 놓기만 한 상태
head(f_merge_gene)

T_sample = colnames(f_merge_gene)[1:36]

f_merge_gene_36 = f_merge_gene[, c(T_sample)]

dim(f_merge_gene_36)#36506 36

cnt.zero = rowSums(f_merge_gene_36==0)
table(cnt.zero)
hist(cnt.zero)

#ㄹㅇ crpc만 뽑자 (hspc 1,4,5,13,27), 1,2,11,14,23 index
colnames(f_merge_gene_36)
f_merge_gene_31 = f_merge_gene_36[,-c(1,2,11,14,23)]
save(f_merge_gene_31,file='./data/f_merge_gene_31.rda')

Fusion = c('08_expected_count',
           '14_expected_count',
           '25_expected_count',
           '29_expected_count',
           '34_expected_count',
           '15_expected_count',
           '03_expected_count',
           '12_expected_count',
           '39_expected_count',
           '11_expected_count',
           '06_expected_count'
)

type = as.data.frame(ifelse(colnames(f_merge_gene_31) %in% Fusion, 'FBOX25_SEPTIN14_INC', 'FBOX25_SEPTIN14_noINC'))
colnames(type)='type'
rownames(type) = colnames(f_merge_gene_31)


preprocess_deg <- function(zero_threshold, count){
  
  cnt.zero = rowSums(count == 0)
  
  threshold = round(ncol(count)*zero_threshold)
  message('threshold is ')
  print(threshold)
  
  zero_count = count[which(cnt.zero < threshold),]
  message('zero_count_dim is')
  print(dim(zero_count))
  

  feature = type
  col = rownames(feature)
  return(list(zero_count, feature, col))
}
library(DESeq2)
DEG <- function(zero_count, coldata){
  dds <-DESeqDataSetFromMatrix(countData = round(zero_count),
                               colData = coldata,
                               design= ~type,
                               tidy=F)
  dds <- DESeq(dds)
  print(table(is.na(dds)))
  dds_na = na.omit(dds)
  
  vst <- vst(dds,blind = F)
  vst <- assay(vst)
  vst <- as.data.frame(vst)
  
  vst <- (t(scale(t(as.matrix(vst)))))
  
  
  res <- results(dds_na)
  print(table(is.na(res)))
  
  res_na <- na.omit(res)
  return(list(res_na, vst))
}

for (i in c(0.1,0.3,0.5,1)){ 
  R = preprocess_deg(i, f_merge_gene_31)
  zero = as.data.frame(R[1])
  
  col = as.data.frame(R[3])
  
  colnames(zero) = col[[1]]
  
  feature= as.data.frame(R[2])
  
  write.table(zero,file=paste0('./data/Fbox25_septin14_preprocess_zero',i,'.txt'), sep='\t')
  
  D=DEG(zero, feature)
  colnames(zero) == rownames(feature)
  res_na=as.data.frame(D[1])
  vst=as.data.frame(D[2]) #test
  colnames(vst) = col[[1]]
  dim(vst)
  write.table(vst, file=paste0('./data/Fbox25_septin14_for_heatmap_',i,'.txt'), sep='\t')
  
  write.table(res_na,file=paste0('./data/Fbox25_septin14_res_',i,'.txt'), sep='\t')
  
  final_res=res_na[res_na$pvalue<0.05 & abs(res_na$log2FoldChange)>=1.5,]
  write.table(final_res, file= paste0('./data/Fbox25_septin14_res_',i,'_p0.05_log2_1.5.txt'), sep='\t')
  
  up_final_res = final_res[final_res$log2FoldChange >=1.5,]
  down_final_res = final_res[final_res$log2FoldChange <=-1.5,]
  
  print('up_final_res')
  print(dim(up_final_res))
  
  print('down_final_res')
  print(dim(down_final_res))
  
  
  write.table(up_final_res, file= paste0('./data/Fbox25_septin14_up_res_',i,'_p0.05_log2_1.5.txt'), sep='\t')
  write.table(down_final_res, file= paste0('./data/Fbox25_septin14_down_res_',i,'_p0.05_log2_1.5.txt'), sep='\t')
  
  if (nrow(up_final_res) >=25 && nrow(down_final_res) >=25){
    up_final_res=up_final_res[order(-up_final_res$log2FoldChange),]
    print(head(up_final_res))
    top_25_gene = rownames(up_final_res[1:25,])
    write.table(rownames(up_final_res[1:30,]), file= paste0('./data/Fbox25_septin14_up_30_',i,'.txt'), sep='\t')
    
    down_final_res=down_final_res[order(down_final_res$log2FoldChange),]
    print(head(down_final_res))
    bottom_25_gene = rownames(down_final_res[1:25,])
    
    write.table(rownames(down_final_res[1:30,]), file= paste0('./data/Fbox25_septin14_down_30_',i,'.txt'), sep='\t')
  }
  else{
    print('니가 해')
  }
  
  final_50_genes= c(top_25_gene, bottom_25_gene)
  print(final_50_genes)
  print(head(vst))
  
  as.data.frame(vst)
  a=vst[final_50_genes,]
  colnames(a) =col$c..01_expected_count....04_expected_count....08_expected_count...
  assign(paste0('final_50_genes_',i),a)
  
}

#
plot_1 =get('final_50_genes_1')

range(plot_1)
colnames(plot_1) = rownames(type)
library(pheatmap)

hist(apply(plot_1,2,as.numeric))

not_fusion = colnames(plot_1[,-which(colnames(plot_1) %in% Fusion)])
s_plot_1 = plot_1[,c(Fusion,not_fusion)]
pheatmap(s_plot_1, cluster_rows = F, cluster_cols = F, annotation_col=feature,
         breaks=seq(-1,1,length.out=100), show_colnames = F)










