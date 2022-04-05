data("GTEX_54_selected") #sample로 저장
data("f_merge_gene") #f_merge_gene으로 저장

dim(f_merge_gene)#36506 984 biomart로 gene symbol로 바꿔 놓기만 한 상태
head(f_merge_gene)

T_sample = colnames(f_merge_gene)[1:36]

#extract 56 people from gtex data and merge t sample
54+36 #90

f_merge_gene_90 = f_merge_gene[, c(T_sample,sample)]

dim(f_merge_gene_90)#36506 90

cnt.zero = rowSums(f_merge_gene_90==0)
table(cnt.zero)
hist(cnt.zero)

36506*0.9 #32855
#전체 10%,30%,50%,100%

#preprocess_deg
library(preprocessCore)
preprocess_deg <- function(zero_threshold, count){

  cnt.zero = rowSums(count == 0)
  
  threshold = round(ncol(count)*zero_threshold)
  message('threshold is ')
  print(threshold)
  
  zero_count = count[which(cnt.zero < threshold),]
  message('zero_count_dim is')
  print(dim(zero_count))
  
  type = c(rep('T',36),rep('N',54))
  sample =colnames(zero_count)
  
  feature = as.data.frame(cbind(sample,type))
  rownames(feature) = sample
  
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
  R = preprocess_deg(i, f_merge_gene_90)
  colnames(f_merge_gene_90)
  zero = as.data.frame(R[1])
  
  col = as.data.frame(R[3])
  
  colnames(zero) = col$c..01_expected_count....04_expected_count....08_expected_count...
  
  feature= as.data.frame(R[2])
  
  write.table(zero,file=paste0('./data/F0405_preprocess_zero',i,'.txt'), sep='\t')
  
  D=DEG(zero, feature)
  colnames(zero) == rownames(feature)
  res_na=as.data.frame(D[1])
  vst=as.data.frame(D[2]) #test
  colnames(vst) = col$c..01_expected_count....04_expected_count....08_expected_count...
  dim(vst)
  write.table(vst, file=paste0('./data/F0405_for_heatmap_',i,'.txt'), sep='\t')
  
  write.table(res,file=paste0('./data/F0405_res_',i,'.txt'), sep='\t')
  
  final_res=res_na[res_na$pvalue<0.05 & abs(res_na$log2FoldChange)>=1.5,]
  write.table(res, file= paste0('./data/F0405_res_',i,'_p0.05_log2_1.5.txt'), sep='\t')
  
  up_final_res = final_res[final_res$log2FoldChange >=1.5,]
  down_final_res = final_res[final_res$log2FoldChange <=-1.5,]
  
  print('up_final_res')
  print(dim(up_final_res))
  
  print('down_final_res')
  print(dim(down_final_res))
  
  
  write.table(up_final_res, file= paste0('./data/F0405_up_res_',i,'_p0.05_log2_1.5.txt'), sep='\t')
  write.table(down_final_res, file= paste0('./data/F0405_down_res_',i,'_p0.05_log2_1.5.txt'), sep='\t')
  
  if (nrow(up_final_res) >=25 && nrow(down_final_res) >=25){
    up_final_res=up_final_res[order(-up_final_res$log2FoldChange),]
    print(head(up_final_res))
    top_25_gene = rownames(up_final_res[1:25,])
    write.table(rownames(up_final_res[1:30,]), file= paste0('./data/F0405_up_30_',i,'.txt'), sep='\t')
    
    down_final_res=down_final_res[order(down_final_res$log2FoldChange),]
    print(head(down_final_res))
    bottom_25_gene = rownames(down_final_res[1:25,])
    
    write.table(rownames(down_final_res[1:30,]), file= paste0('./data/F0405_down_30_',i,'.txt'), sep='\t')
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

plot_0.1 =get('final_50_genes_0.1')
plot_0.3 =get('final_50_genes_0.3')
plot_0.5 =get('final_50_genes_0.5')
plot_1 =get('final_50_genes_1')

vst = apply(vst,1,as.numeric)
hist(vst)
#boxplot(vst)

anno = as.data.frame(feature[,'type'])
rownames(anno) = rownames(feature)
colnames(anno)[1] = 'type'

range(plot_0.1)
pheatmap(plot_0.1, cluster_rows = F, cluster_cols = F, annotation_col=anno,
         breaks=seq(-3,3,length.out=100), show_colnames = F)
