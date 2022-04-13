#n,t,r구분
setwd('D:/')
data("tcga_prostate_gene_clinical_0222_eset")
pdata= pData(tcga_prostate_gene_clinical_eset)
dim(pdata)#548 105

exp=exprs(tcga_prostate_gene_clinical_eset)
fdata=fData(tcga_prostate_gene_clinical_eset)

r_pdata <- pdata[pdata$biochemical_recurrence=='YES'& substr(pdata$sampleID,14,15)=='01' |substr(pdata$sampleID,14,15)=='06' ,]

dim(r_pdata)#59 105
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
dim(f_pdata)#480 105
dim(f_exp)#38662 480

f_pdata$t_vs_n_group=c(rep('n',52),rep('t',428))
f_pdata$n_vs_p_vs_r_group = c(rep('n',52),rep('p',369),rep('r',59))

feature_data = as.data.frame(rownames(f_exp))
colnames(feature_data)='gene'
rownames(feature_data)= feature_data$gene

featureData=new("AnnotatedDataFrame", data=feature_data) 
phenoData = new("AnnotatedDataFrame", data=f_pdata)

rownames(f_exp) == rownames(feature_data)
colnames(f_exp)== rownames(f_pdata)
tcga_crpc_grouped_eset= new('ExpressionSet', exprs=f_exp, phenoData=phenoData, featureData=featureData)

save(tcga_crpc_grouped_eset,file = './data/tcga_crpc_grouped_eset.rda')


###############
data("tcga_crpc_grouped_eset")
#n vs t
59+369 #428
n_t=f_exp[,c(n_sample,p_sample,r_sample)]
dim(n_t)#38662 111, 480

nt_p =f_pdata[c(n_sample,p_sample,r_sample),]
table(rownames(nt_p)==colnames(n_t))#all true 110
nt_p$
library(preprocessCore)
preprocess_deg <- function(zero_threshold, count){
  
  cnt.zero = rowSums(count == 0)
  
  threshold = round(ncol(count)*zero_threshold)
  message('threshold is ')
  print(threshold)
  
  zero_count = count[which(cnt.zero < threshold),]
  message('zero_count_dim is')
  print(dim(zero_count))
  
  type = c(rep('n',52),rep('t',428))
  sample =colnames(zero_count)
  
  feature = as.data.frame(cbind(sample,type))
  rownames(feature) = sample
  
  col=colnames(zero_count)
  
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
  R = preprocess_deg(i, n_t)
  colnames(n_t)
  zero = as.data.frame(R[1])
  
  col = as.data.frame(R[3])
  
  colnames(zero) = colnames(n_t)
  
  feature= as.data.frame(R[2])
  
  write.table(zero,file=paste0('./data/tcga_nt_0413_preprocess_zero',i,'.txt'), sep='\t')
  
  D=DEG(zero, feature)
  colnames(zero) == rownames(feature)
  res_na=as.data.frame(D[1])
  vst=as.data.frame(D[2]) #test
  colnames(vst) = colnames(n_t)
  dim(vst)
  write.table(vst, file=paste0('./data/tcga_nt_0413_for_heatmap_',i,'.txt'), sep='\t')
  
  write.table(res_na,file=paste0('./data/tcga_nt_0413_res_',i,'.txt'), sep='\t')
  
  final_res=res_na[res_na$pvalue<0.05 & abs(res_na$log2FoldChange)>=1.5,]
  write.table(final_res, file= paste0('./data/tcga_nt_0413_res_',i,'_p0.05_log2_1.5.txt'), sep='\t')
  
  up_final_res = final_res[final_res$log2FoldChange >=1.5,]
  down_final_res = final_res[final_res$log2FoldChange <=-1.5,]
  
  print('up_final_res')
  print(dim(up_final_res))
  
  print('down_final_res')
  print(dim(down_final_res))
  
  
  write.table(up_final_res, file= paste0('./data/tcga_nt_0413_up_res_',i,'_p0.05_log2_1.5.txt'), sep='\t')
  write.table(down_final_res, file= paste0('./data/tcga_nt_0413_down_res_',i,'_p0.05_log2_1.5.txt'), sep='\t')
  
  if (nrow(up_final_res) >=25 && nrow(down_final_res) >=25){
    up_final_res=up_final_res[order(-up_final_res$log2FoldChange),]
    print(head(up_final_res))
    top_25_gene = rownames(up_final_res[1:25,])
    write.table(rownames(up_final_res[1:30,]), file= paste0('./data/tcga_nt_0413_up_30_',i,'.txt'), sep='\t')
    
    down_final_res=down_final_res[order(down_final_res$log2FoldChange),]
    print(head(down_final_res))
    bottom_25_gene = rownames(down_final_res[1:25,])
    
    write.table(rownames(down_final_res[1:30,]), file= paste0('./data/tcga_nt_0413_down_30_',i,'.txt'), sep='\t')
  }
  else{
    print('니가 해')
  }
  
  final_50_genes= c(top_25_gene, bottom_25_gene)
  print(final_50_genes)
  print(head(vst))
  
  as.data.frame(vst)
  a=vst[final_50_genes,]
  colnames(a) =colnames(n_t)
  assign(paste0('final_50_genes_',i),a)
  
}

plot_0.1 =get('final_50_genes_0.1')
plot_0.3 =get('final_50_genes_0.3')
plot_0.5 =get('final_50_genes_0.5')
plot_1 =get('final_50_genes_1')

vst = apply(vst,1,as.numeric)
hist(vst)
range(vst)
mean(vst)
boxplot(vst)
memory.limit(1000000)

anno = as.data.frame(feature[,'type'])
rownames(anno) = rownames(feature)
colnames(anno)[1] = 'type'

range(plot_0.1)
library(pheatmap)
range(plot_0.1)
pheatmap(plot_0.1, cluster_rows = F, cluster_cols = F, annotation_col=anno,
         breaks=seq(-3,3,length.out=100), show_colnames = F)

###############
##########
#n vs r
52+59
nr_f=f_exp[,c(n_sample,r_sample)]
dim(nr_f)#38662 111

nr_p =f_pdata[c(n_sample,r_sample),]
table(rownames(nr_p)==colnames(nr_f))#all true 110
nr_p$n_vs_p_vs_r_group
preprocess_deg <- function(zero_threshold, count){
  
  cnt.zero = rowSums(count == 0)
  
  threshold = round(ncol(count)*zero_threshold)
  message('threshold is ')
  print(threshold)
  
  zero_count = count[which(cnt.zero < threshold),]
  message('zero_count_dim is')
  print(dim(zero_count))
  
  type = c(rep('n',52),rep('r',59))
  sample =colnames(zero_count)
  
  feature = as.data.frame(cbind(sample,type))
  rownames(feature) = sample
  
  col=colnames(zero_count)
  
  return(list(zero_count, feature, col))
}

for (i in c(0.1,0.3,0.5,1)){ 
  R = preprocess_deg(i, nr_f)
  colnames(nr_f)
  zero = as.data.frame(R[1])
  
  col = as.data.frame(R[3])
  
  colnames(zero) = colnames(nr_f)
  
  feature= as.data.frame(R[2])
  
  write.table(zero,file=paste0('./data/tcga_nr_0413_preprocess_zero',i,'.txt'), sep='\t')
  
  D=DEG(zero, feature)
  colnames(zero) == rownames(feature)
  res_na=as.data.frame(D[1])
  vst=as.data.frame(D[2]) #test
  colnames(vst) = colnames(nr_f)
  dim(vst)
  write.table(vst, file=paste0('./data/tcga_nr_0413_for_heatmap_',i,'.txt'), sep='\t')
  
  write.table(res_na,file=paste0('./data/tcga_nr_0413_res_',i,'.txt'), sep='\t')
  
  final_res=res_na[res_na$pvalue<0.05 & abs(res_na$log2FoldChange)>=1.5,]
  write.table(final_res, file= paste0('./data/tcga_nr_0413_res_',i,'_p0.05_log2_1.5.txt'), sep='\t')
  
  up_final_res = final_res[final_res$log2FoldChange >=1.5,]
  down_final_res = final_res[final_res$log2FoldChange <=-1.5,]
  
  print('up_final_res')
  print(dim(up_final_res))
  
  print('down_final_res')
  print(dim(down_final_res))
  
  
  write.table(up_final_res, file= paste0('./data/tcga_nr_0413_up_res_',i,'_p0.05_log2_1.5.txt'), sep='\t')
  write.table(down_final_res, file= paste0('./data/tcga_nr_0413_down_res_',i,'_p0.05_log2_1.5.txt'), sep='\t')
  
  if (nrow(up_final_res) >=25 && nrow(down_final_res) >=25){
    up_final_res=up_final_res[order(-up_final_res$log2FoldChange),]
    print(head(up_final_res))
    top_25_gene = rownames(up_final_res[1:25,])
    write.table(rownames(up_final_res[1:30,]), file= paste0('./data/tcga_nr_0413_up_30_',i,'.txt'), sep='\t')
    
    down_final_res=down_final_res[order(down_final_res$log2FoldChange),]
    print(head(down_final_res))
    bottom_25_gene = rownames(down_final_res[1:25,])
    
    write.table(rownames(down_final_res[1:30,]), file= paste0('./data/tcga_nr_0413_down_30_',i,'.txt'), sep='\t')
  }
  else{
    print('니가 해')
  }
  
  final_50_genes= c(top_25_gene, bottom_25_gene)
  print(final_50_genes)
  print(head(vst))
  
  as.data.frame(vst)
  a=vst[final_50_genes,]
  colnames(a) =colnames(nr_f)
  assign(paste0('final_50_genes_',i),a)
  
}

plot_0.1 =get('final_50_genes_0.1')
plot_0.3 =get('final_50_genes_0.3')
plot_0.5 =get('final_50_genes_0.5')
plot_1 =get('final_50_genes_1')
