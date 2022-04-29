kegg_ref=read.table('./data/kegg_ref.txt',sep='\n')
head(kegg_ref)
dim(kegg_ref)#125 1
kegg_ref = as.data.frame(kegg_ref)
kegg_ref[[1]][3]

kegg = list()
for (i in 1:125){
  j = substr(kegg_ref[[1]][i],1,2)
  if (startsWith(kegg_ref[[1]][i],paste0(j,' '))){
    print(kegg_ref[[1]][i])
  }
  else {
    kegg = append(kegg, kegg_ref[[1]][i])
  }
}


kegg[[2]] 
main_pathway = list()
sub_pathway = list()
for (i in 1:118){
  if ( i%%2 == 0){
    print(kegg[[i]])
    sub_pathway = append(sub_pathway, kegg[[i]])
  }
  else{
    main_pathway = append(main_pathway, kegg[[i]])
  
  }
}



write.table(main_pathway,'./data/kegg_main.txt',sep = '\n')
write.table(sub_pathway,'./data/kegg_sub.txt',sep='\n')


#########
#jupyter랑 엑셀에서 파일 정리하고 여기서는 gene list 불러올 것
getwd()
##
repair_genes = as.data.frame(readxl::read_excel('./data/repair_gene_list.xlsx'))
repair_genes = repair_genes$repair

metabolism_genes = as.data.frame(readxl::read_excel('./data/metabolism_gene_list.xlsx'))
metabolism_genes =metabolism_genes$metabolism

immune_genes = as.data.frame(readxl::read_excel('./data/immune_gene_list.xlsx'))
immune_genes =immune_genes$immune

signal_genes = as.data.frame(readxl::read_excel('./data/signal_gene_list.xlsx'))
signal_genes =signal_genes$signal

#####

exp = read.table(file='./data/F0405_for_heatmap_0.3.txt', sep='\t')
dim(exp)#18008 90
head(exp) # ZRANB2-AS2 이 유전자가 다 na값인데 omit 안된채로 저장돼있음 

gene = exp$V1
exp = exp[,-1]
colnames(exp) = rownames(anno)
rownames(exp) = gene

res = read.table(file='./data/F0405_res_0.3.txt')
dim(res)#18006 6
########
res[repair_genes,]

library(pheatmap)
pheatmap(exp[repair_genes,], cluster_rows = F, cluster_cols = F, annotation_col=anno,
         breaks=seq(-3,3,length.out=100), show_colnames = F)
#metabolism gene은 365개로 너무 많아서 fold change가 몇 이상인것만 뽑아서 그리기로 함
a=res[metabolism_genes,]
intersect(metabolism_genes,rownames(res))

a=a[a$log2FoldChange>=5,]
dim(a)#53
meta_fc5_gene= rownames(a)

range(exp[meta_fc5_gene,])#-1.4 4.4
pheatmap(exp[meta_fc5_gene,], cluster_rows = F, cluster_cols = F, annotation_col=anno,
         breaks=seq(-2,2,length.out=100), show_colnames = F)


range(exp[immune_genes,])#-2 6
length(unique(immune_genes))
pheatmap(exp[immune_genes,], cluster_rows = F, cluster_cols = F, annotation_col=anno,
         breaks=seq(-2,2,length.out=100), show_colnames = F, show_rownames = F)


range(exp[signal_genes,])#-3.3 6.6
length(unique(signal_genes))
pheatmap(exp[signal_genes,], cluster_rows = F, cluster_cols = F, annotation_col=anno,
         breaks=seq(-3,3,length.out=100), show_colnames = F, show_rownames = F)

write.table(anno$type,file='./data/anno.txt')
