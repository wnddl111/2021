setwd('D:/')
# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

if (!("GSVA" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("GSVA", update = FALSE)
}

if (!("qusage" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("qusage", update = FALSE)
}

if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}
library(DESeq2)
library(qusage)
library(GSVA)
library(org.Hs.eg.db)
library(magrittr)

vst = read.table(file='./data/F0405_for_heatmap_0.3.txt', sep='\t')
col = read.table(file='./data2/my_data_colname_90.txt')
head(vst)
dim(vst) #18008
head(col)
colnames(vst)[1]='gene'
colnames(vst)[2:91] = col$V1
gene = vst[,1]
rownames(vst) =gene
#vst= vst[,-1]
#89명으로 만들기 61번 u3zg환자
#vst = vst[,-61]
colnames(vst)
vst =as.matrix(vst)

write.table(vst,file='./data/89_F0405_for_heatmap_0.3.txt',sep='\t')


#install.packages('msigdbr')
hallmark_gene_sets <- msigdbr::msigdbr(
  species='Homo sapiens', #바꿀 수 있음
  category ='H' #only hallmark gene sets
)
head(hallmark_gene_sets)

hallmarks_list <- split(
  hallmark_gene_sets$gene_symbol, #이걸 나눠서 list로 저장
  hallmark_gene_sets$gs_name #이걸 기준으로
)

head(filtered_mapped_matrix)

gsva_results <- gsva(
  vst,
  hallmarks_list,
  method='gsva',
  kcdf='Gaussian',
  min.sz=15,
  max.sz=500,
  mx.diff=TRUE,
  verbose=F
)
#만약 count data면 Poisson을 하면 되고 log2 취한거면 Gaussian

head(gsva_results[,1:10])

gsva_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column('pathway') %>%
  readr::write_tsv('./results/CRPC_gsva_results.tsv')

type = c(rep('T',36),rep('N',53))
sample =col$V1[-61]

feature = as.data.frame(type)
rownames(feature) = sample

write.table(feature,'./data/89_anno.txt',sep = '\t')

pathway_heatmap <- pheatmap(gsva_results, 
                            annotation_col = feature,
                            show_colnames=F, fontsize_row=6,breaks=seq(-0.4,0.4,length.out=100))
pathway_heatmap

png('./plots/CRPC_heatmap.png', width=1000, height=800)
pathway_heatmap
dev.off()

gsva_results


