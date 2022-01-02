library(Seurat)

sample_list_gastric=read.table(file='/data/juyoung/scRNAseq/immunotherapy_gastric/run_cellranger_count/for_R_samplelist_name.txt')
sample_list_liver=read.table(file='/data/juyoung/scRNAseq/immunotherapy_liver/run_cellranger_count/for_R_samplelist_name.txt')

sample_list_gastric <- as.list(sample_list_gastric$V1)
sample_list_liver <- as.list(sample_list_liver$V1)

dir_name='run_cellranger_count'
project_name='immunotherapy_gastric'

setwd('/data/juyoung/scRNAseq/')

for (i in sample_list_gastric){
  assign(paste0(i,'.data'),Read10X(data.dir=paste0('./',project_name,'/',dir_name,'/run_count_',i,'/outs/filtered_feature_bc_matrix')))
  assign(i, CreateSeuratObject(counts=get(paste0(i,'.data')),project = i,min.cells = 3, min.features = 200))
  seurat=get(i)
  saveRDS(seurat,file=paste0('/data/juyoung/scRNAseq/for_local/gastric_',i,'.rds'))
}

project_name='immunotherapy_liver'
for (i in sample_list_liver){
  print(i)
  assign(paste0(i,'.data'),Read10X(data.dir=paste0('./',project_name,'/',dir_name,'/run_count_',i,'/filtered_feature_bc_matrix')))
  assign(i, CreateSeuratObject(counts=get(paste0(i,'.data')),project = i))
  seurat=get(i)
  saveRDS(seurat,file=paste0('/data/juyoung/scRNAseq/for_local/liver_',i,'.rds'))
}
















