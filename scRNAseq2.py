import sys
import os
import csv

#python3 scRNAseq2.py 1058B-10xP /data/juyoung/rawdata/ /data/juyoung/refdata-gex-GRCh38-2020-A /data/juyoung/scRNAseq/test2/

wdir ='/data/juyoung/scRNAseq'
os.chdir(wdir)

sample_names = sys.argv[1]
fastq_dir= sys.argv[2]
transcriptome_dir=sys.argv[3] 
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct#:~:text=transcriptome.%20From%20the-,download%20page,-for%20the%20FASTQ
new_dir = sys.argv[4]

os.makedirs(new_dir,exist_ok=True)
os.chdir(new_dir)
wdir=os.getcwd()
os.makedirs('./run_cellranger_count',exist_ok=True)

sample_names = sample_names.strip('[').strip(']').split(',')


os.chdir('./run_cellranger_count')

count_dir = os.getcwd()


for sample in sample_names:
  ids = sample.split('-')
Id = 'run_count_'+ids[0]+ids[1][-1]
print(Id)
os.system('cellranger count --id '+Id +' --fastqs='+fastq_dir+' --sample='+sample+' --transcriptome='+transcriptome_dir)



