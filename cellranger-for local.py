#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import os
import csv

wdir ='/data/juyoung/scRNAseq'
os.chdir(wdir)

sample_names = sys.argv[1]
fastq_dir= sys.argv[2]
transcriptome_dir=sys.argv[3] 

new_dir = sys.argv[4]

os.makedirs(new_dir,exist_ok=True) #새로운 프로젝트에대한 dir생성
os.chdir(new_dir)
wdir=os.getcwd()
os.makedirs('./run_cellranger_count',exist_ok=True)#count파일을 담을 dir생성

sample_names = sample_names.strip('[').strip(']').split(',') #[,,]이런식으로 입력을 받기 때문에 필요함

os.chdir('./run_cellranger_count')

count_dir = os.getcwd()

for sample in sample_names:
    ids = sample.split('-')#[1058B,10xp]
    Id = 'run_count_'+ids[0]+ids[1][-1]
    os.system('cellranger count --id '+Id +' --fastqs='+fastq_dir+' --sample='+sample+' --transcriptome='+transcriptome_dir)
    
'''    
> example
python3 scRNAseq2.py [1058B-10xP,1058B-10xT] /data/juyoung/rawdata/ /data/juyoung/refdata-gex-GRCh38-2020-A /data/juyoung/scRNAseq/ test2/

> cellranger 공식 홈페이지의 tutorial대로 만듬
> cellranger option
--id : directory가 만들어질건데 거기에 쓸 이름 
--fastqs : raw data가 있는 위치
--sample : fastq파일의 이름의 앞쪽 글자! 다른 샘플파일과 구분할 수 있는 것까지 쓰면 됨 
--transcriptome : cell ranger reference ! 다운은 아래링크에서 받음
* #https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct#:~:text=transcriptome.%20From%20the-,download%20page,-for%20the%20FASTQ


> 최종파일은
프로젝트명 > id dir > outs > filterd ~ > 3개의 파일이 있음 

이 파일은 서버에서 한꺼번에 돌리고 로컬r에서 합치기위한것 

