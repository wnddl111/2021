import os
import sys
from contextlib import contextmanager
import zipfile
import subprocess



fastq1=list(map(str,sys.argv[1].strip().split(','))) #input type is string which pretend list, ex)E010T_1.fastq, E010T2_1.fastq ....
fastq2=list(map(str,sys.argv[2].strip().split(',')))

project_name= sys.argv[3]

wdir='/data/juyoung/bulkRNAseq/'
genome_ref_dir = '/data/juyoung/ref/' 
rsem_ref_dir='/data/juyoung/rsem/'
raw_dir='/data/juyoung/rawdata/'

#match fastq error

if (len(fastq1) != len(fastq2)):
        os.system('echo fastq file is not matched')
        sys.exit()

#make directory
os.chdir(wdir)
os.makedirs('./'+project_name,exist_ok=True)
os.chdir('./'+project_name)
path_list=['./qc','./qc2', './trimming', './mapping_sorting', 'counting']
for path in path_list:
        os.makedirs(path, exist_ok=True)

cdir = os.getcwd()

#fastq file


for i in range(len(fastq1)):
        os.system('fastqc -o '+cdir+'/qc/ -f fastq '+raw_dir+fastq1[i]+' '+raw_dir+fastq2[i])


os.chdir('./qc')#easy to move os.chdir command 
cdir = os.getcwd()

files = [file for file in os.listdir(cdir) if file.endswith('_fastqc.zip')]
all_summary = []

for file in files:
        archive = zipfile.ZipFile(file,'r')
        members = archive.namelist()
        
        fname = [member for member in members if 'summary.txt' in member][0]
        data=archive.open(fname)
        
        for line in data:
                all_summary.append(line)
        data.close()
        archive.close()
        
with open('all_summary.txt', 'wb') as f: #summary.txt is binary file ;(
        for content in all_summary:
                f.write(content)
        os.system('echo "./qc/all_summary.txt -> qc results of all fastqc file"')

for i in range(len(fastq1)):
        project = project_name.replace('#',str(i+1)) # revision is required 
        #trimming
        os.system('sickle pe -f '+raw_dir+fastq1[i]+' -r '+raw_dir+fastq2[i]+' -t sanger -o '+wdir+project_name+'/trimming/trimmed_'+fastq1[i]+' -p '+wdir+project_name+'/trimming/trimmed_'+fastq2[i]+' -s '+wdir+project_name+'/trimming/single_trimmed_'+project+' -q 20 -l 20')
        
        os.system('fastqc -o '+wdir+project_name+'/qc2/ -f fastq '+wdir+project_name+'/trimming/trimmed_'+fastq1[i]+' '+wdir+project_name+'/trimming/trimmed_'+fastq2[i])

os.chdir(wdir+project_name+'/qc2/')# '/' is must required 
cdir=os.getcwd()

files = [file for file in os.listdir(cdir) if file.endswith('_fastqc.zip')]
all_summary = []

for file in files:
        archive = zipfile.ZipFile(file,'r')
        members = archive.namelist()
        fname = [member for member in members if 'summary.txt' in member][0]
        data=archive.open(fname)

        for line in data:
                all_summary.append(line)
        data.close()
        archive.close()
with open('after_trimming_all_summary.txt','wb') as f:
        for content in all_summary:
                f.write(content)
        os.system('echo "./qc2/after_trimming_all_summary.txt -> qc results of all trimmed_fasta_qc"')



for i in range(len(fastq1)):
        project= project_name.replace('#',str(i+1))
        #mapping
        os.system('STAR --runMode alignReads --runThreadN 16 --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 265006 --genomeDir '+genome_ref_dir+' --readFilesIn '+wdir+project_name+'/trimming/trimmed_'+fastq1[i]+' '+wdir+project_name+'/trimming/trimmed_'+fastq2[i]+' --outFileNamePrefix '+wdir+project_name+'/mapping_sorting/star_mapsort_'+project+'_'+' --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM')
        #counting
        os.system('rsem-calculate-expression -p 8 --alignments --paired-end --strandedness reverse --no-bam-output '+wdir+project_name+'/mapping_sorting/star_mapsort_'+project+'_Aligned.toTranscriptome.out.bam '+rsem_ref_dir+' '+wdir+project_name+'/counting/count_'+project)














