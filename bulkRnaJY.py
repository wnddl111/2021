#############
#2021.08.27
#.gz/.fastq -> .gz만 돌아가게수정
#원본은 .fastq를 위한건데 잘 돌아갔음
#.gz들어가니까 if문 for문이 추가되는데 좀 번거로운 것 같기도 하고 앞으로 gz만 써야 용량 관리할 때 효율적일 것 같아서 gz만 되는 걸로 바꿈
##############
import os 
import sys
#from contextlib import contextmanager 필요없었던 것 같음 
import zipfile

#초기설정 
#ex)E010T_1.fastq, E010T2_1.fastq .... (gz만 돌아가게 수정 )
#이때 꼭 1과2의 순서를 맞춰서 입력해야한다. 
fastq1=list(map(str,sys.argv[1].strip().split(','))) 
fastq2=list(map(str,sys.argv[2].strip().split(',')))

#project별로 dir가 생성되기 때문에 참고해서 진행한다 
project_name= sys.argv[3]

#돌리기전에 rsem,star에 대한 index 파일은 이미 생성돼 있어야한다
#star는 ref밑에 있고 rsem은 따로있음 
wdir='/data/juyoung/bulkRNAseq/'
genome_ref_dir = '/data/juyoung/ref/' 
rsem_ref_dir='/data/juyoung/rsem/'
raw_dir='/data/juyoung/rawdata/'

#fastaq 파일이 match되지 않을 때 - > 시스템종료 
if (len(fastq1) != len(fastq2)):
        os.system('echo fastq file is not matched')
        sys.exit()
        
for i in range(len(fastq1)):
        fastq1[i] = fastq1[i][:-3]
        fastq2[i] = fastq2[i][:-3]
        
#make directory
#trimming 전후로 qc를 체크함 
os.chdir(wdir)
os.makedirs('./'+project_name,exist_ok=True)
os.chdir('./'+project_name)
path_list=['./qc','./qc2', './trimming', './mapping_sorting', 'counting']
for path in path_list:
        os.makedirs(path, exist_ok=True)

cdir = os.getcwd()

#qc
for i in range(len(fastq1)):
        os.system('gzip -d '+raw_dir+fastq1[i]+'.gz '+raw_dir+fastq2[i]+'.gz')
        os.system('fastqc -o '+cdir+'/qc/ -f fastq '+raw_dir+fastq1[i]+' '+raw_dir+fastq2[i])
        
os.chdir('./qc')#easy to move os.chdir command 
cdir = os.getcwd()

#fastqc 결과가 모여있는 폴더에서 _fastqc.zip파일을 열고 그 안에서 summary.txt에 해당하는 파일을 연다. 그 후 안에 있는 내용들을 모아 하나의 파일로 정리한다
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
        
with open('all_summary.txt', 'wb') as f: #summary.txt 가 바이너리 파일이라서 쓸 때도 이렇게 써야 한다 
        for content in all_summary:
                f.write(content)
        os.system('echo "./qc/all_summary.txt -> qc results of all fastqc file"')


#trimming        
for i in range(len(fastq1)):
        project = project_name.replace('#',str(i+1)) # '#'을 통해서 1번 2번 3번...환자를 구분한다 
        os.system('cutadapt -a AGATCGGAAGAGC -g AGATCGGAAGAGC -q 30 -m 20 -o '+wdir+project_name+'/trimming/trimmed_'+fastq1[i]+' -p '+wdir+project_name+'/trimming/trimmed_'+fastq2[i]+' '+raw_dir+fastq1[i]+' '+raw_dir+fastq2[i])
        os.system('fastqc -o '+wdir+project_name+'/qc2/ -f fastq '+wdir+project_name+'/trimming/trimmed_'+fastq1[i]+' '+wdir+project_name+'/trimming/trimmed_'+fastq2[i])
        os.system('gzip '+raw_dir+fastq1[i]+' '+raw_dir+fastq2[i])

#trimming 후 qc 진행 
os.chdir(wdir+project_name+'/qc2/')# '/'게 있어야 이동한대 
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

#star로 mapping rsem으로 counting 
for i in range(len(fastq1)):
        project= project_name.replace('#',str(i+1))
        #mapping
        os.system('STAR --runMode alignReads --runThreadN 16 --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 265006 --genomeDir '+genome_ref_dir+' --readFilesIn '+wdir+project_name+'/trimming/trimmed_'+fastq1[i]+' '+wdir+project_name+'/trimming/trimmed_'+fastq2[i]+' --outFileNamePrefix '+wdir+project_name+'/mapping_sorting/star_mapsort_'+project+'_'+' --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM')
        #counting
        os.system('rsem-calculate-expression -p 8 --alignments --paired-end --strandedness reverse --no-bam-output '+wdir+project_name+'/mapping_sorting/star_mapsort_'+project+'_Aligned.toTranscriptome.out.bam '+rsem_ref_dir+' '+wdir+project_name+'/counting/count_'+project)
        #최종 파일 이름 r에서 사용하기 쉽게 변경 
        os.chdir(wdir+project_name+'/counting/')
        result_name = 'count_'+project+'.genes.results'
        os.system('mv '+result_name+' '+result_name+'.txt')
        
        
        
'''설명 
python3 bulkRnaJY.py fastq1-1,fastq2-1(.gz) fastq1-2,fastq2-2(.gz) 프로젝트명_#(#은 필수/1번환자 2번환자 구분하기위함)
######STAR INDEX########
STAR --runThreadN 30 —-runMode genomeGenerate —genmoeDir /data/juyoung/ref —genomeFastaFiles /data/juyoung/ref/gencode~.fa —-sjdbGTFfile /data/juyoung/ref/gencode~.gtf —-sjdbOverhang 100

--runMode: STAR 프로그램의 실행모드 (인덱스 생성은 “genomeGenerate”로 설정함)
--genomeDir: 참조 유전체 데이터 파일의 인덱스를 생성하고 저장할 디렉토리 
--genomeFastaFiles: 참조 유전체 서열 데이터 파일 (FASTA 포맷)
--sjdbGTFfile: gtf파일넣기 -> 왠만하면 무조건 하는 걸 추천한대 
--sjdbOverhang : genomic sequence의 길이 , max(ReadLength)-1 로 구하면되는데 대부분의 case는 100에서 잘 돌아감 

*고대사랩에서는 genmome ref와 star ref구분을 안해줬지만 다음에는 해주는 게 좋을듯하다 

######RSEM INDEX########
GENOME=~/ref/GRCh38/GRCh38.primary_assembly.genome.fa 
GENCODE=~/ref/GRCh38/gencode.v36.primary_assembly.annotation.gtf
OUTDIR=~/ref/rsem/GRCh38_GENCODEv36
rsem-prepare-reference --gtf $GENCODE --bowtie $GENOME $OUTDIR

######GENOME DATA######
Gene code에서 download함 
https://www.gencodegenes.org/human/

#####fastqc###########
conda로 다운 받음
outputfile name을 지정하지 않아도 된다! input파일 이름을 고대로 사용하고 확장자명만 바뀐다

-o : 결과물을 생성하고 저장할 디렉토리
-f : 입력 데이터파일의 종류(bam, sam, fastq 지원) 

결과로 zip파일과 html이 나옴
zip 안에 summary.txt를 보면 품질검사 단계에 따른 최종결과가 나옴(pass/warn/fail)

#####sickle#########

se/pe: 원시 데이터의 타입이 single-end 이면 se, paired-end 이면 pe 로 설정
-f: 입력 데이터의 첫번째 read 파일
-r: 입력 데이터의 두번째 read 파일 (원시 데이터의 타입이 paired-end 인 경우에만 해당됨) 
-t: base quality 인코딩 타입에 따라 설정 (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7),
sanger (CASAVA >=1.8))
-o: 첫번째 read 파일의 결과물을 생성하고 저장할 파일 이름
-p: 두번째 read 파일의 결과물을 생성하고 저장할 파일 이름 (원시 데이터의 타입이
paired-end 인 경우에만 해당됨)
-s: single-end 로 변환되는 결과물을 생성하고 저장할 파일 이름 (원시 데이터의 타입이
paired-end 인 경우에만 해당됨)
-q: window 내 base quality score 의 평균 기준 -l: 최소 read 길이 기준

#####star##########
STAR --runMode alignReads --runThreadN 16 --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 265006 --genomeDir '+genome_ref_dir+' --readFilesIn '+wdir+project_name+'/trimming/trimmed_'+fastq1[i]+' '+wdir+project_name+'/trimming/trimmed_'+fastq2[i]+' --outFileNamePrefix '+wdir+project_name+'/mapping_sorting/star_mapsort_'+project+'_'+' --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

--runMode: STAR 프로그램의 실행모드 (맵핑은 “alignReads”로 설정함)
--runThreadN: 맵핑에 사용할 thread 개수
--outFilterMultimapNmax: 결과물 (BAM 포맷의 맵핑된 데이터)에 기록 가능한 각 read 별 최대
맵핑 허용치 (허용치를 초과할 경우 맵핑 되지 않는 것으로 판단함)
--alignIntronMin: 맵핑 허용 최소 인트론 길이 (본 내용에서는 알려진 인간 유전자의 인트론
길이 분포로부터 0.01% 구간을 설정함)
--alignIntronMax: 맵핑 허용 최대 인트론 길이 (본 내용에서는 알려진 인간 유전자의 인트론
길이 분포로부터 99.9% 구간을 설정함)
--genomeDir: 맵핑에 사용할 참조 유전체 인덱스가 저장된 디렉토리
--readFilesIn: 품질 관리된 입력 데이터 파일 (입력 데이터의 타입이 paired-end 인 경우는
공백으로 첫번째 read 와 두번째 read 파일을 구분하여 입력)
--outSAMtype: 출력 데이터 파일의 종류 (SAM: SAM 파일, BAM Unsorted: 포지션으로 정돈되지
않은 BAM 파일, BAM SortedByCoordinate: 포지션으로 정돈된 BAM 파일) --outFileNamePrefix: 맵핑된 결과물을 생성하고 저장할 파일 이름의 접두사를 지정
--quantMode: With --quantMode TranscriptomeSAM GeneCounts, and get both the Aligned.toTranscriptome.out.bam and ReadsPerGene.out.tab outputs.

**aligned transcriptome 결과파일이 rsem을 돌리기 위해서는 반드시 필요함 


####rsem########
THREADS=8
BAM=~/data/project/bam/star.toTranscriptome.out.bam
REF=~/ref/rsem/GRCh38_GENCODEv36
OUT_PREFIX=~/data/project/rsem/quant
rsem-calculate-expression \
-p $THREADS \
--alignments \
--paired-end \
--strandedness reverse \ #strand specific 프로토콜을 사용했는지 
--no-bam-output \ #bam 파일 만들지 말라는 것 
$BAM \
$REF \
$OUT_PREFIX
-p : 8로 고정 -> 이걸 추천하고 있음

최종 
'''








