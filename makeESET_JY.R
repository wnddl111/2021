
library(GEOquery)
library(stringr)
library(dplyr)
library("biomaRt")

processEXP<- function(geo_file,matrix){ #matrix가 2이상일 때 씀 
  intersect_row=rownames(exprs(geo_file[[1]]))
  for (i in 2:matrix){
    intersect_row=intersect(intersect_row, rownames(exprs(geo_file[[i]])))
  }
  message('\n your matrix num is larger than 1, so need extracting common row')
  message('common row length : ',length(intersect_row))
  
  geo_exp=as.data.frame(exprs(geo_file[[1]])[intersect_row,])
  for (i in 2:matrix){
    geo_exp=cbind(geo_exp,exprs(geo_file[[i]])[intersect_row,])
  }
  return(geo_exp)
}


processPDATA <-function(geo_file,matrix){
  pdata=pData(geo_file[[1]])
  message('if colnames of pdatas are different each other, it dosen\'t work')
  message('so please use same columns..')
  for ( i in 2:matrix){
    pdata=rbind(pdata,pData(geo_file[[i]]))
  }
  return(as.data.frame(pdata))
}

processFDATA <- function(geo_file,matrix){
  f_row=rownames(fData(geo_file[[1]]))
  for ( i in 2:matrix){
    f_row=intersect(f_row,rownames(fData(geo_file[[i]])))
  }
  fdata=fData(geo_file[[1]])[f_row,]
  return(fdata)
}

changeGENEsymbol <-function(fdata,matrix){
  answer=''
  while(answer != 'y'&& answer !='n' && answer != '='){
    message('do you have any column name which represents gene symbol?
                  (if rownames of fdata are gene symbol, just write "=" ) 
                  ex) y or n or = ')
    answer=readline(': ')
  }
  if (answer=='y'){
    g=readline('write that column name : ')
    genesymbol=sapply(1:length(rownames(fdata)), function(a) trimws(strsplit(fdata[a,g],"///")[[1]][1]))
    return(genesymbol)
  }else if(answer == 'n'){
    message('if you dont have any gene_symbol, write your type of id...
             ex) entrezgene_id (write a replaceable real column name!!!)
             ')
    col=readline(': ')
    message('For changing your id, Choose one matches with yours.. 
    ensembl_gene_id
    ensembl_gene_id_version
    ensembl_transcript_id
    ensembl_peptide_id
    hgnc_symbol
    entrezgene
    entrezgene_id')
    b=readline(': ')
    
    i=0
    while(i<5){
      hsmart=useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
      
      mygenes=fdata$col
      
      mapping=getBM(
        attributes = c(b, hgnc_symbol), 
        filters = b,
        values = mygenes,
        mart = hsmart
      )
      i=i+1
      if (length(mapping)>1){
        return(mapping)
        
        break
      }
      else{
        cat('we try to get gene symbol, but fails..
          because of biomart server, so try again..!')
      }
    }
    
  } 
  else {
    genesymbol=fdata[,g]
    return(genesymbol)
  }

  
  
}
 
makeESET<- function(geo_num,res,con,matrix){
  cat('。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。  。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・ ❅・。・。 ゜・ 。 。 ゜M a  k ❅e。 ❅・。E・S。 E゜T。。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅')
  cat('\n')
  geo_file=getGEO(geo_num)
  
  res_num=strsplit(res,',')
  con_num=strsplit(con,',')
  
  
  if (matrix >1){
    geo_exp=processEXP(geo_file,matrix)
    message('expression dimension(just finished "common row" process) : ',length(rownames(geo_exp)),' ',length(colnames(geo_exp)))
    pdata=processPDATA(geo_file,matrix)
    fdata=processFDATA(geo_file,matrix)
  }
  
  else{
    geo_exp=as.data.frame(exprs(geo_file[[1]]))
    message('expression dimension(just use one array case) : ', length(rownames(geo_exp)),' ',length(colnames(geo_exp)))
    pdata=pData(geo_file[[1]])
    fdata=fData(geo_file[[1]])
  }
  
  
  
  for (i in 1:length(colnames(geo_exp))){
    colnames(geo_exp)[i]=str_sub(colnames(geo_exp)[i],-2,-1)
  }
  for (i in 1:length(rownames(pdata))){
    rownames(pdata)[i]=str_sub(rownames(pdata)[i],-2,-1)
  }

  res_sub=select(geo_exp,unlist(intersect(colnames(geo_exp), unlist(res_num))))
  con_sub=select(geo_exp,unlist(intersect(colnames(geo_exp), unlist(con_num))))
  
  res_pdata=pdata[unlist(intersect(rownames(pdata), unlist(res_num))),]
  con_pdata=pdata[unlist(intersect(rownames(pdata), unlist(con_num))),]
  
  message('\n')
  cat('❅。 ❅ 。 。 ゜。・゜・ 。 Final GROUP DIM❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。')
  cat('\n')
  message('Res dim is ',  length(rownames(res_sub)),' ',length(colnames(res_sub)))
  message('Con dim is ',  length(rownames(con_sub)),' ',length(colnames(con_sub)))
  
  exp=cbind(con_sub,res_sub)
  pdata=rbind(con_pdata,res_pdata)
  
  genesymbol=changeGENEsymbol(fdata,matrix) 
  
  a_exp=aggregate(exp,list(genesymbol),mean,na.rm=T)
  fdata=fdata[rownames(a_exp),]
  g=a_exp[,1]
  
  a_exp=a_exp[,-1]
  rownames(a_exp)=g
  rownames(fdata)=g
  
  
  featureData=new("AnnotatedDataFrame", data=fdata)
  phenoData = new("AnnotatedDataFrame", data=pdata)
  eset= new('ExpressionSet', exprs=a_exp, phenoData=phenoData, featureData=featureData)
  
  message('WAIT! final eset alreay ready!')
  message('Expression dim is ',  length(rownames(exprs(eset))),' ',length(colnames(exprs(eset))))
  message('Pdata dim is ',  length(rownames(pData(eset))),' ',length(colnames(pData(eset))))
  message('Fdata dim is ',  length(rownames(fData(eset))),' ',length(colnames(fData(eset))))
  
  cat('。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅・。・。 ゜・ 。 。 ゜。・゜ ❅。 ❅')
  return(eset)
  
}

res_num='21,22,26'
con_num='19,20,25'
#eset=makeESET('GSE158494',res_num,con_num,1)
eset=makeESET('GSE36135',res_num,con_num,2)

res_num=readline('res sample의 끝 두자리를 쓰세요(여러개라면 , 로 구분) : ')
con_num=readline('con sample의 끝 두자리를 쓰세요(여러개라면 , 로 구분) : ')
matrix_num=readline('series matrix 개수를 쓰세요 : ')
#19 20 25 sen
#21 22  26 res
GSE36135=getGEO('GSE36135')
exprs(GSE36135$`GSE36135-GPL570_series_matrix.txt.gz`)
GSE
length(intersect(rownames(exprs(GSE36135$`GSE36135-GPL570_series_matrix.txt.gz`)),
          rownames(exprs(GSE36135$`GSE36135-GPL571_series_matrix.txt.gz`))))

strsplit('15,16,17',',')




a <-function(i){
  print('a야')
  return(i)
}
(fData(GSE36135$`GSE36135-GPL570_series_matrix.txt.gz`))$'Gene Symbol'

b<-function(i){
  a(i)
}

b(1)





