library(limma)

preprocessESET <- function(eset){
  exprs(eset)=apply(exprs(eset),2,as.numeric)
  message('your expression data dimension is ')
  cat(length(rownames(exprs(eset))),length(colnames(exprs(eset))),'\n')
  cnt.zero=(rowSums(exprs(eset)==0))
  cnt.zero=as.data.frame(cnt.zero)
  
  cnt_exp=cbind(exprs(eset),cnt.zero)
  
  for (cnt in cnt_exp$'cnt.zero'){
    
    if(cnt==length(sampleNames(eset))){
      message('yor expression data has all 0 value..')
      eset=eset[-which(cnt_exp$'cnt.zero'==cnt),]

      message('tada!i removed it *----*','\n')
      message('Now, your expression data dimension is ')
      cat(length(rownames(exprs(eset))),length(colnames(exprs(eset))),'\n')
      break
    }
    
  }

  if (max(exprs(eset))>20){
    message('your expression data needs log2 transformation..')
    exprs(eset)=log2(exprs(eset)+1)
    cat(min(exprs(eset)),max(exprs(eset)),'\n')
    message('tada!i changed it *----*')

  }
  else{
    message('your data dosen\'t need log2 transformation')
  }
  
  
  message('i will use limma package normalizationfunction')
  exprs(eset)=normalizeQuantiles(exprs(eset))
  boxplot(exprs(eset))
}

#example
preprocessESET(eset)
