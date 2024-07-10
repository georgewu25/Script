#!/bin/usr/Rscript

library(data.table)

OR3<-function(querys,terms_list,background,min.term.size=0,max.term.size=Inf,overlap_column=TRUE,verbose=FALSE){
  if(is.list(querys)){
    dt<-Reduce(rbind,lapply(names(querys),
                            function(q)OR3(querys = querys[[q]],
                                           terms_list = terms_list,
                                           background = background,
                                           min.term.size = min.term.size,
                                           max.term.size = max.term.size,
                                           overlap_column = overlap_column,
                                           verbose = verbose )[,query:=q]))
    
    return(dt[,query.:=query][,.SD,.SDcols=c(ncol(dt),1:(ncol(dt)-1))])
  }else{
    #Common elements between query and background
    queryf<-intersect(querys,background)
    #Commom elements between each list and background
    terms_listf<-lapply(terms_list,function(x)intersect(x,background))
    
    #Restrain on term list size
    res_or<-data.table(term=names(terms_listf),term.size=sapply(terms_listf,length))
    res_or<-res_or[term.size<=max.term.size]
    n_terms<-nrow(res_or)
    if(verbose)message(length(terms_listf)-n_terms, " terms were filtered due to term.size above the limit of ",max.term.size," genes")
    res_or<-res_or[term.size>=min.term.size]
    n_terms<-nrow(res_or)
    if(verbose)message(n_terms-nrow(res_or), " terms were filtered due to term.size below the limit of ",min.term.size," genes")
    
    res_or[,n.query:=length(queryf)]
    res_or[,n.overlap:=sum(queryf%in%terms_listf[[term]]),by="term"]
    
    res_or[,pct.query.overlap:=n.overlap/n.query]
    res_or[,precision:=pct.query.overlap]
    
    res_or[,pct.term.overlap:=n.overlap/term.size]
    
    res_or[,background_size:=length(background)]
    
    res_or[,pct.term.background:=term.size/background_size] 
    
    #Hypergeometric test
    res_or[,pval:=phyper(q=n.overlap-1, 
                         m=term.size, 
                         n=background_size-term.size, 
                         k=n.query, 
                         lower.tail=FALSE),
           by="term"]
    
    #Adjusted p value
    res_or[,padj:=p.adjust(pval,method = 'BH')]
    res_or[,fold.enrichment:=pct.query.overlap/pct.term.background]
    if(overlap_column==TRUE){
      res_or[,genes.overlap:=paste(queryf[queryf%in%terms_listf[[term]]],collapse="|"),by="term"]
    }
    if(verbose)message(nrow(res_or[padj<0.05])," terms enriched in your genes of interest with padj<0.05")
    return(res_or)
  }
}



