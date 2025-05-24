
library(igraph)
library(ggraph)
library(stringr)
library(ggplot2)
library(data.table)
library(rlang, lib = "/projectnb/tcwlab/LabMember/mwu/R/")

#Useful Functions####
#get sim matrix

LeadingEdges<-function(res_fgsea){
  if(all(c('term','n.overlap','genes.overlap')%in%colnames(res_fgsea))){
    l_genes<-str_extract_all(res_fgsea$genes.overlap,'[A-Za-z0-9]+')
    l_genes<-lapply(l_genes, function(x)x[x!='c'])
    names(l_genes)<-res_fgsea$pathway
    return(l_genes)
    
  }else{
    l_genes<-str_extract_all(res_fgsea$leadingEdge,'[A-Za-z0-9]+')
    l_genes<-lapply(l_genes, function(x)x[x!='c'])
    names(l_genes)<-res_fgsea$pathway
    return(l_genes)
  }
  
  
}

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}


get_similarity_matrix <- function(leading_edge_list) {
  
  n <- length(leading_edge_list)
  ids<-names(leading_edge_list)
  w <- matrix(NA, nrow=n, ncol=n)
  colnames(w) <- rownames(w) <- names(leading_edge_list)
  
  for (i in seq_len(n-1)) {
    for (j in (i+1):n) {
      w[i,j] <- overlap_ratio(leading_edge_list[ids[i]], leading_edge_list[ids[j]])
    }
  }
  return(w)
}


# get graph of sim
get_igraph <- function(res_fgsea, simmat,leading_edge_list,
                       pathway_names=NULL, min_edge=0.2, col.var=NULL) {
  if(any(duplicated(res_fgsea$pathway)))stop('error: duplicated pathways')
  
  if(is.null(pathway_names)){
    pathway_names<-res_fgsea$pathway
  }
  wd <- reshape2::melt(simmat[pathway_names,pathway_names])
  wd <- wd[wd[,1] != wd[,2],]
  # remove NA
  wd <- wd[!is.na(wd[,3]),]
  
  g <- graph.data.frame(wd[, -3], directed=FALSE)
  E(g)$width <- sqrt(wd[, 3] * 5) 
  
  
  
  # Use similarity as the weight(length) of an edge
  E(g)$weight <- wd[, 3]
  g <- delete.edges(g, E(g)[wd[, 3] < min_edge])
  
  res_fgseaf<-res_fgsea[V(g)$name,on='pathway']
  #idx <- unlist(sapply(V(g)$name, function(x) which(x == res_fgseaf$pathway)))
  cnt <- sapply(leading_edge_list, length)
  
  V(g)$size <- cnt[V(g)$name]
  
  if(!is.null(col.var)){
    colVar <- res_fgseaf[V(g)$name, on='pathway'][[col.var]]
    V(g)$colvar <- colVar
    
  }
  
  return(g)
}


#plot the graphs
add_category_nodes <- function(p,col.var,cols=c('blue','white','red'),cols_lims=NULL) {
  
  if(!any('discrete'%in%cols)){
    locol=cols[1]
    if(length(cols)==3)midcol=cols[2]
    
    hicol=ifelse(length(cols)==3,cols[3],cols[2])
    
  }
  
  p<-p + ggnewscale::new_scale_fill() +geom_point(shape = 21, aes_(x =~ x, y =~ y, fill =~ colvar,
                                                                   size =~ size)) +
    scale_size_continuous(name = "number of genes",
                          range = c(3, 8) )
  
  if('discrete'%in%cols){
    p<-p+scale_fill_discrete(name=col.var ) 
    
  }else if(length(cols)==3){
    p<-p+scale_fill_gradient2(low = locol,mid=midcol, high = hicol,name=col.var,
                              guide = guide_colorbar(),
                              limits=cols_lims,midpoint = 0) 
    
  }else{
    p<-p+scale_fill_continuous(low = locol, high = hicol,name=col.var,
                               guide = guide_colorbar(),
                               limits=cols_lims)
    
    
  }
  
  p<-p+theme(legend.title = element_text(size = 10),
             legend.text  = element_text(size = 10)) +
    theme(panel.background = element_blank()) 
  
  return(p)
}
add_node_label <- function(p,label.size=label.size,max.overlaps=10) {
  
  p <- p + geom_node_text(aes_(label=~name),
                          size = label.size, repel=TRUE,
                          max.overlaps=max.overlaps)
  
  return(p)
}

#main function####

emmaplot<-function(res_fgsea,
                   pathway_names=NULL, 
                   col.var="NES",
                   show_pathway_of=NULL,
                   min_edge=0.2,
                   label.size=2.5,
                   cols=c('blue','white','red'),
                   max.overlaps=10,cols_lims=NULL){
  require('ggrepel')
  if(!'data.table'%in%class(res_fgsea)){
    res_fgsea<-data.table(res_fgsea)
  }
  
  if(!is.null(cols_lims)&length(cols)>2){
    if(!(any(cols_lims<0)&any(any(cols_lims>=0)))){
      cols=c('white','red')
    }
  }
  
  
  if(all(c('term','n.overlap')%in%colnames(res_fgsea))){
    res_fgsea[,pathway:=term]
    res_fgsea[,NES:=fold.enrichment]
    res_fgsea[,size:=n.overlap]
    
  }
  
  if(!'pathway'%in%colnames(res_fgsea)){
    stop('expected format: fgsea results or overrepresentation results (OR3) format')
  }
  
  
  if(all(table(res_fgsea[[col.var]])>1)){
    
    res_fgsea[,(col.var):=lapply(.SD,as.character),.SDcols=col.var]
    
    cols='discrete'
  }
  
  
  
  if(is.null(pathway_names))pathway_names=res_fgsea[order(pval)]$pathway
  
  lelist<-LeadingEdges(res_fgsea[pathway%in%pathway_names])
  
  if(!is.null(show_pathway_of)){
    
    lelist<-lelist[sapply(lelist, function(leadingedges)any(show_pathway_of%in%leadingedges))]
    if(length(lelist)>0){
      pathway_names<-names(lelist)
      
    }else{
      stop('This gene are not found in any leading edges of the given pathways')
    }
  }
  
  if(length(lelist)>1){
    
    simat<-get_similarity_matrix(lelist)
    
    g <- get_igraph(res_fgsea = res_fgsea,
                    pathway_names = pathway_names,
                    simmat = simat,
                    leading_edge_list = lelist,
                    min_edge = min_edge,
                    col.var = col.var
    )
    
    
    p <- ggraph(g, layout='nicely')
    
    p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                            colour='darkgrey')
    ## add dot
    p <- add_category_nodes(p = p,col.var =col.var,cols=cols,cols_lims=cols_lims)
    
    ## add node label
    
    p <- add_node_label(p = p,label.size=label.size,max.overlaps=max.overlaps)
    
  }else{
    p <- ggplot(res_fgsea[pathway%in%pathway_names][,x:=1][,y:=1],aes_string(x='x',y='x'))+
      geom_point(aes_string(size= "size",col=col.var))+
      geom_text_repel(aes(label=pathway))+
      scale_color_gradient2(low = cols[1],high = cols[length(cols)],
                            midpoint = 0,limits=c(-abs(as.numeric(as.vector(res_fgsea[pathway%in%pathway_names][,..col.var]))),
                                                  abs(as.numeric(as.vector(res_fgsea[pathway%in%pathway_names][,..col.var])))))+
      theme_graph()
    
  }
  
  
  if(!is.null(show_pathway_of)){
    if(length(show_pathway_of)>1){
      return(p+ggtitle(paste('Enriched pathways for selected genes')))
      
    }else{
      return(p+ggtitle(paste('Enriched pathways with', show_pathway_of)))
      
    }
    
  }else{
    return(p)
    
  }
}


#Annexe functions####

GetPathwaysLinks<-function(res_fgsea,
                           pathway_names=NULL,
                           show_pathway_of=NULL,
                           min_edge=0.2){
  require('ggrepel')
  
  if(is.null(pathway_names))pathway_names=res_fgsea[order(pval)]$pathway
  
  lelist<-LeadingEdges(res_fgsea[pathway%in%pathway_names])
  
  if(!is.null(show_pathway_of)){
    
    lelist<-lelist[sapply(lelist, function(leadingedges)any(show_pathway_of%in%leadingedges))]
    if(length(lelist)>0){
      pathway_names<-names(lelist)
      
    }else{
      stop('This gene are not found in any leading edges of the given pathways')
    }
  }
  
  if(length(lelist)>1){
    
    simat<-get_similarity_matrix(lelist)
    
    g <- get_igraph(res_fgsea = res_fgsea,
                    pathway_names = pathway_names,
                    simmat = simat,
                    leading_edge_list = lelist,
                    min_edge = min_edge
    )
    pathways_links<-data.table(get.edgelist(g))
    setnames(pathways_links,new = c('pathway1','pathway2'))
    pathways_links[,weight:=simat[pathway1,pathway2],by=c('pathway1','pathway2')]
    
    return(pathways_links)
    
  }else{
    stop('only one pathway gave')
    
  }
}


ClusterPathways<-function(x,resolution=1,method='louvain',min_edge=NULL,weights=NULL){
  require(igraph)
  
  
  if(!'weight'%in%colnames(x)){
    if(is.null(min_edge)){
      min_edge=0.2
    }
    x<-GetPathwaysLinks(x,min_edge = min_edge)
    
  }
  
  graph<-graph.data.frame(data.frame(x),directed = F)
  
  
  if(method=='louvain'){
    cl<-cluster_louvain(graph, weights = weights, resolution = resolution)
    
  }else{
    stop('no other method than louvain implemented yet')
  }
  
  
  return(data.table(pathway=cl$names,cluster=cl$membership))
}


GetRepresentativePathways<-function(res_fgsea,group.by=NULL,pathway_names=NULL,return_full=TRUE,padj.thr=0.05,resolution=1){
  
  res<-copy(res_fgsea)
  if(max(res$padj)>padj.thr){
    message('removing pathways not passing padj threshold ',padj.thr)
    message(nrow(res[padj<padj.thr]),'/',nrow(res),' pathways conserved')
    res<-res[padj<padj.thr]
    
  }
  
  if(!is.null(pathway_names))
    res<-res[pathway%in%pathway_names]
  
  if(is.null(group.by)){
    group.by='group'
    res<-res[,(group.by):=1]
  }
  
  if(length(group.by)>1){
    res<-res[,group:=apply(.SD,1,function(x)paste(x,collapse = '_')),.SDcols=group.by]
    group.by='group'
  }
  
  res<-rbindlist(lapply(unlist(unique(res[,..group.by])), function(g){
    message('finding representative pathways for ',g)
    res1<-unique(res[g,on=group.by][order(padj)],by='pathway')
    if(nrow(res1)>1){
      pathways_links<-GetPathwaysLinks(res1)
      
      pathways_cluster<-ClusterPathways(pathways_links,
                                        method='louvain',resolution = resolution,weights = NA)
      
      
      pathway_stat<-data.table(pathway=union(pathways_links$pathway1,pathways_links$pathway2))
      
      pathway_stat[,n.link:=nrow(pathways_links[pathway1==pathway|pathway2==pathway]),by='pathway']
      
      pathways_links<-merge(pathways_links,copy(pathways_cluster)[,pathway1:=pathway][,cluster1:=cluster][,.(pathway1,cluster1)])
      pathways_links<-merge(pathways_links,copy(pathways_cluster)[,pathway2:=pathway][,cluster2:=cluster][,.(pathway2,cluster2)],by='pathway2')
      pathway_stat<-merge(pathway_stat,pathways_cluster)
      pathway_stat[,n.link.in.cluster:=nrow(pathways_links[(pathway1==pathway&cluster2==cluster)|(pathway2==pathway&cluster1==cluster)]),by='pathway']
      
      res1<-merge(res1,pathway_stat,all.x = TRUE)
      
      res1[is.na(cluster),cluster:=0]
      
      #central one + the one with bigger NES, padj by cluster
      res1[,top.central:=rank(n.link.in.cluster)==max(rank(n.link.in.cluster)),by='cluster']
      res1[,top.NES:=rank(abs(NES))==max(rank(abs(NES))),by='cluster']
      res1[,top.pval:=rank(pval)==min(rank(pval)),by='cluster']
      
      res1[,tops.cluster:=top.central|top.NES|top.pval,by='cluster']
      return(res1)
      
    }else{
      return(res1)
    }
    
  }))
  
  if(return_full)
    return(res)
  else
    return(res[(tops.cluster)])
  
}