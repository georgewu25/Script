#!/bin/usr/Rscript

#DEG Plot, takes index of the deseq2_name_list
plot_deg <- function(x) {
  deg_df <- data.frame(
    Comparison = character(length(x)),
    DEG_count = integer(length(x)),
    Up = integer(length(x)),
    Down = integer(length(x))
  )
  
  for (i in 1:length(x)) {
    
    index <- x[i]
    comparison <- deseq2_conf_list[index]
    name <- deseq2_name_list[index]
    
    dir <- paste0(deseq2_root, comparison)
    subdir <- list.dirs(dir)[2]
    df <- read.csv(paste0(subdir, "/", "Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
    
    #p_value threshold
    df_subset <- df[df$padj <= 0.1, ]
    
    deg_df$Comparison[i] <- name
    deg_df$DEG_count[i] <- nrow(df_subset)
    deg_df$Up[i] <- nrow(df_subset[df_subset$log2FoldChange > 0, ])
    deg_df$Down[i] <- nrow(df_subset[df_subset$log2FoldChange < 0, ])
  }
  
  #Melt data
  plot_df <- tidyr::pivot_longer(deg_df, cols = c(Up, Down), names_to = "DEG", values_to = "value")
  plot_df$DEG <- factor(plot_df$DEG, levels = c("Up", "Down"))
  
  # Change the order of the x-axis groups
  plot_df$Comparison <- factor(plot_df$Comparison, levels = deseq2_name_list[x])
  
  ggplot(plot_df, aes(x = Comparison, y = DEG_count, fill = DEG)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "", x = "", y = "Number of DEGs (FDR<0.1)") +
    scale_fill_manual(values = c("red", "blue"), name = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(color = "black", linewidth = 0.2)) +
    scale_x_discrete(expand = c(0.5, 0)) +
    scale_y_continuous(expand = c(0, 0))
}


#Pathway Plot, takes index of the deseq2_name_list
pathway_bar_plot <- function(x) {
  df <- readRDS(paste0(data_dir, deseq2_conf_list[x], "/", deseq2_conf_list[x], "_c2.cp.v2023.1.Hs.symbols.gmt_result.rds"))
  
  df_sorted <- df[order(df$NES), ]
  top_neg_pathway <- head(df_sorted, 5)
  top_pos_pathway <- tail(df_sorted, 5)
  
  top_pathways <- rbind(top_neg_pathway, top_pos_pathway)
  
  ggplot(top_pathways, aes(x = NES, y = pathway, fill = -log10(padj))) +
    geom_col() +
    scale_fill_gradientn(colors = c("yellow", "red")) +
    theme_light() +
    theme(
      panel.margin = unit(.05, "lines"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.y = element_text(size = 8)
    ) +
    labs(title = deseq2_name_list[x], x = "NES", y = "")
}


#Plot category NES level between genotypes across time
#Takes index of deseq_conf_list and time factor
plot_subtype_time <- function(x, t) {
  
  #Contains GSEA result for the four astrocyte subtypes
  df_all <- data.frame()
  
  for (i in x) {
    
    comparison <- deseq2_conf_list[i]
    
    gsea_res <- readRDS(paste0(data_dir, comparison, "/Astrocyte_subtype_signatures/Astrocyte_subtype_signature_pathway_stats.rds"))
    
    gsea_res$Comparison <- comparison
    
    df_all <- rbind(df_all, gsea_res)
  }
  #Specify the time factor for each comparison
  df_all$time <- rep(t, each = 4)
  
  df_all_subset <- df_all[, c(1,3,6,10)]
  
  ggplot(df_all_subset, aes(x = time, y = NES, color = pathway)) +
    geom_line() +
    geom_point(aes(size = -log10(padj))) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      panel.margin = unit(.05, "lines"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_x_continuous(breaks = unique(df_all_subset$time)) +
    guides(color = guide_legend(order = 1), size = guide_legend(order = 2)) +
    labs(x = "Time (hr)", y = "NES", color = "Category")
  
  ggplot(df_all_subset, aes(x = time, y = NES)) +
    geom_line() +
    geom_point(aes(size = -log10(padj))) +
    scale_size_continuous(name = "-log10(padj)", range = c(4, 10)) +
    facet_grid(~ pathway, space = "free_x") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      panel.margin = unit(.05, "lines"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_x_continuous(breaks = unique(df_all_subset$time)) +
    scale_size_continuous(range = c(2,7)) +
    guides(color = guide_legend(order = 1), size = guide_legend(order = 2)) +
    labs(x = "Time (hr)", y = "NES", color = "Category")
}


#Plot category NES level across time and genotype
#x is the index for the deseq_list, t is the time label
bubble_plot <- function(x, t) {
  #Contains GSEA result for the four astrocyte subtypes
  df_all <- data.frame()
  
  for (i in x) {
    
    comparison <- deseq2_conf_list[i]
    
    gsea_res <- readRDS(paste0("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Pathway/", comparison, "/Astrocyte_subtype_signatures/Astrocyte_subtype_signature_pathway_stats.rds"))
    
    gsea_res$Comparison <- comparison
    
    df_all <- rbind(df_all, gsea_res)
  }
  
  df_all$Genotype <- sapply(df_all$Comparison, function(x) {
    unlist(strsplit(x, "_"))[1]
  })
  
  df_all$time <- rep(rep(t, each = 4),2)
  
  df_all_subset <- df_all[, c(1,3,6,10,11)]
  
  df_all_subset$time <- factor(df_all_subset$time, levels = c("8 hr", "24 hr", "48 hr"))
  
  ggplot(df_all_subset, aes(x = pathway, y = NES, size = -log10(padj), color = Genotype)) +
    geom_point() + 
    scale_size_continuous(name = "-log10(padj)", range = c(4, 10)) +
    facet_grid(rows = vars(time), space = "free_y") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    guides(color = guide_legend(order = 1), size = guide_legend(order = 2)) +
    scale_size_continuous(range = c(3, 8)) +
    labs(x = "Category", y = "NES")
}


# Takes a vector of indecies of deseq config files.
# The first indices will be the reference, where all the pathways that are overlapped will be selected based on that.
get_overlap_pathway <- function(ref, target, filename){
  
  outdir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Pathway/overlapped_pathway/"
  
  #Get the pathways in the reference dataset
  df_ref <- pathway_all_df[pathway_all_df$Comparison == deseq2_name_list[ref], ]
  
  #Take a maximum of 10 top enriched pathways
  df_ref_sorted <- df_ref[order(abs(df_ref$NES), decreasing = T), ]
  df_ref_sorted <- df_ref_sorted[1:min(nrow(df_ref_sorted), 10), ]
  
  #Get the pathways in the target dataset
  df_target <- pathway_all_df[pathway_all_df$Comparison %in% deseq2_name_list[target] & pathway_all_df$pathway %in% df_ref_sorted$pathway, ]
  
  #Save the overlapped pathways
  overlapped_pathway <- rbind(df_ref_sorted, df_target)
  
  overlapped_pathway$Comparison <- factor(overlapped_pathway$Comparison, levels = c(deseq2_name_list[ref], deseq2_name_list[target]))
  
  heatmap_df <- overlapped_pathway[, c(1,6,9)]
  
  pathway_matrix <- as.data.frame(reshape(heatmap_df, idvar = "pathway", timevar = "Comparison", direction = "wide"))
  
  colnames(pathway_matrix)[2:ncol(pathway_matrix)] <- sapply(colnames(pathway_matrix)[2:ncol(pathway_matrix)], function(x) {
    unlist(strsplit(x, "\\."))[2]
  })
  
  rownames(pathway_matrix) <- pathway_matrix$pathway
  
  pval_df <- overlapped_pathway[, c(1,3,9)]
  
  pval_matrix <- as.data.frame(reshape(pval_df, idvar = "pathway", timevar = "Comparison", direction = "wide"))
  
  pval_matrix[, -1] <- lapply(pval_matrix[, -1], function(col) {
    sapply(col, function(x) {
      if (x < 0.1 & x > 0.05) {
        return("\u2731")
      } else if (x < 0.05) {
        return("\u2731\u2731")
      } else {
        return(" ")
      }
    })
  })
  
  color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(50)
  
  pheatmap::pheatmap(as.matrix(pathway_matrix[,-1]), cluster_cols = FALSE, cluster_rows = T, color =  color, angle_col = "45", main = "NES", fontsize_row = 7, fontsize_col = 10, na_col = "white", cellwidth = 50, display_numbers = as.matrix(pval_matrix[,-1]), number_color = "white")
}


get_opposite_pathway <- function(ref, target) {
  
  outdir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Pathway/overlapped_pathway/"
  
  #Get the pathways in the reference dataset
  df_ref <- pathway_all_df[pathway_all_df$Comparison == deseq2_name_list[ref], ]
  #Get the pathways in the target dataset
  df_target <- pathway_all_df[pathway_all_df$Comparison %in% deseq2_name_list[target], ]
  
  df_merged <- merge(df_ref, df_target, by = "pathway")
  
  return(df_merged)
}

#Emma Plot
library(igraph)
library(ggraph)
library(stringr)
library(ggplot2)
library(data.table)

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

plot_opposite_pathway <- function(ref, target, filename) {
  
  outdir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Pathway/overlapped_pathway/"
  
  #Get the pathways in the reference dataset
  df_ref <- pathway_all_df[pathway_all_df$Comparison == deseq2_name_list[ref], ]
  #Get the pathways in the target dataset
  df_target <- pathway_all_df[pathway_all_df$Comparison %in% deseq2_name_list[target], ]
  
  df_merged <- merge(df_ref, df_target, by = "pathway")
  
  #Take a maximum of 80 most oppositely enriched pathways 
  df_merged_sorted <- df_merged[order(abs(df_merged$NES.x - df_merged$NES.y), decreasing = T),]
  df_merged_sorted <- df_merged_sorted[1:min(nrow(df_merged_sorted), 40),]
  
  top_pathways <- df_ref[df_ref$pathway %in% df_merged_sorted$pathway, ]
  
  emmaplot(top_pathways, label.size = 3)
  ggsave(paste0(outdir, filename, "_opposite_pathway_emma_plot.pdf"), height = 10, width = 15)
}

