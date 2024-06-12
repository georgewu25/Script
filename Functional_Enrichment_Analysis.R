
#Gene-module
module_gene_df <- read.csv("/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/module_genes.csv")

astro_subtype <- read.csv("/RefData/Astrocytes_signatures/astrocyte_signatures_trans_human.csv")

colnames(astro_subtype)

astro_subtype <- astro_subtype[, c(3,7)]
astro_subtype <- astro_subtype[astro_subtype$Symbol.human != "", ]

#Gene heatmap for each subtype
astro_subtype[duplicated(astro_subtype$Symbol.human) | duplicated(astro_subtype$Symbol.human, fromLast = TRUE) ,]

for (i in unique(astro_subtype$type)) {
  print(nrow(astro_subtype[astro_subtype$type == i, ]))
}

subtype_gene_module_df <- merge(module_gene_df, astro_subtype, by.x = "gene_symbol", by.y = "Symbol.human", all = FALSE)

library(data.table)


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




pathways<-fread('/MSigDB/all_CPandGOs_gene_and_genesets.csv.gz')

pathways_infos<-fread('/MSigDB/all_CPandGOs_genesets_metadata.csv.gz')
setnames(pathways_infos,old = 'pathway','term')

#test pathways <2k et >5k
pathwaysf<-pathways[pathway.size>5&pathway.size<2000]
length(unique(pathwaysf$pathway))#12k

#Gene-module
module_gene_df <- read.csv("/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/module_genes.csv")

#rm non annotated or unassigned genes
module_gene_df<-module_gene_df[!(is.na(module_gene_df$gene_symbol)|module_gene_df$gene_symbol=='')&module_gene_df$module!='grey', ] #grey is for unassigned genes

res_enr<-rbindlist(lapply(split(pathwaysf, by='subcat'),function(msigdbf)OR3(split(module_gene_df$gene_symbol,module_gene_df$module),
                                                                             terms_list = split(msigdbf$gene,msigdbf$pathway),
                                                                             background =module_gene_df$gene_symbol)))


#add subcategory and pathway size info
res_enr<-merge(res_enr[padj<0.1&n.overlap>5],unique(pathways_infos,by='term'),by='term')[order(query,term,pval)]

#write.csv(res_enr,'/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/functional_enrichment_padj0.1_overlap5.csv', row.names = FALSE, quote = FALSE)


#Functional Enrichment for only astrocyte subtype signatures
ast_enr <- OR3(split(module_gene_df$gene_symbol,module_gene_df$module), terms_list = split(astro_subtype$Symbol.human,astro_subtype$type), background =module_gene_df$gene_symbol)

#write.csv(ast_enr,'/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/astrocyte_subtype_functional_enrichment_padj0.1_overlap5.csv', row.names = FALSE, quote = FALSE)

library(dplyr)
library(ggplot2)

ast_enr <- read.csv('/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/astrocyte_subtype_functional_enrichment_padj0.1_overlap5.csv')

color_scale <- unique(ast_enr$query.)

ggplot(ast_enr, aes(x = term, y = n.overlap, fill = query.)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Module Overlapped Genes Stack Bar Plot", x = "Astrocyte Subtype", y = "Number of overlapped genes") +
  scale_fill_manual(values = color_scale, name = "Modules") +
  theme_minimal()

# Enrichment heatmap for all subtypes
library(reshape2)

# Melt the data
df_melt <- melt(ast_enr, id.vars = c("query.", "term"), measure.vars = c("fold.enrichment", "padj"))

# Spread the data back into a wide format with fold.enrichment and padj as separate columns
df_wide <- dcast(df_melt, query. + term ~ variable)

# Create a new column for asterisks based on padj values
df_wide <- df_wide %>%
  dplyr::mutate(significance = ifelse(padj < 0.1, "*", ""))

# Plot the heatmap
ggplot(df_wide, aes(x = term, y = query.)) +
  geom_tile(aes(fill = fold.enrichment), color = "white") +
  scale_fill_gradient(low = "blue", high = "red", name = "Fold_Enrichment") +
  geom_text(aes(label = round(fold.enrichment, 2)), color = "black") +
  geom_text(aes(label = significance), color = "white", vjust = 2, size = 5) +
  theme_minimal() +
  labs(title = "Module Overlapped Pathway Fold Enrichment Heatmap", x = "Astrocyte Subtype", y = "Module")