
#Takes file index in the deseq2_conf_liststr
run_gsea <- function(x){
  comparison <- deseq2_conf_list[x]
  dir <- paste0(deseq2_root, comparison)
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
  
  #Rank data frame by stat
  df_sorted <- df[order(df$stat),]
  ranks <- setNames(df_sorted$stat, df_sorted$gene_name)
  
  for (i in 1:nrow(pathway_df)) {
    pathways <- gmtPathways(paste0(gmt_dir, pathway_df$gmt[i]))
    
    #GSEA
    fgseaRes <- fgsea(pathways, ranks, scoreType='std',nPermSimple = 10000)
    
    
    #Leading edges
    gsea_genes <- data.frame(leadingEdge = sapply(fgseaRes$leadingEdge, paste, collapse = ","))
    rownames(gsea_genes) <- fgseaRes$pathway
    
    #Pathways
    topPathwaysUp <- fgseaRes[NES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown <- fgseaRes[NES < 0][head(order(pval), n=10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    gsea_subdir <- paste0(gsea_root, comparison, "/")
    if(!file.exists(gsea_subdir)) {
      dir.create(gsea_subdir, mode="0755", recursive=TRUE)
    }
    
    #Name each GSEA result by the comparison name
    gsea_result_filename <- paste0(gsea_subdir, comparison, "_", pathway_df$gmt[i], "_result.rds")
    
    #Save the data frame
    saveRDS(fgseaRes, gsea_result_filename)
  }
}

#Define pathways of interest
pathway_df <- data.frame(
  name="C2_CP",
  gmt="c2.cp.v2023.1.Hs.symbols.gmt",
  desc="C2 CP",
  category="C2")


#Two-way plot for uptake experiment
for (i in 1:length(deseq2_conf_list)) {
  run_gsea(i)
}