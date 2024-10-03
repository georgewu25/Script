#Get all pathways from each gmt for each comparison
pathway_all_df <- data.frame()

for (i in 1:length(deseq2_conf_list)) {
  comparison <- deseq2_conf_list[i]
  gsea_dir <- paste0(gsea_root, comparison, "/")
  
  temp <- readRDS(paste0(gsea_dir, comparison, "_c2.cp.v2023.1.Hs.symbols.gmt_result.rds"))
  temp_df <- data.frame()
  temp_df <- rbind(temp_df, temp)
  temp_df$Comparison <- deseq2_name_list[i]
  
  pathway_all_df <- rbind(pathway_all_df, temp_df)
}

# Takes a vector of indecies of deseq config files.
# x will be the indecies of the datasets, y will be the indecies for endpoint visualization
get_overlap_pathway_time <- function(x, y, colname, filename){
  
  outdir <- "/output/Pathway/overlapped_pathway/"
  
  overlapped_df <- data.frame()
  
  for (i in 1:length(x)) {
    df_subset <- pathway_all_df[pathway_all_df$Comparison %in% deseq2_name_list[x[i]], ]
    
    overlapped_df <- rbind(overlapped_df, df_subset)
  }
  
  #Get pathways that are present in all data sets
  pathway_count <- table(overlapped_df$pathway)
  
  overlapped_pathway <- names(pathway_count)[which(pathway_count == length(x))]
  
  overlapped_pathway_df <- overlapped_df[overlapped_df$pathway %in% overlapped_pathway, ]
  
  #Get significant p values
  p_count <- table(overlapped_pathway_df$pathway[overlapped_pathway_df$padj < 0.05])
  
  overlapped_pathway_df$sig_p <- sapply(overlapped_pathway_df$pathway, function(x) p_count[x])
  
  #Remove pathways that are insignificant in all conditions
  overlapped_pathway_df <- overlapped_pathway_df[!is.na(overlapped_pathway_df$sig_p),]
  
  #Rank the pathways by standard deviations across conditions
  overlapped_pathway_df$sd <- sapply(overlapped_pathway_df$pathway, function(x) {
    sd(overlapped_pathway_df[overlapped_pathway_df$pathway == x, ]$NES)
  })
  
  overlapped_pathway_df_sorted <- overlapped_pathway_df[order(overlapped_pathway_df$sd, decreasing = F), ]
  
  heatmap_df <- overlapped_pathway_df_sorted[, c(1,6,9)]
  
  heatmap_df$Comparison <- factor(heatmap_df$Comparison, levels = c(deseq2_name_list[x]))
  
  pathway_matrix <- as.data.frame(reshape(heatmap_df, idvar = "pathway", timevar = "Comparison", direction = "wide"))
  
  colnames(pathway_matrix)[2:ncol(pathway_matrix)] <- colname
  
  rownames(pathway_matrix) <- pathway_matrix$pathway
  
  #Take the top 10 most and least opposite pathways for heatmap
  top_pathways <- rbind(head(pathway_matrix, 10), tail(pathway_matrix, 10))
  
  
  pval_df <- overlapped_pathway_df_sorted[overlapped_pathway_df_sorted$pathway %in% top_pathways$pathway, c(1,3,9)]
  
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
  
  pheatmap::pheatmap(as.matrix(top_pathways[,-1]), cluster_cols = FALSE, cluster_rows = T, color =  color, angle_col = "0", main = "", fontsize_row = 7, fontsize_col = 10, na_col = "white", cellwidth = 25, display_numbers = as.matrix(pval_matrix[,-1]), number_color = "white")
  
  #Take the top 100 most and least opposite pathways for Emma Plot
  top_pathways2 <- rbind(head(pathway_matrix, 50), tail(pathway_matrix, 50))
  
  top_pathways2_df <- pathway_all_df[pathway_all_df$Comparison == deseq2_name_list[y] & pathway_all_df$pathway %in% top_pathways2$pathway, ]
  
  #emmaplot(top_pathways2_df, label.size = 3)
  #ggsave(paste0(outdir, filename, "_opposite_pathway_emma_plot.pdf"), height = 10, width = 15)
}


get_overlap_pathway <- function(ref, target, pathways, filename){
  
  outdir <- "/output/Pathway/overlapped_pathway/"
  
  #Get the pathways in the reference dataset
  df_ref <- pathway_all_df[pathway_all_df$Comparison == deseq2_name_list[ref], ]
  
  df_ref_sorted <- df_ref[df_ref$pathway %in% pathways, ]
  
  #Get the pathways in the target data set
  df_target <- pathway_all_df[pathway_all_df$Comparison %in% deseq2_name_list[target] & pathway_all_df$pathway %in% df_ref_sorted$pathway, ]
  
  #Save the overlapped pathways
  overlapped_pathway <- rbind(df_target, df_ref_sorted)
  
  overlapped_pathway$Comparison <- factor(overlapped_pathway$Comparison, levels = c(deseq2_name_list[target], deseq2_name_list[ref]))
  
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
  #colnames(pathway_matrix)[c(2,3)] <- c("APOE 44 vs APOE 33 D0h", "APOE 44 vs APOE 33 D24h")
  pheatmap::pheatmap(as.matrix(pathway_matrix[,-1]), cluster_cols = FALSE, cluster_rows = T, color =  color, angle_col = "45", main = "", fontsize_row = 7, fontsize_col = 10, na_col = "white", cellwidth = 25, display_numbers = as.matrix(pval_matrix[,-1]), number_color = "white")
}