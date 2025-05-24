get_deg <- function(deseq_dir, index, pval) {
  
  deg_df <- data.table(
    Comparison = character(length(index)),
    DEG_count = integer(length(index)),
    Up = integer(length(index)),
    Down = integer(length(index))
  )
  
  for (i in 1:length(index)) {
    
    file <- index[i]
    comparison <- basename(list.dirs(deseq_dir, recursive = F))[file]
    
    df <- read.csv(paste0(deseq_dir, "/", comparison, "/Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
    
    #p_value threshold
    df_subset <- df[!is.na(df$padj) & df$padj <= pval, ]
    
    deg_df[i, Comparison := comparison]
    deg_df[i, DEG_count := nrow(df_subset)]
    deg_df[i, Up := sum(df_subset$log2FoldChange > 0)]
    deg_df[i, Down := sum(df_subset$log2FoldChange < 0)]
  }
  
  #Melt data
  plot_df <- tidyr::pivot_longer(deg_df, cols = c(Up, Down), names_to = "DEG", values_to = "value")
  plot_df$DEG <- factor(plot_df$DEG, levels = c("Up", "Down"))
  
  return(plot_df)
}


plot_deg <- function(df, bar_width, x_spacing, y_max, title) {
  
  ggplot(df, aes(x = Comparison, y = value, fill = DEG)) +
    geom_bar(stat = "identity", position = "stack", width = bar_width) +
    facet_grid(~ Comparison) +
    scale_fill_manual(values = c("red", "blue"), name = "") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          axis.line = element_line(color = "black", linewidth = 0.2),
          #panel.margin = unit(0, "lines"),
          strip.text = element_text(size = 14, face = "bold"),
          strip.background = element_rect(fill = "lightblue", color = "black"),
    ) +
    scale_x_discrete(expand = c(x_spacing, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
    labs(title = title, x = "", y = "DEG Count")
}


plot_volcano <- function(deseq_dir, index, fc, padj) {
  
  comparison <- basename(list.dirs(deseq_dir, recursive = F))[index]
  
  df <- read.csv(paste0(deseq_dir, "/", comparison, "/Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
  
  df$expression <- ifelse(abs(df$log2FoldChange) > fc & df$padj < padj, "FDR & Log2FC", 
                          ifelse(abs(df$log2FoldChange) < fc & df$padj < padj, "FDR",
                                 ifelse(abs(df$log2FoldChange) > fc, "Log2FC", "Not Sig")))
  
  #Remove NAs
  df <- df[!is.na(df$expression),]
  df$expression <- factor(df$expression, levels = c("FDR & Log2FC", "FDR", "Log2FC", "Not Sig"))
  
  top_gene_df <- df[df$expression == "FDR & Log2FC", ]
  top_gene_df <- top_gene_df[order(abs(top_gene_df$log2FoldChange), decreasing = T), ]
  
  x_min <- ifelse(nrow(top_gene_df) > 2, min(min(top_gene_df$log2FoldChange) * 1.05, -5), -5)
  x_max <- ifelse(nrow(top_gene_df) > 2, max(max(top_gene_df$log2FoldChange) * 1.05, 5), 5)
  
  y_min <- ifelse(nrow(top_gene_df) > 2, 0, -log10(median(df$padj, na.rm = TRUE)))
  y_max <- ifelse(nrow(top_gene_df) > 2, -log10(min(df$padj, na.rm = TRUE)) * 1.2, -log10(min(df$padj, na.rm = TRUE)) * 1.2)
  
  ggplot(df, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(color = expression),size = 2/5) +
    geom_hline(yintercept = -log10(padj), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-fc, fc), linetype = "dashed", color = "black") +
    theme_minimal()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", color = "black")) +
    scale_color_manual(
      values = c("FDR & Log2FC" = "red", 
                 "FDR" = "green", 
                 "Log2FC" = "blue", 
                 "Not Sig" = "grey50")
    )  +
    guides(colour = guide_legend(override.aes = list(size=2.5))) +
    geom_label_repel(data = top_gene_df[1:5, ],
                     mapping = aes(log2FoldChange, -log(padj,10), label = gene_name),
                     size = 3) +
    scale_x_continuous(expand = c(0, 0), limits = c(x_min, x_max)) +
    scale_y_continuous(expand = c(0, 0), limits = c(y_min, y_max)) +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.line = element_line(linewidth = 0.5)) +
    labs(color= paste0("Log2FC Cutoff: ", fc, "\nFDR Cutoff: ", padj, "\n\n"), size = 13) +
    labs(title = comparison,
         x = expression("Log"[2]*"FC"),
         y = expression("-Log"[10]*"FDR"))
}



run_gsea <- function(deseq_dir, index, pathways){
  
  comparison <- basename(list.dirs(deseq_dir, recursive = F))[index]
  
  df <- read.csv(paste0(deseq_dir, "/", comparison, "/Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
  
  #Rank data frame by stat
  df_sorted <- df[order(df$stat),]
  ranks <- setNames(df_sorted$stat, df_sorted$gene_name)
    
  #GSEA
  fgseaRes <- fgsea(pathways, ranks, scoreType='std',nPermSimple = 10000)

}



pathway_bar_plot <- function(deseq_dir, gsea_dir, index, gmt, p_val, title, filename) {
  
  comparison <- basename(list.dirs(deseq_dir, recursive = F))[index]
  
  df <- readRDS(paste0(gsea_dir, "/", comparison, "_", gmt, "_gsea.rds"))
  
  df <- df %>%
    dplyr::filter(padj <= p_val)
  
  if (nrow(df) == 0) {
    print("No pathway match criteria")
    return(NULL)
  } else{
    
    up <- df %>%
      dplyr::filter(NES >= 0 & padj <= 0.05) %>%
      arrange(-abs(NES))
    
    if (nrow(up) == 0) {
      print("No up-regulated pathways match criteria")
      top_up <- data.frame()
    } else{
      top_up <- up[1:min(nrow(up), 5), ]
    }
    
    down <- df %>%
      dplyr::filter(NES <= 0 & padj <= 0.05) %>%
      arrange(-abs(NES))
    
    if (nrow(down) == 0) {
      print("No down-regulated pathways match criteria")
      top_down <- data.frame() 
    } else{
      top_down <- down[1:min(nrow(down), 5), ]
    }
    
    top_pathways <- do.call(rbind, Filter(Negate(is.null), list(top_up, top_down)))
    
    if (nrow(top_pathways) == 0 ) {
      
      return(NULL)
    } else{
      ggplot(top_pathways, aes(x = NES, y = pathway, fill = -log10(padj))) +
        geom_col() +
        scale_fill_gradientn(colors = c("yellow", "red")) +
        theme_light() +
        theme(
          panel.margin = unit(.05, "lines"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.text.y = element_text(size = 8)
        ) +
        labs(title = title, x = "NES", y = "")
      ggsave(paste0(gsea_dir, "/", filename), height = 7, width = 8)
      
    }
  }
}
