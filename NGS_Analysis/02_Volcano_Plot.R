#Volcano Plot
library(ggrepel)

plot_volcano <- function(x, fc, padj) {
  
  comparison <- deseq2_conf_list[x]
  dir <- paste0(deseq2_root, comparison)
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
  
  filename <- deseq2_name_list[x]
  
  df$expression <- ifelse(abs(df$log2FoldChange) > fc & df$padj < padj, "FDR & Log2FC", ifelse(abs(df$log2FoldChange) > fc, "Log2FC", "Not Sig"))
  
  df <- df[!is.na(df$expression),]
  
  df$expression <- factor(df$expression, levels = c("FDR & Log2FC", "Log2FC", "Not Sig"))
  
  top_gene_df <- df[df$expression == "FDR & Log2FC", ]
  
  top_gene_df <- top_gene_df[order(abs(top_gene_df$log2FoldChange), decreasing = T), ]
  
  x_min <- ifelse(nrow(top_gene_df) > 2, min(top_gene_df$log2FoldChange) * 1.05, -5)
  x_max <- ifelse(nrow(top_gene_df) > 2, max(top_gene_df$log2FoldChange) * 1.05, 5)
  
  y_min <- ifelse(nrow(top_gene_df) > 2, 0, -log10(median(df$padj, na.rm = TRUE)))
  y_max <- ifelse(nrow(top_gene_df) > 2, -log10(min(df$padj, na.rm = TRUE)) * 1.05, -log10(min(df$padj, na.rm = TRUE)) * 1.05)
  
  ggplot(df, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(color = expression),size = 2/5) +
    geom_hline(yintercept = -log10(padj), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-fc, fc), linetype = "dashed", color = "black") +
    theme_minimal()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", color = "black")) +
    scale_color_manual(values = c("red", "green", "gray50"), drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size=2.5))) +
    geom_label_repel(data = top_gene_df[1:30],
                     mapping = aes(log2FoldChange, -log(padj,10), label = gene_name),
                     size = 4) +
    scale_x_continuous(expand = c(0, 0), limits = c(x_min, x_max)) +
    scale_y_continuous(expand = c(0, 0), limits = c(y_min, y_max)) +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.line = element_line(linewidth = 0.5)) +
    labs(color= paste0("Log2FC Cutoff: ", fc, "\nFDR Cutoff: ", padj, "\n\n"), size = 13) +
    labs(title = filename,
         x = expression("Log"[2]*"FC"),
         y = expression("-Log"[10]*"FDR"))
}