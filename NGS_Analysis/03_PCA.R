
#PCA
library(tidyverse)
#BiocManager::install("ggfortify")
library(ggfortify)

#Takes file index in the deseq2_conf_liststr
make_pca <- function(x){
  comparison <- deseq2_conf_list[x]
  dir <- paste0(deseq2_root, comparison)
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "deseq2_normalized_counts.csv"))
  
  #Assign genes as rownames
  rownames(df) <- df[, 1]
  
  matrix_data <- as.matrix(df[, -1])
  
  # Perform the PCA
  pca_result <- prcomp(t(matrix_data))
  
  pc_eigenvalues <- pca_result$sdev^2
  
  pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), variance = pc_eigenvalues) %>% 
    #Percent variance
    mutate(pct = variance/sum(variance)*100) %>% 
    #Cumulative variance explained
    mutate(pct_cum = cumsum(pct))
  
  #Percentage explained by PC Visualization
  prct_pc <- pc_eigenvalues %>% 
    ggplot(aes(x = PC)) +
    geom_col(aes(y = pct)) +
    geom_line(aes(y = pct_cum, group = 1)) + 
    geom_point(aes(y = pct_cum)) +
    labs(x = "Principal component", y = "Percent variance explained", title = paste0(comparison, " Top PCs")) +
    theme_minimal()
  
  output_subdir <- paste0("/output/PCA/", comparison, "/")
  if(!file.exists(output_subdir)) {
    dir.create(output_subdir, mode="0755", recursive=TRUE)
  }
  
  #PC file name
  plot1 <- paste0(output_subdir, comparison, "_", "PC_percentage.pdf")
  
  #Save as pdf
  ggsave(plot1, plot = prct_pc, height = 10, width = 10)
  
  #Shorten the sample name to fit the label
  substrings_to_replace <- c("TCW", "\\.Ast")
  for (i in substrings_to_replace) {
    row.names(pca_result$x) <- gsub(i, "", row.names(pca_result$x))
  }
  
  #Get PC scores
  pc_scores <- as_tibble(pca_result$x, rownames = "Sample")
  
  #Sample Cluster Visualization
  pca_plot <- autoplot(pca_result, data = pc_scores, colour = 'Sample', size = 5, main = paste0(comparison, "_PCA")) + theme_minimal()
  
  #PCA Cluster file name
  plot2 <- paste0(output_subdir, comparison, "_", "PCA_cluster.pdf")
  
  #Save as pdf
  pdf(plot2, height = 10, width = 10)
  print(pca_plot)
  dev.off()
  
  #Get PC dimensions
  pc_loadings <- as_tibble(pca_result$rotation, rownames = "gene")
  
  #Get top 10 genes by PCs
  top_genes <- pc_loadings %>% 
    # Select the PCs of interest
    dplyr::select(gene, PC1, PC2) %>%
    pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
    dplyr::group_by(PC) %>% 
    dplyr::arrange(desc(abs(loading))) %>%
    # take the 10 top rows
    dplyr::slice(1:10) %>% 
    # pull the gene column as a vector
    dplyr::pull(gene) %>% 
    # ensure only unique genes are retained
    unique()
  
  top_loadings <- pc_loadings %>% 
    filter(gene %in% top_genes)
  
  gene_pca <- ggplot(data = top_loadings) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
                 arrow = arrow(length = unit(0.1, "in"))) +
    geom_text(aes(x = PC1, y = PC2, label = gene),
              nudge_y = 0.005, size = 3) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    labs(title = paste0(comparison, " PCA"), x = "PC1", y = "PC2") +
    theme_minimal()
  
  pdf(paste0(output_subdir, comparison, "_PCA_Eigengene.pdf"), height = 10, width = 10)
  print(gene_pca)
  dev.off()
}

for (i in 1:length(deseq2_conf_list)) {
  make_pca(i)
}