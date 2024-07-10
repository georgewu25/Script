#!/bin/usr/Rscript

deseq2_root <- "/mwu/Project/Astrocyte_Ab/2X100/results/DESEQ2/"
gsea_root <- "/mwu/Project/Astrocyte_Ab/2X100/results/GSEA/"
if(!file.exists(gsea_root)) {
      dir.create(gsea_root, mode="0755", recursive=TRUE)
    }

deseq2_conf_liststr <- "Degrade_33_8vs24hr,Degrade_33_8vs48hr,Degrade_33vs44_24hr,Degrade_33vs44_48hr,Degrade_44_8vs24hr,Degrade_44_8vs48hr,Uptake_33_AbvsCtrl,Uptake_33vs44_Ab,Uptake_33vs44_Ctrl,Uptake_44_AbvsCtrl"

deseq2_conf_list <- unlist(strsplit(deseq2_conf_liststr, ","))


#Volcano Plot
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

#Takes data frame and plot title
volcano <- function(df, title){
  p <- EnhancedVolcano(df,
  lab = df$gene_name,
  x = 'log2FoldChange',
  y = 'pvalue',
  title = title,
  pointSize = c(ifelse(df$log2FoldChange>3, 3, 1)),
  pCutoff = 0.1 ,
  FCcutoff = 0.5,
  cutoffLineType = 'twodash',
  cutoffLineCol = 'black',
  cutoffLineWidth = 0.8,
  hline = c(0.1, 1e-10),
  hlineCol = 'black',
  hlineType = 'longdash',
  hlineWidth = 1,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  legendLabels=c('Not sig.','Log (base 2) FC','p-value', 'p-value & Log (base 2) FC'),
  legendPosition = "right",
  legendLabSize = 12,
  legendIconSize = 5.0,
  caption = bquote(~Log[2]~ "fold change cutoff, 0.5; p-value cutoff, 0.1"))
}

for (i in 1:10) {
  comparison <- deseq2_conf_list[i]
  dir <- paste0(deseq2_root, comparison)
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "Results_GTFAnnotated_NoGeneIDDuplicates.csv"))

  filename <- paste0("/mwu/Project/Astrocyte_Ab/2X100/results/output/volcano_plot", comparison, "_Volcano Plot.pdf")
  pdf(filename, height = 10, width = 10)
  print(volcano(df, comparison))
  dev.off()
}



#PCA
library(tidyverse)
#BiocManager::install("ggfortify")
library(ggfortify)
library(biomaRt)

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
    labs(x = "Principal component", y = "Percent variance explained", title = paste0(comparison, " Top PCs"))
  
  output_subdir <- paste0("/mwu/Project/Astrocyte_Ab/2X100/results/output/PCA/", comparison, "/")
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
  pca_plot <- autoplot(pca_result, data = pc_scores, colour = 'Sample')
  
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
    group_by(PC) %>% 
    arrange(desc(abs(loading))) %>%
    # take the 10 top rows
    slice(1:10) %>% 
    # pull the gene column as a vector
    dplyr::pull(gene) %>% 
    # ensure only unique genes are retained
    unique()
  
  #Convert Ensemble_gene_id to gene symbol
  ID_convert <- function(x){
    
    ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
    
    getBM(filters="ensembl_gene_id", attributes="hgnc_symbol", values=x, mart=ensembl)
  }
  
  temp <- 0
  #Remove the version at the end of the Ensemble ID
  for (i in 1:length(top_genes)) {
      temp[i] <- unlist(strsplit(top_genes[i], "\\."))[1]
  }
  
  genes <- ID_convert(temp)
  
  write.table(genes, paste0(output_subdir, comparison, "Top_genes_by_PCs"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")
}

for (i in 1:10) {
  make_pca(i)
}
