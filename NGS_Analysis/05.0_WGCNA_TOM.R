#!/usr/bin/Rscript

#WGCNA
#BiocManager::install("WGCNA")
library("ggplot2") # For PCA plot
library("WGCNA") # For network analysis
library("DESeq2") # For data transformation
#BiocManager::install("sva")
library("sva") # For batch/covariate correction
library("GSA") # To extract gene sets
library("erer") # For list writing
library("car") # For Levene's homogeneity test
library("bestNormalize") # for data transformation
library("RColorBrewer") # For PCA color pallette


plot_WGCNA_threshold <- function() {
  df <- read.table("/results/FeatureCount/featureCounts_clean.txt", header = TRUE)
  
  outdir <- "/output/WGCNA/"
  
  #Filter genes based on expression levels.
  df_filtered <- df[rowSums(edgeR::cpm(df[,-1])>1)>= (ncol(df)-1)*0.95,]
  
  # Filter genes using WGCNA function collapseRows, which finds rows with duplicate gene names and collapses them using two methods. For gene names duplicated exactly 2 times, the gene id with the maximum mean counts across samples is selected. For gene names duplicated 3 or more times, the gene id with the highest connectivity according to a signed weighted correlation network adjacency matrix among the corresponding rows is selected. 
  df_filtered_by_gene <- collapseRows(df_filtered[, -1], df_filtered[, 1], rownames(df_filtered), method="MaxMean", connectivityBasedCollapsing=TRUE, connectivityPower=1)
  
  df_filtered_by_gene <- df_filtered[df_filtered_by_gene$selectedRow, ]
  
  # Normalize the filtered expression data using a variance stabilizing transformation, a DESq2 function:
  df_matrix <- data.matrix(df_filtered_by_gene[, -1])
  df_matrix <- ceiling(df_matrix)
  
  df_matrix <- varianceStabilizingTransformation(df_matrix)
  
  
  df_final <- as.data.frame(df_matrix) # Turn matrix back into data frame.
  df_final <- cbind(df_filtered_by_gene[, 1], df_final) # Add the gene names back to the matrix
  colnames(df_final)[1] <- "gene_name"
  
  
  # Filter data based on variance, choosing the top 90% most variable genes:
  df_final$variance <- apply(df_final[,- 1], 1, var) # Create a column of calculated variances.
  df_final_var <- df_final[df_final$variance >= quantile(df_final$variance, c(.10)), ] # Use the 10% quantile as the cut-off.
  df_final_var$variance <- NULL # Remove the variance column.
  
  # Adjust for covariates using the "comBat" function:
  
  meta <- read.csv("/Metadata.csv")
  
  modcombat <- model.matrix(~1 +RIN, data=meta)
  df_matrix_combat <- data.matrix(df_final_var[,-1])
  
  #Corrected for batch effect
  df_combat <- ComBat(dat=df_matrix_combat, batch = meta$Individual, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  # Transpose corrected data so the rows are samples and columns are genes
  df_final_t <- as.data.frame(t(df_combat))
  names(df_final_t) <- df_final_var[, 1]
  
  #Save data frame
  write.csv(df_final_t, paste0(outdir, "normalized_expression.csv"))
  
  h_tree <- hclust(dist(df_final_t), method="average")
  
  pca_res <- prcomp(df_final_t, scale.=FALSE, center=FALSE)# Perform PCA
  pca_df <- as.data.frame(pca_res$x) # Extract results
  pca_prct <- round(pca_res$sdev / sum(pca_res$sdev) * 100, 2) # Calculate PC %s
  pca_prct <- paste0(colnames(pca_df), " (", paste(as.character(pca_prct), "%)")) # Create PC axis labels
  
  pca_plot <- ggplot(data=pca_df, aes(x=PC1, y=PC2, color=rownames(df_final_t))) + 
    geom_point(size=2) + xlab(pca_prct[1]) + ylab(pca_prct[2]) + 
    theme(legend.title=element_blank()) + ggtitle("PCA")
  
  #Save as pdf
  pdf(paste0(outdir,  "pca.pdf"), height = 10, width = 10)
  print(pca_plot)
  dev.off()
  
  # Test for soft threshold
  powers <- c(1:30)
  sft <- pickSoftThreshold(df_final_t, dataIsExpr=TRUE, powerVector=powers, verbose=5, corFnc="bicor", corOptions="use='p'")
  
  plot_df <- sft$fitIndices
  
  sft_r2 <- ggplot(plot_df, aes(x = Power, y = -sign(slope)*SFT.R.sq)) +
    geom_text(aes(label = powers), color = "red", size = 3) +
    geom_hline(yintercept = 0.90, color = "red") +
    labs(x = "Soft threshold (power)", y = "Scale-free topology fit (R^2)", 
         title = paste0("Scale independence"))
  
  sft_k <- ggplot(plot_df, aes(x = Power, y = mean.k.)) +
    geom_text(aes(label = powers), color = "red", size = 3) +
    labs(x = "Soft threshold (power)", y = "Mean connectivity", 
         title = paste0("Mean Connectivity"))
  
  #Save as pdf
  pdf(paste0(outdir, "soft_scale.pdf"), height = 10, width = 10)
  print(sft_r2)
  dev.off()
  
  #Save as pdf
  pdf(paste0(outdir, "mean_connectivity.pdf"), height = 10, width = 10)
  print(sft_k)
  dev.off()
}

plot_WGCNA_threshold()

#PCA after normalization
test <- read.csv(paste0(outdir, "normalized_expression.csv"))

h_tree <- hclust(dist(test), method="average")

pca_res <- prcomp(test[, -1], scale.=T, center=T)# Perform PCA

pca_df <- as.data.frame(pca_res$x) # Extract results
pca_prct <- round(pca_res$sdev / sum(pca_res$sdev) * 100, 2) # Calculate PC %s
pca_prct <- paste0(colnames(pca_df), " (", paste(as.character(pca_prct), "%)")) # Create PC axis labels

genotype_color <- substr(test[, 1], start = 5, stop = 7)

ggplot(data=pca_df, aes(x=PC1, y=PC2, color=genotype_color)) + 
  geom_point(size=2) + xlab(pca_prct[1]) + ylab(pca_prct[2]) + 
  theme(legend.title=element_blank()) + ggtitle("PCA")


 wgcna_tom <- function() {
   wgcna_dir <- "/output/WGCNA/"
   
   #Read the normalized expression data
   df <- read.csv("/output/WGCNA/normalized_expression.csv")
   
   outdir <- "/output/WGCNA/deepsplit3/"
   
   TOM <-  blockwiseModules(df[, -1], maxBlockSize=100000, minBlockSize=0, minModuleSize=30,
                            corType="bicor", maxPOutliers=0.10, pearsonFallback="individual",
                            power=sft, networkType="signed", TOMType="signed", reassignThreshold=1E-8,
                            mergeCutHeight=0.3, deepsplit=3, numericLabels=TRUE, verbose=4)
   
   saveRDS(TOM, paste0(outdir, "TOM_deepsplit3.rds"))
 }
 
 wgcna_tom()
