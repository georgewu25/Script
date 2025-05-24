plot_WGCNA_threshold <- function(df, outdir, metadata) {
  
  #Remove rows with missing protein IDs
  df <- df %>%
    dplyr::filter(!is.na(protein_Id))
  
  df_wide <- df %>%
    select(protein_Id, Sample, transformedIntensity) %>%
    pivot_wider(names_from = Sample, values_from = transformedIntensity)
  
  #Filter proteins by presence (>0.1 intensity in 95% of samples)
  keep <- rowSums(df_wide[, -1] > quantile(df$transformedIntensity, probs = 0.1, na.rm = TRUE)
, na.rm = TRUE) >= ceiling(ncol(df_wide) * 0.95)
  df_filtered <- as.data.frame(df_wide[keep, ])
  
  rownames(df_filtered) <- df_filtered$protein_Id
  df_filtered <- df_filtered[, -1]

  # Filter data based on variance, choosing the top 90% most variable genes:
  df_filtered$variance <- apply(df_filtered[,- 1], 1, var) # Create a column of calculated variances.
  df_var <- df_filtered[df_filtered$variance >= quantile(df_filtered$variance, 0.10, na.rm = TRUE), ] # Use the 10% quantile as the cut-off.

  df_var$variance <- NULL # Remove the variance column.
  
  #Match samples
  rownames(metadata) <- metadata$abundance_id_final
  meta <- metadata[colnames(df_var), , drop = FALSE]
  
  # Adjust for batch effect using the "comBat" function:
  modcombat <- model.matrix(~1, data = meta)

  df_combat <- ComBat(dat=data.matrix(df_var), batch = metadata$Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  #Covariate Correction
  #Convert categorical variables
  meta$Sex  <- factor(meta$Sex, levels = c("female","male"))
  
  #covs <- model.matrix(~ PMI..hr. + Race + Sex, data = meta)
  covs <- model.matrix(~ PMI..hr. + Race + Sex + CERAD + Braak + Diagnosis + Age, data = meta)
  covs <- covs[, colnames(covs) != "(Intercept)", drop=FALSE]
  
  df_final <- removeBatchEffect(df_combat, covariates = covs)
  
  # Transpose corrected data so the rows are samples and columns are genes
  df_final_t <- as.data.frame(t(df_final))

  #Save data frame
  if (!dir.exists(outdir)) {
    dir.create(outdir, mode = "0755", recursive = T)
  }
  write.csv(df_final_t, paste0(outdir, "wgcna_matrix.csv"), row.names = T, quote = F)
  
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



wgcna_tom <- function(df, outdir, sft, deepsplit_n) {
  
  TOM <-  blockwiseModules(df[, -1], maxBlockSize=100000, minBlockSize=0, minModuleSize=30,
                            corType="bicor", maxPOutliers=0.10, pearsonFallback="individual",
                            power=sft, networkType="signed", TOMType="signed", reassignThreshold=1E-8,
                            mergeCutHeight=0.3, deepsplit=deepsplit_n, numericLabels=TRUE, verbose=4)
   
   outdir <- paste0(outdir, "/deepsplit", deepsplit_n)
   if(!file.exists(outdir)) {
      dir.create(outdir, mode="0755", recursive=TRUE)
  }
   
   saveRDS(TOM, paste0(outdir, "/TOM_deepsplit", deepsplit_n, ".rds"))
}


plot_wgcna <- function(tom, df, outdir) {
  
  # Label modules using colors, grey color is reserved for unassigned genes.
  module_color <- labels2colors(tom$colors, zeroIsGrey=TRUE)
  names(module_color)<-names(tom$colors)
  
  # Calculate module eigengenes, which are the first PCs.
  ME <- tom$MEs
  # Assign colors to MEs
  ME_color <- moduleEigengenes(df[, -1], module_color)$eigengenes
  # Order MEs by hierarchical clustering
  ME <- orderMEs(ME_color)
  
  # Plot the cluster tree of module eigengenes
  plotEigengeneNetworks(ME, "Astrocytes, with KOs \nModule eigengene correlation dendrogram", plotDendrograms=TRUE, plotHeatmaps=FALSE, plotAdjacency=FALSE, excludeGrey=TRUE, marDendro=c(1,4,4,1))
  
  # Merge modules based on enrichment validation
  module_merged <- mergeCloseModules(df[, -1], module_color, cutHeight=0.6, verbose=4)
  
  # Combine old module colors with merged module colors for plotting
  color_combined <- cbind(module_color, module_merged$colors)
  
  # Plot dendrogram
  pdf(paste0(outdir, "Module_Dendrogram.pdf"), width = 10, height = 10)
  plotDendroAndColors(tom$dendrograms[[1]], main=paste("Cluster Module Dendrogram"), colors=module_color, groupLabels="", hang=0.03, dendroLabels=FALSE, addGuide=TRUE, guideHang=0.05)
  dev.off()
  
  # Plot dendrogram with both original module color and merged colors
  pdf(paste0(outdir, "Module_Dendrogram_combined_colors.pdf"), width = 10, height = 10)
  plotDendroAndColors(tom$dendrograms[[1]], main=paste("Cluster Module Dendrogram with original and merged colors"), colors=color_combined, groupLabels=c("Module colors", "Merged colors"), hang=0.03, dendroLabels=FALSE, addGuide=TRUE, guideHang=0.05)
  dev.off()
  
  # Plot dendrogram with only the merged module color
  pdf(paste0(outdir, "Module_Dendrogram_merged_colors.pdf"), width = 10, height = 10)
  plotDendroAndColors(tom$dendrograms[[1]], main=paste("Cluster Module Dendrogram with merged colors"), colors=module_merged$colors, groupLabels=c("Merged colors"), hang=0.03, dendroLabels=FALSE, addGuide=TRUE, guideHang=0.05)
  dev.off()
  
  # Plot the eigengene correlation from original module
  pdf(paste0(outdir, "Module_Eigengene_Correlation_Dendrogram.pdf"), width = 10, height = 10)
  plotEigengeneNetworks(ME, "Module Eigengene Correlation Dendrogram", plotDendrograms=TRUE, plotHeatmaps=FALSE, plotAdjacency=FALSE, excludeGrey=TRUE, marDendro=c(1,4,4,1))
  dev.off()
  
  # Plot the eigengene correlation from merged module
  pdf(paste0(outdir, "Merged_Module_Eigengene_Correlation_Dendrogram.pdf"), width = 10, height = 10)
  plotEigengeneNetworks(module_merged$newMEs, "Merged Module Eigengene Correlation Dendrogram", plotDendrograms=TRUE, plotHeatmaps=FALSE, plotAdjacency=FALSE, excludeGrey=TRUE, marDendro=c(1,4,4,1))
  dev.off()
  
  # Plot the eigengene heatmap
  pdf(paste0(outdir, "Module_Eigengene_Correlation_Heat map.pdf"), width = 10, height = 10)
  plotEigengeneNetworks(ME, "Module Eigengene Correlation Heatmap", plotHeatmaps=TRUE, plotDendrograms=FALSE, plotAdjacency=FALSE, excludeGrey=FALSE, marHeatmap=c(8,10,3,3), xSymbols=names(ME), ySymbols=names(ME), printAdjacency = TRUE)
  dev.off()
  
  # Plot the eigengene heatmap for merged module
  pdf(paste0(outdir, "Merged_Module_Eigengene_Correlation_Heatmap.pdf"), width = 10, height = 10)
  plotEigengeneNetworks(module_merged$newMEs, "Merged Module Eigengene Correlation Heatmap", plotHeatmaps=TRUE, plotDendrograms=FALSE, plotAdjacency=FALSE, excludeGrey=FALSE, marHeatmap=c(8,10,3,3), xSymbols=names(module_merged$newMEs), ySymbols=names(module_merged$newMEs), printAdjacency = TRUE)
  dev.off()
}



process_tom <- function(file_main_dir, file_sub_dir, filename) {
  
  tom <- readRDS(paste0(file_main_dir, file_sub_dir, filename))
  
  wgcna_mtx <- read.csv(paste0(file_main_dir, "/wgcna_matrix.csv"))
  
  module_color <- labels2colors(tom$colors, zeroIsGrey=TRUE)
  names(module_color)<-names(tom$colors)
  
  # Calculate module eigengenes, which represent the expression profile of their respective module. They are the first PCs.
  ME <- tom$MEs
  # Assign colors to MEs
  ME_color <- moduleEigengenes(wgcna_mtx[, -1], module_color)$eigengenes
  rownames(ME_color) <- wgcna_mtx[,1]
  # Order MEs by hierarchical clustering
  ME <- orderMEs(ME_color)

  df_filtered_t <- as.data.frame(t(wgcna_mtx[, -1]))
  colnames(df_filtered_t) <- wgcna_mtx[, 1]
  df_filtered_t$gene_name <- rownames(df_filtered_t)
  
  # Save gene id and module color infromation from merged module into a table
  modules_df<-data.table(gene_id=names(module_color),module=module_color)
  
  # Merge the ME group with the expression data on gene_id
  modules_df<-merge(modules_df, df_filtered_t, by.x = 'gene_id', by.y = "gene_name")
  
  # Order by ME group
  modules_df<-modules_df[order(module)]
 
  write.csv(modules_df, paste0(file_main_dir, file_sub_dir, "/module_genes_df.csv"), row.names = FALSE, quote = FALSE)
}



run_enrichment <- function(module_df, file_main_dir, outdir) {
  
  res_enr<-rbindlist(lapply(split(pathwaysf, by='subcat'),function(msigdbf)OR3(split(module_df$gene_id,module_df$module),
                                                                            terms_list = split(msigdbf$gene,msigdbf$pathway),
                                                                            background =module_df$gene_id)))

  #add subcategory and pathway size info
  res_enr_final<-merge(res_enr[padj<0.1&n.overlap>5],unique(pathways_infos,by='term'),by='term')[order(query,term,pval)]


  write.csv(res_enr_final,paste0(file_main_dir, outdir, '/functional_enrichment_padj0.1_overlap5.csv'), row.names = FALSE, quote = FALSE)
  
}



plot_trait_heatmap <- function(eigen_gene, trait_data, expression_data, filename) {

  # Correlate modules with traits using Spearman's correlation function:
  module_trait_cor <- cor(eigen_gene, trait_data, use="p", method="spearman") # Obtain correlation value
  
  module_trait_p <- corPvalueStudent(module_trait_cor, nrow(expression_data)) # Obtain Student asymptotic p-values. 
  
  text_p <- paste(signif(module_trait_cor, 2), "\n(", round(signif(module_trait_p, 2), 3), ")", sep="")
  text_value <- signif(module_trait_cor, 2)
  dim(text_value) <- dim(module_trait_cor)
  
  outdir <- paste0(wgcna_outdir, '/Total_Proteome/deepsplit2')
  
  # Create figure with * instead of p-values and no correlation values:
  # p<0.001 '***', p<0.01 '**', p<0.05 '*', p<0.1 '.' 
  module_trait_p_mark <- module_trait_p
  module_trait_p_mark[module_trait_p_mark < 0.001] <- "***"
  module_trait_p_mark[module_trait_p_mark < 0.01 & module_trait_p_mark >= 0.001] <- "**"
  module_trait_p_mark[module_trait_p_mark < 0.1 & module_trait_p_mark >= 0.01] <- "*"
  module_trait_p_mark[module_trait_p_mark > 0.1] <- ""
  
  module_trait_p_mark <- sapply(as.data.frame(module_trait_p), function(col) {
    sapply(col, function(x) {
      if (x < 0.001) {
        return("\u2731\u2731\u2731")
      } else if (x < 0.01 & x >= 0.001) {
        return("\u2731\u2731")
      } else if (x < 0.1 & x >= 0.01) {
        return("\u2731")
      } else {
        return(" ")
      }
    })
  })
  
  color=colorRampPalette(c("dodgerblue", "white", "firebrick1"))(50)
  
  pheatmap::pheatmap(module_trait_cor, cluster_cols = FALSE, cluster_rows = T, color =  color, angle_col = "45", main = filename, fontsize_row = 7, fontsize_col = 10, na_col = "white", cellwidth = 25, display_numbers = as.matrix(module_trait_p_mark), number_color = "black", breaks = seq(-1, 1, length.out = 51))
}

#Numberical variables
num_vars<-colnames(trait_df_final)[c(2:5)]

#Categorical variables
cat_vars <- setdiff(colnames(trait_df_final), num_vars)

# Convert categorical variables to dummy variables
for (i in cat_vars) {
  trait_df_final[, colnames(trait_df_final) == i] <- as.numeric(factor(trait_df_final[, colnames(trait_df_final) == i])) -1
}

rownames(trait_df_final) <- metadata_final$abundance_id_final

# Match sample names
rows_matched <- match(total_p_wgcna_mtx[, 1], rownames(trait_df_final))
trait_data_final <- trait_df_final[rows_matched, ]

module_emma <- function(module_pathway, subdir) {
  
  #Change colnames to meet Emmaplot parameters
  colnames(module_pathway)[1] <- "pathway"
  colnames(module_pathway)[13] <- "NES"
  
  module_pathway$split_genes <- vector("list", nrow(module_pathway))
  
  # Loop through each row and split the string, storing the result in the new column
  for (i in 1:nrow(module_pathway)) {
    genes <- module_pathway$genes.overlap
    module_pathway$split_genes[[i]] <- strsplit(genes, "\\|")[[1]]
  }
  
  colnames(module_pathway)[19] <- "leadingEdge"
  
  
  for (i in 1:length(unique(module_pathway$query))) {
    df <- module_pathway[module_pathway$query. == unique(module_pathway$query)[i], ]
    
    df_filtered <- df[df$padj <= 0.05, ]
    
    df_sorted <- df_filtered[order(abs(df_filtered$NES), decreasing = T), ]
    
    df_sorted <- df_sorted[!duplicated(df_sorted$pathway), ]
    
    #Top 50 pathways
    n.pathway <- min(nrow(df_sorted), 50)
    
    outdir <- paste0(wgcna_outdir, subdir, "/Modules/")
    if (!dir.exists(outdir)) {
      dir.create(outdir, mode = "0755", recursive = T)
    }
    
    if (n.pathway > 2) {
      
      emmaplot(df_sorted[1:n.pathway, ])
      ggsave(paste0(outdir, unique(module_pathway$query)[i], "_Emma.pdf"), height = 10, width = 20)
    } else{
    print(paste0("No Pathways for ", unique(module_pathway$query)[i]))
    }
  }
}


ggplot(enr_res_plot_df, aes(x=-log10(padj), y=term, fill=query.))+
  geom_col() +
  facet_grid(rows = vars(query.), scales = "free") +
  scale_fill_manual(values = c(red="red")) +
  geom_vline(xintercept = -log10(0.1), linetype = "solid", color = "red") + 
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 30),
    axis.title.x = element_text(size = 20),
    panel.margin = unit(.05, "lines"),
    panel.border = element_blank(),
    strip.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.line.x.top = element_blank()
  ) +
  guides(fill = FALSE) +
  scale_x_continuous(expand = c(0, 0), breaks = c(seq(0, 25, by = 5)), limits = c(0, 6)) +
  labs(title = "", x = "-Log10(FDR)", y = "", fill = "")