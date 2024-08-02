#!/bin/usr/Rscript

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
  
  meta <- read.csv("Metadata.csv")
  
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


#PCA after normalization
df <- read.csv(paste0(outdir, "normalized_expression.csv"))

h_tree <- hclust(dist(df), method="average")

pca_res <- prcomp(df[, -1], scale.=T, center=T)# Perform PCA

pca_df <- as.data.frame(pca_res$x) # Extract results
pca_prct <- round(pca_res$sdev / sum(pca_res$sdev) * 100, 2) # Calculate PC %s
pca_prct <- paste0(colnames(pca_df), " (", paste(as.character(pca_prct), "%)")) # Create PC axis labels

genotype_color <- substr(df[, 1], start = 5, stop = 7)

ggplot(data=pca_df, aes(x=PC1, y=PC2, color=genotype_color)) + 
  geom_point(size=2) + xlab(pca_prct[1]) + ylab(pca_prct[2]) + 
  theme(legend.title=element_blank()) + ggtitle("PCA")


wgcna_tom <- function() {
  wgcna_dir <- "output/WGCNA/"
  
  #Read the normalized expression data
  df <- read.csv("output/WGCNA/normalized_expression.csv")
  
  outdir <- "output/WGCNA/deepsplit3/"
  
  TOM <-  blockwiseModules(df[, -1], maxBlockSize=100000, minBlockSize=0, minModuleSize=30,
                           corType="bicor", maxPOutliers=0.10, pearsonFallback="individual",
                           power=sft, networkType="signed", TOMType="signed", reassignThreshold=1E-8,
                           mergeCutHeight=0.3, deepsplit=3, numericLabels=TRUE, verbose=4)
  
  saveRDS(TOM, paste0(outdir, "TOM_deepsplit3.rds"))
}


tom <- readRDS("output/WGCNA/deepsplit2/TOM_deepsplit2.rds")

df <- read.csv("output/WGCNA/normalized_expression.csv")

outdir <- "/output/WGCNA/deepsplit2/"

# Label modules using colors, grey color is reserved for unassigned genes.
module_color <- labels2colors(tom$colors, zeroIsGrey=TRUE)
names(module_color)<-names(tom$colors)

# Calculate module eigengenes, which represent the expression profile of their respective module. They are the first PCs.
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
pdf(paste0(outdir, "Module_Eigengene_Correlation_Heatmap.pdf"), width = 10, height = 10)
plotEigengeneNetworks(ME, "Module Eigengene Correlation Heatmap", plotHeatmaps=TRUE, plotDendrograms=FALSE, plotAdjacency=FALSE, excludeGrey=FALSE, marHeatmap=c(8,10,3,3), xSymbols=names(ME), ySymbols=names(ME), printAdjacency = TRUE)
dev.off()

# Plot the eigengene heatmap for merged module
pdf(paste0(outdir, "Merged_Module_Eigengene_Correlation_Heatmap.pdf"), width = 10, height = 10)
plotEigengeneNetworks(module_merged$newMEs, "Merged Module Eigengene Correlation Heatmap", plotHeatmaps=TRUE, plotDendrograms=FALSE, plotAdjacency=FALSE, excludeGrey=FALSE, marHeatmap=c(8,10,3,3), xSymbols=names(module_merged$newMEs), ySymbols=names(module_merged$newMEs), printAdjacency = TRUE)
dev.off()

#Create module-eigene data frame
#Transform the expression data
df_filtered_t <- as.data.frame(t(df_filtered[, -1]))
colnames(df_filtered_t) <- df_filtered[, 1]
df_filtered_t$gene_name <- rownames(df_filtered_t)

library(data.table)

# Save gene id and module color infromation from merged module into a table
modules_df<-data.table(gene_id=names(module_color),module=module_color)

# Merge the ME group with the expression data on gene_id
modules_df<-merge(modules_df, df_filtered_t, by.x = 'gene_id', by.y = "gene_name")

# Order by ME group
modules_df<-modules_df[order(module)]

#Convert to gene ID to EnsembleID
#Save the gene ids to csv for online conversion using Biotools
write.csv(modules_df$gene_id, paste0(outdir, "genes_to_convert.csv"), row.names = FALSE, quote = FALSE)


gene_symbol <- read.csv(paste0(outdir, "gene_names_converted.csv"))

#Keep the gene ids if conversion to gene symbols failed
empty_rows <- gene_symbol$gene_symbol == ""
gene_symbol$gene_symbol[empty_rows] <- gene_symbol$gene_id[empty_rows]

modules_df$gene_symbol <- gene_symbol$gene_symbol


plot_trait_heatmap <- function(eigen_gene, trait_data, expression_data, filename) {
  
  # Correlate modules with traits using Spearman's correlation function:
  module_trait_cor <- cor(eigen_gene, trait_data, use="p", method="spearman") # Obtain correlation value
  
  module_trait_p <- corPvalueStudent(module_trait_cor, nrow(expression_data)) # Obtain Student asymptotic p-values. 
  
  text_p <- paste(signif(module_trait_cor, 2), "\n(", round(signif(module_trait_p, 2), 3), ")", sep="")
  text_value <- signif(module_trait_cor, 2)
  dim(text_value) <- dim(module_trait_cor)
  
  outdir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/deepsplit2/"
  
  pdf(paste0(outdir, filename, ".pdf"), width = 10, height = 10)
  par(mar=c(9, 10, 3, 3))
  labeledHeatmap(Matrix=module_trait_cor,
                 xLabels=names(trait_data),
                 yLabels=names(eigen_gene),
                 ySymbols=names(eigen_gene),
                 colorLabels=TRUE,
                 colors=blueWhiteRed(30),
                 textMatrix=text_value,
                 setStdMargins=FALSE,
                 cex.text=0.6,
                 zlim=c(-1,1),
                 main=paste0(filename, " Module-trait relationships"))
  dev.off()
  
  # Create figure with * instead of p-values and no correlation values:
  # p<0.001 '***', p<0.01 '**', p<0.05 '*', p<0.1 '.' 
  module_trait_p_mark <- module_trait_p
  module_trait_p_mark[module_trait_p_mark < 0.001] <- "***"
  module_trait_p_mark[module_trait_p_mark < 0.01 & module_trait_p_mark >= 0.001] <- "**"
  module_trait_p_mark[module_trait_p_mark < 0.05 & module_trait_p_mark >= 0.01] <- "*"
  module_trait_p_mark[module_trait_p_mark < 0.1 & module_trait_p_mark >= 0.05] <- "."
  module_trait_p_mark[module_trait_p_mark > 0.1] <- ""
  
  pdf(paste0(outdir, filename, "_p_value.pdf"), width = 10, height = 10)
  par(mar=c(9, 10, 3, 3))
  labeledHeatmap(Matrix=module_trait_cor,
                 xLabels=names(trait_data),
                 yLabels=names(eigen_gene),
                 ySymbols=names(eigen_gene),
                 colorLabels=TRUE,
                 colors=blueWhiteRed(30),
                 textMatrix=module_trait_p_mark,
                 setStdMargins=FALSE,
                 cex.text=1.4,
                 zlim=c(-1,1),
                 main=paste0(filename, " Module-trait relationships with p_value"))
  dev.off()
}

trait_data <- read.csv("Metadata.csv", header=TRUE)

# At 24 hr, it is assumed to have no more uptake
trait_data$Treatment[17:32] <- "N"

colnames(trait_data)[14] <- "Uptake"
colnames(trait_data)[15] <- "Degradation"

# Reformat the sample name
trait_data <- trait_data[, c(1,2,3,5,6,7,10,11,14,15)]

#Numberical variables
num_vars<-colnames(trait_data)[c(2,6,8,10,11)]

#Categorical variables
cat_vars<-colnames(trait_data)[c(3,4,7,9)]

# Convert categorical variables to dummy variables
for (i in cat_vars) {
  trait_data[, colnames(trait_data) == i] <- as.numeric(factor(trait_data[, colnames(trait_data) == i])) -1
}

#Convert individuals to dummy variables due to four levels
dummy <- model.matrix(~ as.factor(trait_data$Individual) - 1)
colnames(dummy) <- levels(as.factor(trait_data$Individual))

# Remove original individual column from the dataframe
trait_data <- trait_data[, -5]

# Combine the original dataframe with the dummy variables dataframe
trait_data_final <- cbind(trait_data, dummy)

# Change sample names to match
trait_data_final$Sample_Name <- gsub("-", ".", trait_data_final$Sample_Name)

# Match the samples
rows_matched <- match(expression_df[, 1], trait_data_final$Sample_Name)

trait_data_final <- trait_data_final[rows_matched, ]

rownames(trait_data_final) <- trait_data_final$Sample_Name

trait_data_final <- trait_data_final[, -1]

#Correlation with all metadata
plot_trait_heatmap(ME, trait_data_final, expression_df, "Trait_analysis")



#Get ME data for important modules
ME_of_interest <- c("MEsalmon", "MEblue", "MEgreenyellow", "MEpurple", "MEturquoise", "MEtan", "MEyellow")
ME_data <- cbind(expression_df[,1], ME[ME_of_interest])
colnames(ME_data)[1] <- "Sample"

#Change sample labels
ME_data[, 1] <- gsub("\\.", "-", ME_data[, 1])

for (i in 1:nrow(eigenegene_comparison_df)) {
  
  df1 <- ME_data[grep(eigenegene_comparison_df$key1[i], ME_data[,1]),]
  
  df1_ordered <- df1[order(factor(substring(df1$Sample, 1, 4), levels = c("TCW1","TCW2","TCW3","TCW4"))),]
  df1_ordered$condition <- eigenegene_comparison_df$condition1[i]
  
  df2 <- ME_data[grep(eigenegene_comparison_df$key2[i], ME_data[,1]),]
  
  df2_ordered <- df2[order(factor(substring(df2$Sample, 1, 4), levels = c("TCW1","TCW2","TCW3","TCW4"))),]
  df2_ordered$condition <- eigenegene_comparison_df$condition2[i]
  
  df_final <- rbind(df1_ordered, df2_ordered)
  
  melted_df <- melt(df_final, id.vars = c("Sample","condition"), variable.name = "Module", value.name = "Expression")
  
  melted_df$Sample <- factor(melted_df$Sample, levels = df_final$Sample)
  
  melted_df$label <- substr(melted_df$Sample, 1, 4)
  
  melted_df$condition <- factor(melted_df$condition, levels = c(eigenegene_comparison_df$condition1[i], eigenegene_comparison_df$condition2[i]))
  
  ggplot(melted_df, aes(x = label, y = Module, fill = Expression)) +
    geom_tile() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(colors = c("red", "blue")) +
    facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size = 12)
    ) +
    xlab(NULL) + ylab(NULL) +
    ggtitle(eigenegene_comparison_df$comparison_name[i])
  ggsave(paste0(outdir, eigenegene_comparison_df$comparison_name[i], "_heatmap.pdf"), height = 10, width = 15)
}

run_ME_stat_test <- function(df, ME_color) {
  df = ME_data_subset
  ME_color = "MEsalmon"
  shapiro.results <- vector(mode="list")
  model <- lm(df[[ME_color]] ~ df$group)
  model_res <- residuals(model)
  shapiro <- shapiro.test(model_res) # Test residuals for normality.
  shapiro.results[[color]] <- shapiro$p.value
  
  levene.results <- vector(mode="list")
  levene <- leveneTest(df[[ME_color]] ~ df$group)
  levene.results[[color]] <- levene$`Pr(>F)`
  
  t.test.results <- vector(mode="list")
  t.test <- t.test(df[[ME_color]] ~ df$group)
  t.test.results[[color]] <- t.test$p.value
}

for (color in ME_of_interest) {
  run_ME_stat_test(ME_data_subset, color)
}



#Gene-module
module_gene_df <- read.csv(paste0(outdir, "module_genes.csv"))

library(data.table)

pathways<-fread('/projectnb/tcwlab/MSigDB/all_CPandGOs_gene_and_genesets.csv.gz')

pathways_infos<-fread('/projectnb/tcwlab/MSigDB/all_CPandGOs_genesets_metadata.csv.gz')
setnames(pathways_infos,old = 'pathway','term')

#test pathways <2k et >5k
pathwaysf<-pathways[pathway.size>5&pathway.size<2000]
length(unique(pathwaysf$pathway))#12k

#rm non annotated or unassigned genes
module_gene_df<-module_gene_df[!(is.na(module_gene_df$gene_symbol)|module_gene_df$gene_symbol==''), ]

res_enr<-rbindlist(lapply(split(pathwaysf, by='subcat'),function(msigdbf)OR3(split(module_gene_df$gene_symbol,module_gene_df$module),
                                                                             terms_list = split(msigdbf$gene,msigdbf$pathway),
                                                                             background =module_gene_df$gene_symbol)))

#add subcategory and pathway size info
res_enr_final<-merge(res_enr[padj<0.1&n.overlap>5],unique(pathways_infos,by='term'),by='term')[order(query,term,pval)]


write.csv(res_enr_final,paste0(outdir, 'functional_enrichment_padj0.1_overlap5.csv'), row.names = FALSE, quote = FALSE)



#Change colnames to meet Emmaplot parameters
colnames(res_enr_copy)[1] <- "pathway"
colnames(res_enr_copy)[13] <- "NES"

res_enr_copy$split_genes <- vector("list", nrow(res_enr_copy))

# Loop through each row and split the string, storing the result in the new column
for (i in 1:nrow(res_enr_copy)) {
  genes <- res_enr_copy[i, 14]
  res_enr_copy$split_genes[[i]] <- strsplit(genes, "\\|")[[1]]
}

colnames(res_enr_copy)[19] <- "leadingEdge"


for (i in 1:length(unique(res_enr_copy$query))) {
  df <- res_enr_copy[res_enr_copy$query. == unique(res_enr_copy$query)[i], ]
  
  df_filtered <- df[df$padj <= 0.05, ]
  
  df_sorted <- df_filtered[order(abs(df_filtered$NES), decreasing = T), ]
  
  df_sorted <- df_sorted[!duplicated(df_sorted$pathway), ]
  
  #Top 50 pathways
  
  n.pathway <- min(nrow(df_sorted), 50)
  
  emmaplot(df_sorted[1:n.pathway, ])
  ggsave(paste0(outdir, unique(res_enr_copy$query)[i], "_Emma.pdf"), height = 10, width = 20)
}

#Module-Pathway Enrichment
module_pathway_list <- list(blue, brown, green, magenta, pink, red, salmon, turquoise, yellow)
colors <- list('blue', 'brown', 'green', 'magenta', 'pink', 'red', 'salmon', 'turquoise', 'yellow')

module_pathway_df <- data.frame()
for (i in 1:length(unique(res_enr_final$query.))) {
  
  df_subset <- res_enr_final[res_enr_final$query. == unique(res_enr_final$query.)[i] & res_enr_final$term %in% unlist(module_pathway_list[i]), ]
  
  module_pathway_df <- rbind(module_pathway_df, df_subset)
}

ggplot(module_pathway_df, aes(x=-log10(padj), y=term, fill=query.))+
  geom_col() +
  facet_grid(rows = vars(query.), scales = "free") +
  scale_fill_manual(values = colors) +
  geom_vline(xintercept = -log10(0.1), linetype = "solid", color = "red") + 
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.margin = unit(.05, "lines"),
    panel.border = element_blank(),
    strip.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2)
  ) +
  guides(fill = FALSE) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, max(-log10(module_pathway_df$padj)), by = 5)) +
  labs(title = "", x = "-Log10(FDR)", y = "", fill = "")



#Takes module specific data
eigene_plot <- function(df, filename) {
  df=module_gene_df[module_gene_df$module == "salmon", ]
  plot_df <- data.frame()
  
  for (i in c("Ast.8hr", "oAb.24hr", "oAb.48hr")) {
    
    expression_df_subset <- df[, grep(i, colnames(df))]
    
    temp <- data.frame()
    
    for (j in c("E44", "E33")) {
      
      expression_df_subset2 <- expression_df_subset[, grep(j, colnames(expression_df_subset))]
      
      df_melted <- expression_df_subset2 %>%
        pivot_longer(cols = everything(), 
                     names_to = "sample", 
                     values_to = "expression")
      
      df_melted$genotype <- j
      df_melted$time <- i
      
      temp <- rbind(temp, df_melted)
      
    }
    
    plot_df <- rbind(plot_df, temp)
    
  }
  
  plot_df$genotype <- ifelse(plot_df$genotype == "E44", "APOE44", "APOE33")
  plot_df$genotype <- factor(plot_df$genotype, levels = c("APOE44", "APOE33"))
  
  plot_df$time <- ifelse(plot_df$time == "Ast.8hr", "8hr", ifelse(plot_df$time == "oAb.24hr", "24hr", "48hr"))
  plot_df$time <- factor(plot_df$time, levels = c("8hr", "24hr", "48hr"))
  
  # Perform the Mann-Whitney U test
  outdir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/deepsplit2/Module_Gene_Analysis/"
  
  ggplot(plot_df, aes(x = time, y = expression, fill = genotype)) +
    geom_boxplot() +
    #stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test", label.y = 1.2, size = 5) +  # Mann-Whitney U test
    labs(x = "Time", y = "Eigengene Expression", fill = "Genotype", title = paste0("Eigengene_in_Module_", filename)) +
    theme_minimal()
  
  ggsave(paste0(outdir, "Module_", filename, "_Eigengene_Boxplot.pdf"), height = 10, width = 15)
}


#Gene Expression level in module salmon
eigene_plot(module_gene_df[module_gene_df$module == "salmon", ], "Salmon")


#Hub Genes
expression_df <- read.csv("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/normalized_expression.csv")

chooseTopHubInEachModule(expression_df[,-1], module_color)

connectivity_allClusters <- intramodularConnectivity.fromExpr(expression_df[,-1], module_color, 
                                                                     corFnc = "bicor", corOptions = "use = 'p'",
                                                                     weights = NULL,
                                                                     distFnc = "dist", distOptions = "method = 'euclidean'",
                                                                     networkType = "signed", power = 9,
                                                                     scaleByMax = FALSE,
                                                                     ignoreColors = if (is.numeric(colors)) 0 else "grey",
                                                                     getWholeNetworkConnectivity = TRUE)

rownames(connectivity_allClusters) <- colnames(expression_df[,-1])
connectivity_allClusters$module_color <- module_color
connectivity_allClusters$gene <- row.names(connectivity_allClusters)
connect <- connectivity_allClusters[order(connectivity_allClusters$module_color,-connectivity_allClusters$kWithin),]

module_yellow <- connect[connect$module_color == "yellow", ]

allClusters_kME <- signedKME(expression_df[,-1],ME,corFnc = "bicor")

ggplot(module_yellow, aes(x = kWithin)) +
  geom_histogram(fill = "lightblue") +
  labs(title = "Gene Intramodular Connectivity Distribution in Epigenetic Module", x = "Connectivity", y = "Gene Count") +
  theme_minimal()


ggplot(common_pathways_df, aes(x = -log10(padj), y = term, fill = query.)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("yellow" = "yellow", "red" = "red"), guide = guide_legend(reverse = TRUE)) +
  geom_vline(xintercept = -log10(0.1), linetype = "solid", color = "red") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2)
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, max(-log10(common_pathways_df$padj)), by = 5)) +
  labs(x = "-log10(FDR)", y = "", title = "Pathway Enrichment Between Two Epigenetic Modules", fill = "Module")


