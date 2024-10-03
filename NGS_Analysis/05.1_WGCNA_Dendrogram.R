tom <- readRDS("/output/WGCNA/deepsplit2/TOM_deepsplit2.rds")

df <- read.csv("/output/WGCNA/normalized_expression.csv")

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
pdf(paste0(outdir, "Module_Eigengene_Correlation_Heat map.pdf"), width = 10, height = 10)
plotEigengeneNetworks(ME, "Module Eigengene Correlation Heatmap", plotHeatmaps=TRUE, plotDendrograms=FALSE, plotAdjacency=FALSE, excludeGrey=FALSE, marHeatmap=c(8,10,3,3), xSymbols=names(ME), ySymbols=names(ME), printAdjacency = TRUE)
dev.off()

# Plot the eigengene heatmap for merged module
pdf(paste0(outdir, "Merged_Module_Eigengene_Correlation_Heatmap.pdf"), width = 10, height = 10)
plotEigengeneNetworks(module_merged$newMEs, "Merged Module Eigengene Correlation Heatmap", plotHeatmaps=TRUE, plotDendrograms=FALSE, plotAdjacency=FALSE, excludeGrey=FALSE, marHeatmap=c(8,10,3,3), xSymbols=names(module_merged$newMEs), ySymbols=names(module_merged$newMEs), printAdjacency = TRUE)
dev.off()