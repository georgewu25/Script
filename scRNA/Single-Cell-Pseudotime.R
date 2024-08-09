library(dplyr)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(scales)
library(cowplot)
library(Seurat) #Seurat v.3.0
library(RColorBrewer)
library(plotly)
BiocManager::install("monocle")
library(monocle3)
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
BiocManager::install("slingshot")
library(slingshot)

#LOAD SEURAT OBJECT
load(file = paste(dir, "Seurat_Objects/SingleCell_Data", sep = ""))

#Slingshot
#reference: https://rnabioco.github.io/cellar/posts/2021-10-26-class-3/, 
#https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html

pca_coords <- Embeddings(object = Seurat, reduction = "harmony")
clusters <- Seurat$harmony_FCT

#random start
slingshot <- slingshot(data = pca_coords, clusterLabels = clusters)
print(slingshot)

#start at NPCs
slingshot_2 <- slingshot(data = pca_coords, clusterLabels = clusters,
                             start.clus = 'NPC 1 (Inhibitory NPC')

#print(slingshot_2)

#Add Slingshot Pseudotime Metadata to Seurat Object
slingshot_pseudotime <- data.frame(slingPseudotime(slingshot_2))
Seurat <- AddMetaData(Seurat, metadata = slingshot_pseudotime)

#visualize each pseudotime curve
FeaturePlot(Seurat, features = "curve1", cols = viridis::magma(20))

FeaturePlot(Seurat, features = "curve2", cols = viridis::magma(20))

FeaturePlot(Seurat, features = "curve3", cols = viridis::magma(20))

#visualize the psuedotime across the timepoints for the cell types (for each curve)
meta_data <- Seurat[[]]
meta_data$harmony_FCT <- factor(meta_data$harmony_FCT, levels = unique(list(meta_data$harmony_FCT)))

# Subest to only cells with value for curve 1
meta_data <- meta_data %>%
  dplyr::filter(!is.na(curve1))

density_plot <- ggplot2::ggplot(data = meta_data,
                                ggplot2::aes(x = curve1,
                                             y = harmony_FCT,
                                             fill = harmony_FCT)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_manual(values = sample_colors)

print(density_plot)

#genes that correlate with pseudotime (30 min run time): only test the top 2000 variable features

variable_features_pattern <- paste0("^", VariableFeatures(Seurat), "$")
genes_use <- grep(paste0(variable_features_pattern, collapse="|"),
                  rownames(Seurat))

# fit negative binomial GAM
SCE <- fitGAM(counts = GetAssayData(object = Seurat,
                                        slot = "counts"),
                  sds = slingshot_2,
                  genes = genes_use)

saveRDS(SCE, paste("Harmony/Specific/...."))

#run an association test to identify what genes correlate best with pseudotime.
pseudotime_genes <- associationTest(Seurat, lineages = TRUE)
head(pseudotime_genes)


# Rank by p_value 1
topgenes <- rownames(pseudotime_genes[order(pseudotime_genes$pvalue_1), ])[1:100]

# Get the information for curve 1
cell_info <- Seurat[["curve1"]]

cell_info <- cell_info %>%
  dplyr::filter(!is.na(curve1))

# Get the information for all cells
heatdata <- GetAssayData(object = Seurat, slot = "data")

heatdata <- heatdata[rownames(heatdata) %in% topgenes,
                     colnames(heatdata) %in% rownames(cell_info)]

# Rank by pseudotime order 
heatdata <- heatdata[ , order(cell_info$curve1)]

sample_info <- Seurat[["harmony_FCT"]]
sample_info <- sample_info[colnames(heatdata) , ]

harmony_FCT <- unique(Seurat $harmony_FCT)
sample_colors <- hue_pal()(length(harmony_FCT))
cell_type_order <- match(levels(Seurat@meta.data[["harmony_FCT"]]), metaLevels("harmony_FCT", Seurat))
names(sample_colors) <- "harmony_FCT"

color_list <- list("harmony_FCT" = sample_colors)

heatmap_scale <- t(scale(t(as.matrix(heatdata)), scale = TRUE))

# Colors for heatmap (from the ArchR package)
blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                         "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")
                         
heatmap_scale <- ifelse(heatmap_scale > 2.5, 2.5, heatmap_scale)
heatmap_scale <- ifelse(heatmap_scale < -2.5, -2.5, heatmap_scale)

pheatmap(heatmap_scale, cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = sample_info,
         annotation_colors = color_list, color = blueYellow,
         border_color = NA, clustering_method = "complete")

dimred <- Seurat@reductions$umap@cell.embeddings
clustering <- Seurat$harmony_FCT
counts <- as.matrix(data@assays$RNA@counts[data@assays$RNA@var.features, ])

#Monocle
#reference: http://cole-trapnell-lab.github.io/monocle-release/docs/#differential-expression-analysis

#Take genes that are expressed in at least 5% of cells
top_genes <- detectGenes(Seurat, min_expr = 0.1)
fData(top_genes)$use_for_ordering <- fData(top_genes)$num_cells_expressed > 0.05 * ncol(top_genes)


clustering_DEG_genes <-
  differentialGeneTest(top_genes[HSMM_expressed_genes,],
                       fullModelFormulaStr = '~Cluster',
                       cores = 1)

#Take the top 1000 genes
HSMM_ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

HSMM_myo <-
  setOrderingFilter(HSMM_myo,
                    ordering_genes = HSMM_ordering_genes)

HSMM_myo <-
  reduceDimension(HSMM_myo, method = 'DDRTree')

HSMM_myo <-
  orderCells(HSMM_myo)

HSMM_myo <-
  orderCells(HSMM_myo, root_state = GM_state(HSMM_myo))

plot_cell_trajectory(HSMM_myo, color_by = "Hours")



