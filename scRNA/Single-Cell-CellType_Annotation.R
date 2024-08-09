library(dplyr)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(scales)
library(cowplot)
library(Seurat) 
library(SingleR)
library(RColorBrewer)

load(file = paste(prev,"/harmonized_clustered", sep =""))

#Labeled clustering (after Harmony)
DimPlot(Seurat_SCT, reduction = 'harmony_res_0.5_umap', group.by ='seurat_clusters', label = TRUE)

Seurat_SCT$harmony_clusters <- Seurat_SCT$seurat_clusters

Idents(Seurat_SCT) <- Seurat_SCT$harmony_FCT

DotPlot(
  Seurat_SCT,
  features = x,
  cols = c("light grey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 10,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA) +ggtitle("ependymal")

#Marker Featureplot

Final_Marker_List = list(Final_Markers_NPC, Final_Markers_Neuron, Final_Markers_Astrocyte, Final_Markers_OPC, Final_Markers_VLMC)
Final_Marker_List_Names = c("NPC_Markers", "Neuron_Markers", "Astrocyte_Markers", "OPC_Markers", "VLMC_Markers")

for (x in Final_Marker_List) {
  for (y in x) {
    
    name = Final_Marker_List_Names[count]
    pdf(file = paste(dir, "/Feature_Plots/", name, "_", y, "_feature_plot.pdf", sep = ""), width = 15, height = 15)
    print(FeaturePlot(Seurat_SCT, reduction = "harmony_res_0.5_umap", features = y, label = TRUE, label.size = 6))
    dev.off()
  }
}


feats <- list(Final_Markers_NPC, Final_Markers_Neuron[1:13], Final_Markers_Neuron[14:25], Final_Markers_Astrocyte, Final_Markers_OPC, Final_Markers_VLMC)
marker_names <- list("NPC", "Neuron_1", "Neuron_2", "Astrocyte", "OPC", "VLMC")
count = 0

for (y in feats) {
  
  count = count + 1
  name = marker_names[count]
 
  #dotplot
  pdf(paste(dir, "/Dot_Plots/", name, "_dot_plot.pdf", sep = ""), width = 25, height = 15)
  
   print(DotPlot(
    Seurat_SCT,
    features = y,
    cols = c("light grey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 10,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA))
  
  dev.off()
} 

#Marker Stacked Violin Plot

#code source: https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)),  
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

#Generate Stacked Violin Plot of Expression

#Micrcoglia Check
StackedVlnPlot(obj = Seurat_SCT, features = Microglia_Markers[6:11])

print(FeaturePlot(Seurat_SCT, reduction = "harmony_res_0.5_umap", features = "TMEM119", label = TRUE, label.size = 6))

#DEGs in Every Cluster
#Use FindMarkers to get the DEGs in each cluster
#Marker genes for each cluster are identified using the Wilcoxon rank sum test implemented in the FindMarkers function, 
#with a logFC threshold of 0.25 and Bonferroni correction for multiple testing
#Marker genes for each cluster were manually compared to known marker genes in the PangloDB Database

Idents(Seurat_SCT) <- Seurat_SCT$harmony_clusters

for (x in rev(unique(Seurat_SCT$harmony_clusters))){
  
  #get positive DEGs in each cluster
  cluster.markers <- FindMarkers(TCW, ident.1 = x, min.pct = 0.25, test.use = "wilcox", only.pos = TRUE)
  
  write.csv(cluster.markers %>% slice_max(n=2000, order_by = avg_log2FC), paste(dir, "/Cluster_Markers/cluster_",x,"_markers.csv", sep = ""))
}


#SINGLER ANNOTATIONS: Automatically map query data to reference atlas

library(scRNAseq)
library(SingleR)
install.packages("tidyverse")
library(tidyverse)
library(Seurat)

load(paste(ref_dir, "/sc_dev_cortex_geschwind/raw_counts_mat.rdata", sep = ""))
raw.counts.mat <- as.matrix(raw_counts_mat)
reffy <- CreateSeuratObject(raw_counts_mat)

#load metadata
ref_metadata <- read.csv(file = paste(ref_dir, "/sc_dev_cortex_geschwind/cell_metadata.csv", sep = ""))
rownames(ref_metadata) <- ref_metadata[,1]
ref_metadata[,1] <- NULL

#append metadata to reference object
reffy <- AddMetaData(reffy, ref_metadata)
reffy <- reffy[,!is.na(reffy$Cluster)]

#convert to a single cell experiment
reffo <- as.SingleCellExperiment(reffy)
metadata(reffo)<- ref_metadata

#normalize data
library(scuttle)
reffo <- logNormCounts(reffo)
reffo@metadata$Cluster <- as.matrix(reffo@metadata$Cluster)

#perform SingleR labeling
results <- SingleR(test = as.SingleCellExperiment(TCW), ref = reffo, labels = reffo@colData$Cluster)
save(results, file = "Polioudakis_SingleR_Results.R")

#append resulting cell type labels to metadata
Seurat_SCT$Polioudakis_SingleR <- results$labels

#save UMAP visualizations
DimPlot(Seurat_SCT, reduction = "harmony_res_0.5_umap", group.by = "Polioudakis_SingleR", label = TRUE, label.size = 6)

#save annotated seurat object
save(Seurat_SCT, file = paste(dir, "/TCW_SingleR_Annotated.R", sep = ""))

#Reference Atlas: UCSC Speir reference

labels <- fread(paste(ref_dir, "/exprMatrix_Speir.tsv", sep = ""), nrows = 1)
labels <- labels[,-1]
col_list <- sample(names(labels), 20000) 
col_list <- as.list(col_list)
col_list <- append(col_list, "gene", 0)
col_list <- unlist(col_list)
exp_matrix <- fread(paste(ref_dir,"/exprMatrix_Speir.tsv",sep =""), header = TRUE, select = col_list)

reffy <- as.data.frame(exp_matrix)
row.names(reffy) <- reffy[,1]
reffy[,1] <- NULL
reffy <- as.data.frame(reffy)

#get gene list
gene_list = rownames(reffy)
type(gene_list)
new_gene_list = list()
for (x in gene_list) {
  x = gsub("\\|.*","",x)
  new_gene_list <- append(new_gene_list, x)
}

rownames(reffy) <- new_gene_list

reffy <- CreateSeuratObject(reffy)

ref_metadata <- read_tsv(file = paste(ref_dir, "/meta_Speir.tsv", sep = ""))
ref_metadata <- as.matrix(ref_metadata)
rownames(ref_metadata) <- ref_metadata[,1]
ref_metadata <- ref_metadata[,-1]
ref_metadata <- as.data.frame(ref_metadata)

reffy <- AddMetaData(reffy, ref_metadata)
reffy <- reffy[,!is.na(reffy$Cluster)]

reffo <- as.SingleCellExperiment(reffy)
metadata(reffo)<- ref_metadata 

library(scuttle)
reffo <- logNormCounts(reffo)
reffo@metadata$Cluster <- as.matrix(reffo@metadata$Cluster)

results <- SingleR(test = as.SingleCellExperiment(TCW), ref = reffo, labels = reffo@colData$Type)
save(results, file = "TCW_Speir_SingleR_Subtype_results") 

Seurat_SCT$Speir_Subtype <- results$labels
Seurat_SCT$Speir_Type <- results$labels


#save UMAP visualizations
pdf(file = paste(dir, "/SingleR/Speir_Subtype_Labelled_UMAP.pdf", sep = ""), width = 10, height = 10)
DimPlot(Seurat_SCT, reduction = "harmony_res_0.5_umap", group.by = "Speir_Subtype", label = TRUE, label.size = 6)
dev.off()


#Without Harmony 

DimPlot(Seurat_SCT, reduction = "harmony_40_nn_0.4res_umap", group.by = "harmony_FCT", pt.size = .1, label = TRUE) 
DimPlot(Seurat_SCT, reduction = "umap", group.by = "no_harmony_FCT", pt.size = .1, label = TRUE) + NoLegend()

#Specific Cell Types

Seurat_SCT$Final_Cell_Type <- Seurat_SCT$harmony_clusters
Idents(Seurat_SCT)<- Seurat_SCT$Final_Cell_Type

newIdent <- "Astrocyte 3"
names(newIdent) <- "Astrocyte 1 1"
Seurat_SCT <- RenameIdents(object = Seurat_SCT, newIdent)

Seurat_SCT$Final_Cell_Type <- Seurat_SCT@active.ident
Seurat_SCT <- Seurat_SCT[, !(TCW$Final_Cell_Type %in% "12")]

#General Cell Types
Seurat_SCT$Gen_Final_Cell_Type <- Seurat_SCT$Final_Cell_Type
Idents(Seurat_SCT) <- Seurat_SCT$Final_Cell_Type 

newIdent <- "GPC"
names(newIdent) <- "NPC 2 (Cycling NPC + GPC)"
Seurat_SCT <- RenameIdents(object = TCW, newIdent)

Seurat_SCT$Final_Cell_Type <- Seurat_SCT@active.ident

#Pick cell type colors
FCT_colors <- c("#E42237","#EE7D47", "#F4A652", "#F5CA63", "#DE40AA",
                     "#D97FCA", "#B590D2", "#A2B2DB", "#249BC0", "#219CDE","#008BCC",
                     "#DA9DC6", "#254289") 
 
Gen_FCT_colors <- c("#E42237", "#EE7D47", "#F4A652", "#F5CA63", 
                    "#DE40AA", "#249BC0", "#DA9DC6", "#254289") 
  
#legend order
FCT_order = rev(c("NPC", "Inhibitory Neuron", "Excitatory Neuron", "Mature Excitatory Neuron", "GPC",
                  "Astrocyte 1", "Astrocyte 2", "Astrocyte 3", "Astrocyte 4", "Astrocyte 5", "Astrocyte 6","OPC + Oligodendrocyte", "VLMC"))

gen_FCT_order <- rev(c("NPC", "Inhibitory Neuron", "Excitatory Neuron", "Mature Excitatory Neuron",
                       "GPC", "Astrocyte","OPC + Oligodendrocyte", "VLMC"))

png(file = paste(dir, "/Final_Cell_Type_UMAP.png", sep = ""), width = 1000, height = 1000)
DimPlot(Seurat_SCT, reduction = "harmony_res_0.5_umap", group.by = "Final_Cell_Type", label = TRUE, label.size = 7, order = FCT_order, cols = FCT_colors) + 
  ggtitle("UMAP of Cell Subtypes in hiPSC MCC") +
  theme(plot.title = element_text(size = 30, face = "bold", hjust=0.6),
        axis.text.x=element_text(angle= 0, hjust=1, size= 10),
        axis.text.y=element_text(angle= 0, hjust=1, size= 10),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +  
  guides(color = guide_legend(override.aes = list(size=6), ncol=1) )
dev.off()

png(file = paste(dir, "/Gen_Final_Cell_Type_UMAP.png", sep = ""), width = 1000, height = 1000)
DimPlot(Seurat_SCT, reduction = "harmony_res_0.5_umap", group.by = "Gen_Final_Cell_Type", label = TRUE, label.size = 7, order = gen_FCT_order, cols = Gen_FCT_colors) + 
  ggtitle("UMAP of Cell Types in hiPSC MCC") +
  theme(plot.title = element_text(size = 30, face = "bold", hjust=0.6),
        axis.text.x=element_text(angle= 0, hjust=1, size= 10),
        axis.text.y=element_text(angle= 0, hjust=1, size= 10),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +  
  guides(color = guide_legend(override.aes = list(size=6), ncol=1) )
dev.off()


#BY Individual
png(file = paste(dir, "/Final_Individual_UMAP.png", sep = ""), width = 1000, height = 1000)
DimPlot(Seurat_SCT, reduction = "harmony_res_0.5_umap", group.by = "APOE_Genotype", label = TRUE, label.size = 7, cols = Individual_colors) + 
  ggtitle("UMAP of Individual Distrubtion") +
  theme(plot.title = element_text(size = 30, face = "bold", hjust=0.6),
        axis.text.x=element_text(angle= 0, hjust=1, size= 10),
        axis.text.y=element_text(angle= 0, hjust=1, size= 10),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +  
  guides(color = guide_legend(override.aes = list(size=6), ncol=1) )
dev.off()


#With Harmony

DimPlot(Seurat_SCT, reduction = "harmony_40_nn_0.4res_umap", group.by = "", pt.size = .1, label = TRUE) 
DimPlot(Seurat_SCT, reduction = "no_harmony_50_nn_0.4_res", group.by = "no_harmony_FCT", pt.size = .1, label = TRUE) 

Idents(Seurat_SCT)<- Seurat_SCT$harmony_FCT

newIdent <- "NPC 3"
names(newIdent) <- "13"
Seurat_SCT <- RenameIdents(object = Seurat_SCT, newIdent)

#update metadata
Seurat_SCT$harmony_FCT <- Seurat_SCT@active.ident

DimPlot(Seurat_SCT, reduction = 'umap', group.by = 'harmony_FCT', pt.size = .1, label = TRUE) 


#General Annotation

Seurat_SCT$general_harmony_FCT <- Seurat_SCT$harmony_FCT

Idents(Seurat_SCT) <- Seurat_SCT$general_harmony_FCT

newIdent <- "NPC"
names(newIdent) <- "NPC 3"
Seurat_SCT <- RenameIdents(object = Seurat_SCT, newIdent)

Seurat_SCT$general_harmony_FCT <- Seurat_SCT@active.ident

DimPlot(Seurat_SCT, reduction = "harmony_40_nn_0.4res_umap", group.by = "general_harmony_FCT", pt.size = .1, label = TRUE) 
DimPlot(Seurat_SCT, reduction = "harmony_40_nn_0.4res_umap", group.by = "harmony_FCT", pt.size = .1, label = TRUE)
DimPlot(Seurat_SCT, reduction = "harmony_40_nn_0.4res_umap", group.by = "individual", pt.size = .1, label = TRUE)

Seurat_SCT$general_no_harmony_FCT <- Seurat_SCT$no_harmony_FCT

Idents(Seurat_SCT) <- Seurat_SCT$general_no_harmony_FCT

newIdent <- "NPC"
names(newIdent) <- "NPC 3"
Seurat_SCT <- RenameIdents(object = TCW, newIdent)

Seurat_SCT$general_no_harmony_FCT <- Seurat_SCT@active.ident

DimPlot(Seurat_SCT, reduction = "no_harmony_50_nn_0.4_res_umap", group.by = "general_no_harmony_FCT", pt.size = .1, label = TRUE) 
DimPlot(Seurat_SCT, reduction = "no_harmony_50_nn_0.4_res_umap", group.by = "no_harmony_FCT", pt.size = .1, label = TRUE) 
DimPlot(Seurat_SCT, reduction = "no_harmony_50_nn_0.4_res_umap", group.by = "individual", pt.size = .1, label = TRUE) 


#Barplots
#BiocManager::install("dittoSeq")
library(dittoSeq)

#Color Palettes!
Individual_colors = c("grey", "#bd7ebe")

#retrieve harmony FCT colors
require(scales)
identities <- levels(Seurat_SCT$Final_Cell_Type)
FCT_palette <- hue_pal()(length(identities))

#retrieve harmony FCTorder
#FCT_cluster_order <- match(levels(Seurat_SCT@meta.data[["Final_Cell_Type"]]), metaLevels("Final_Cell_Type", Seurat_SCT))

FCT_order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)

#retrieve harmony FCT colors
require(scales)
identities <- levels(Seurat_SCT$Gen_Final_Cell_Type)
Gen_FCT_palette <- hue_pal()(length(identities))

#retrieve harmony FCTorder
#Gen_FCT_cluster_order <- match(levels(TCW@meta.data[["Gen_Final_Cell_Type"]]), metaLevels("Gen_Final_Cell_Type", Seurat_SCT))

Gen_FCT_order <- c(1, 2, 3, 4, 5, 6, 7, 8)
###

#official id retrieve colors
require(scales)
identities <- unique(Seurat_SCT$Individual)
Ind_color_palette <- hue_pal()(length(identities))

#official id retrieve order
Ind_cluster_order <- match(unique(Seurat_SCT@meta.data[["Individual"]]), metaLevels("Individual", Seurat_SCT))

###

#Harmony APOE distribution split by cell type
png(file = paste(dir, "/Bar_Plots/Barplot_APOE_Genotype_Distribution_by_FCT.png", sep = "") , width = 5, height = 5)
dittoBarPlot(object = Seurat_SCT,
             var = "Final_Cell_Type", 
             group.by = "APOE_Genotype",
             color.panel = FCT_colors,
             var.labels.reorder = FCT_order,
             y.breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +ggtitle(expression(paste("Cell Subtype Proportion Changes in ", italic('APOE'), " 44")))
dev.off()

#Harmony Cell Type split by APOE Genotype
png(file = paste(dir, "/Bar_Plots/Barplot_Cell_Type_Distribution_by_APOE_Genotype.png", sep = "") , width = 6, height = 5)
dittoBarPlot(object = Seurat_SCT,
             var = "APOE_Genotype", 
             group.by = "Final_Cell_Type",
             color.panel = APOE_colors,
             y.breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +ggtitle("Barplot of APOE Genotype Distribution Across Cell Types")
dev.off()


#Cell Type split by Individual
png(file = paste(dir, "/Bar_Plots/Barplot_Individual_Distribution_by_Cell_Type.png", sep = "") , width = 6, height = 5)
dittoBarPlot(object = Seurat_SCT,
             var = "Individual", 
             group.by = "Final_Cell_Type",
             color.panel = Individual_colors,
             y.breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +ggtitle("Barplot of Individual Distribution Across Cell Types")
dev.off()


#WILCOXON RANK SUM TEST for statistically sig cell type proportion differences
unique(TCW$hash.ID)
library(data.table)

wilcoxon_table <- data.table()
new_row <- data.table("sample" = "example", 
                      "cell_type" = "cell", 
                      "APOE 33" = 0,
                      "APOE 44" = 0)
wilcoxon_table <- rbindlist(list(wilcoxon_table, new_row))

for (sample in unique(TCW$hash.ID)){
  for (geno in unique(TCW$APOE_Genotype)){
    
    sub <- subset(TCW, hash.ID %in% sample)  
    
    total_cells = length(colnames(sub))

    if (geno %in% unique(sub$APOE_Genotype)){
      
      sub <- subset(sub, APOE_Genotype %in% geno)
     
    for (cell in unique(sub$Gen_Final_Cell_Type)){
      
      if (cell %in% unique(sub$Gen_Final_Cell_Type)){
                    
        set <- subset(sub, (Gen_Final_Cell_Type %in% cell))  
        set <- colnames(set) #cells
        
        if (identical(geno, "APOE 33")){
          new_row <- data.table("sample" = sample, 
                                 "cell_type" = cell, 
                                 "APOE 33" = length(set)/total_cells,
                                 "APOE 44" = 0)
        }
        
        else{
            new_row <- data.table("sample" = sample, 
                                  "cell_type" = cell,
                                  "APOE 33" = 0,
                                  "APOE 44" = length(set)/total_cells)

        }
      
        wilcoxon_table <- rbindlist(list(wilcoxon_table, new_row))
      
    }
    }
    }
    
  }
}

mat_ex <- wilcoxon_table[wilcoxon_table$cell_type %in% "Mature Excitatory Neuron",] 
neu_ex <- wilcoxon_table[wilcoxon_table$cell_type %in% "Excitatory Neuron",] 
neu_in <- wilcoxon_table[wilcoxon_table$cell_type %in% "Inhibitory Neuron",] 
npc <- wilcoxon_table[wilcoxon_table$cell_type %in% "NPC",] 
gpc <- wilcoxon_table[wilcoxon_table$cell_type %in% "GPC",]
astro <- wilcoxon_table[wilcoxon_table$cell_type %in% "Astrocyte",] 
VLMC <- wilcoxon_table[wilcoxon_table$cell_type %in% "VLMC",]
