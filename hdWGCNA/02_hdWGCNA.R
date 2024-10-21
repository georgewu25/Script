library(Seurat)
library(tidyverse)
library(cowplot)
#library(patchwork)
library(WGCNA)
library(hdWGCNA, lib = "../mwu/R/")
library(glmGamPoi, lib = "../mwu/R/")


pfc_seurat_obj_apoe_final <- readRDS("/output/Post_mortem_brain_validation/PFC_427/Seurat_Final.rds")

#Run Harmony
pfc_seurat_obj_apoe_final <- RunHarmony(pfc_seurat_obj_apoe_final, c("Study", "projid"))

K=75
M=60
Reduction = "harmony"
network_type = "signed"

hdwgcna_obj <- SetupForWGCNA(
  pfc_seurat_obj_apoe_final,
  gene_select = "fraction",
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "PFC_Ast_WGCNA"
)

#Constructing MetaCells
hdwgcna_obj <- MetacellsByGroups(
  hdwgcna_obj,
  group.by = c("major_cell_type", "projid"),
  ident.group = 'major_cell_type',
  reduction = Reduction,
  assay = "SCT",
  slot = "data",
  k = K, #nearest neighbors
  min_cells = 100, #minimum number of cells in a metacell group
  max_shared = M, #maximum number of shared cells
  target_metacells = 1000,
  max_iter = 5000,
  verbose = FALSE
)

hdwgcna_obj <- NormalizeMetacells(hdwgcna_obj)


#Using full Seurat object since Metacell is not performed due to using subset and processed data
hdwgcna_obj_final <- SetDatExpr(
  hdwgcna_obj,
  group_name = "Ast",
  group.by= "major_cell_type",
  assay = 'SCT',
  slot = 'data'
)

hdwgcna_obj_final <- TestSoftPowers(
  hdwgcna_obj_final,
  networkType = network_type
)

# plot_list <- PlotSoftPowers(hdwgcna_obj_final)
# wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(hdwgcna_obj_final)

#Construct co-expression network:
hdwgcna_obj_tom <- ConstructNetwork(
  hdwgcna_obj_final,
  tom_name = "PFC_Ast",
  consensus = FALSE,
  overwrite_tom = TRUE,
  wgcna_name = NULL,
  blocks = NULL,
  maxBlockSize=100000,
  minBlockSize=0,
  minModuleSize = 30,
  corType = "bicor",
  maxPOutliers=0.10,
  networkType = network_type,
  TOMType = network_type,
  TOMDenom = "min",
  scaleTOMs = TRUE,
  deepSplit = 4,
  pamStage = FALSE,
  saveConsensusTOMs = TRUE,
  verbose=4
)
