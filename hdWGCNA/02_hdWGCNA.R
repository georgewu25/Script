library(Seurat)
library(tidyverse)
library(cowplot)
#library(patchwork)
library(WGCNA)
library(hdWGCNA, lib = "../mwu/R/")
library(glmGamPoi, lib = "../mwu/R/")


pfc_seurat_obj_apoe_final <- readRDS("/output/Post_mortem_brain_validation/PFC_427/Seurat_Final.rds")

#Selected after testing various parameters
K=c(25, 50, 75, 100)
M=c(10, 30, 60, 90)
network_type = "signed"

hdwgcna_obj <- SetupForWGCNA(
  pfc_seurat_obj_apoe_final,
  gene_select = "fraction",
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "PFC427_Ast_WGCNA"
)

#Constructing MetaCells
for (i in K) {
  for (j in M) {
    hdwgcna_obj <- MetacellsByGroups(
      hdwgcna_obj,
      group.by = c("cell_type_high_resolution", "projid"),
      ident.group = 'cell_type_high_resolution',
      reduction = "harmony",
      assay = "SCT",
      slot = "data",
      k = i, #nearest neighbors
      #min_cells = m, #minimum number of cells in a metacell group
      max_shared = j, #maximum number of shared cells
      #target_metacells = t,
      verbose = FALSE
    )
    
    seurat_obj <- NormalizeMetacells(hdwgcna_obj)
    seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
    seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
    seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='projid')
    seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:15)
    
    pdf(paste0("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Post_mortem_brain_validation/PFC_427/Metacell_Testing/", i, "k_", j, "m_", "default_t.pdf"), width = 15, height = 15)
    print(DimPlotMetacells(seurat_obj, group.by='cell_type_high_resolution') + umap_theme() + ggtitle(paste0(i, "k, ", j, "m, ", "default t")))
    dev.off()
  }
}

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
