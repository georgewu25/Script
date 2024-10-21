library(Seurat)

pfc_Seurat <- readRDS("/output/Post_mortem_brain_validation/PFC_427/Astrocytes.rds")

#Seurat Metadata processing
pfc_Seurat_metadata <- pfc_Seurat@meta.data

metadata <- read.csv("/output/Post_mortem_brain_validation/ROSMAP_clinical.csv")

pfc_Seurat_metadata_filtered <- metadata %>%
  dplyr::filter(projid %in% metadata_subset_filtered$projid)

pfc_Seurat_metadata_filtered$read <- rownames(pfc_Seurat_metadata_filtered)

pfc_Seurat_metadata_final <- merge(pfc_Seurat_metadata_filtered, metadata, by = "projid")
rownames(pfc_Seurat_metadata_final) <- pfc_Seurat_metadata_filtered$read

pfc_Seurat_metadata_final$major_cell_type <- "Ast"

pfc_count_mtx <- pfc_Seurat@assays$RNA@counts
dim(pfc_count_mtx) # 33538 X 254721

pfc_count_mtx_filtered <- pfc_count_mtx[, colnames(pfc_count_mtx) %in% rownames(pfc_Seurat_metadata_final)]
dim(pfc_count_mtx_filtered) #33538 X 14765


#Recreate the Seurat object
options(Seurat.object.assay.version = "v5")
pfc_seurat_obj <- CreateSeuratObject(counts = pfc_count_mtx_filtered, meta.data = pfc_Seurat_metadata_final)

#Need to see feature distribution between genotypes
pfc_seurat_obj_apoe <- pfc_seurat_obj
Idents(pfc_seurat_obj_apoe) <- "apoe_genotype"

VlnPlot(pfc_seurat_obj_apoe, features = "nCount_RNA", pt.size = 0) + NoLegend() +
  labs(x = "APOE Genotype") +
  ggtitle("PFC nCount_RNA")

VlnPlot(pfc_seurat_obj_apoe, features = "nFeature_RNA", pt.size = 0) + NoLegend() +
  labs(x = "APOE Genotype") +
  ggtitle("PFC nFeature_RNA")

FeatureScatter(pfc_seurat_obj_apoe, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#SCT Normalization
options(future.globals.maxSize = 50 * 1024^3)
pfc_seurat_obj_apoe_normalized <- SCTransform(pfc_seurat_obj)

pfc_seurat_obj_apoe_final <- FindVariableFeatures(pfc_seurat_obj_apoe_normalized)


#Data Scaling
all.genes <- rownames(pfc_seurat_obj_apoe_final)
pfc_seurat_obj_apoe_final <- ScaleData(pfc_seurat_obj_apoe_final, features = all.genes)

#PCA
pfc_seurat_obj_apoe_final <- RunPCA(pfc_seurat_obj_apoe_final, features = VariableFeatures(object = pfc_seurat_obj_apoe_final))


DimPlot(pfc_seurat_obj_apoe_final, reduction = "pca") + NoLegend()

DimHeatmap(pfc_seurat_obj_apoe_final, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(pfc_seurat_obj_apoe_final)

#Cell Clustering
#pfc_seurat_obj_apoe_final <- FindNeighbors(pfc_seurat_obj_apoe_final, dims = 1:10)
#pfc_seurat_obj_apoe_final <- FindClusters(pfc_seurat_obj_apoe_final, resolution = 0.5)

#UMAP
pfc_seurat_obj_apoe_final <- RunUMAP(pfc_seurat_obj_apoe_final, dims = 1:10)

#Run Harmony
pfc_seurat_obj_apoe_final <- RunHarmony(pfc_seurat_obj_apoe_final, c("Study", "projid"))
