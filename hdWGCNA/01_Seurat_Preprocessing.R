library(Seurat)

pfc_Seurat <- readRDS("/output/Post_mortem_brain_validation/PFC_427/Astrocytes.rds")

#Seurat Metadata processing
pfc_Seurat_metadata <- pfc_Seurat@meta.data

#ROSMAP Subset
metadata_subset <- rosmap_metadata %>%
  dplyr::filter(projid %in% pfc_Seurat_metadata$projid)

#Plot APOE genotype distribution
metadata_plot_df <- metadata_subset %>%
  dplyr::filter(!is.na(apoe_genotype)) %>%
  mutate(apoe_genotype = factor(apoe_genotype, levels = c(22,23,24,33,34,44)))

ggplot(metadata_plot_df, aes(x = apoe_genotype)) +
  geom_bar() +
  theme_minimal() +
  labs(x = "APOE Genotype", y = "Individual Count", fill = "", title = "Individual Distribution")

#Filter the metadata to contain only APOE 33, 34, 44
metadata_subset_filtered <- metadata_subset %>%
  dplyr::filter(apoe_genotype %in% c(33, 34, 44))

pfc_Seurat_metadata_filtered <- pfc_Seurat_metadata %>%
  dplyr::filter(projid %in% metadata_subset_filtered$projid)

pfc_Seurat_metadata_filtered$read <- rownames(pfc_Seurat_metadata_filtered)

pfc_Seurat_metadata_final <- merge(pfc_Seurat_metadata_filtered, metadata_subset, by = "projid")
rownames(pfc_Seurat_metadata_final) <- pfc_Seurat_metadata_filtered$read


pfc_count_mtx <- pfc_Seurat@assays$RNA@counts
dim(pfc_count_mtx) # 33538 X 149558

#Filter count matrix for only 33, 34 44 individuals
pfc_count_mtx_filtered <- pfc_count_mtx[, colnames(pfc_count_mtx) %in% rownames(pfc_Seurat_metadata_final)]
dim(pfc_count_mtx_filtered) #33538 X 128285



########Create the Seurat object
#Set the Seurat to V3 for downstream compatibility 
options(Seurat.object.assay.version = "v3")
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


#SCT Normalization
options(future.globals.maxSize = 50 * 1024^3)
pfc_seurat_obj_apoe_normalized <- SCTransform(pfc_seurat_obj)

pfc_seurat_obj_apoe_final <- FindVariableFeatures(pfc_seurat_obj_apoe_normalized)


#Data Scaling
all.genes <- rownames(pfc_seurat_obj_apoe_final)
pfc_seurat_obj_apoe_final <- ScaleData(pfc_seurat_obj_apoe_final, features = all.genes)

#PCA
pfc_seurat_obj_apoe_final <- RunPCA(pfc_seurat_obj_apoe_final, features = VariableFeatures(object = pfc_seurat_obj_apoe_final))

#Run Harmony
pfc_seurat_obj_apoe_final <- RunHarmony(pfc_seurat_obj_apoe_final, c("Study", "projid"))

ElbowPlot(pfc_seurat_obj_apoe_final)

#UMAP
pfc_seurat_obj_apoe_final <- RunUMAP(
  pfc_seurat_obj_apoe_final, 
  reduction = "harmony",
  assay = "SCT",
  dims = 1:20,
  n.neighbors = 100, 
  min.dist = 0.8, 
  seed.use = 2525
)

DimPlot(pfc_seurat_obj_apoe_final, reduction = "umap", group.by = "apoe_genotype") + ggtitle("")


pfc_seurat_obj_apoe_final <- FindNeighbors(
    pfc_seurat_obj_apoe_final, 
    reduction = "harmony", 
    dims = 1:20
)

pfc_seurat_obj_apoe_final <- FindClusters(pfc_seurat_obj_apoe_final, resolution = 0.02)

unique(pfc_seurat_obj_apoe_final$cell_type_high_resolution) #Levels: Ast CHI3L1, Ast DPP10, Ast GRM3

new_cluster_names <- c(
    "0" = "Ast CHI3L1",
    "1" = "Ast DPP10",
    "2" = "Ast GRM3"
)
pfc_seurat_obj_apoe_final <- RenameIdents(pfc_seurat_obj_apoe_final, new_cluster_names)
pfc_seurat_obj_apoe_final$renamed_clusters <- Idents(pfc_seurat_obj_apoe_final)

DimPlot(pfc_seurat_obj_apoe_final, reduction = "umap", group.by = "renamed_clusters") + ggtitle("")

VlnPlot(pfc_seurat_obj_apoe_final, features = "nFeature_RNA", pt.size = 0) + NoLegend() +
  labs(x = "") +
  ggtitle("nFeature_RNA")

VlnPlot(pfc_seurat_obj_apoe_final, features = "nCount_RNA", pt.size = 0) + NoLegend() +
  labs(x = "") +
  ggtitle("nCount_RNA")
