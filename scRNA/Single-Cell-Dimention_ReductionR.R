#Normalization, Scaling, PCA, & Preliminary Clustering
load(file = "Seurat_Final")

Seurat <- CellCycleScoring(Seurat_Final, s.features = cc.genes.updated.2019$s.genes,
                        g2m.features = cc.genes.updated.2019$g2m.genes,
                        set.ident = TRUE,
                        search=TRUE)

Seurat[["CC.Difference"]] <- Seurat$G2M.Score - Seurat$S.Score

all.genes <- rownames(Seurat)
Seurat <- NormalizeData(Seurat, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")

#find variable features
Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)

#scale only mt
Seurat_scaled <- ScaleData(Seurat, assay = "RNA", vars.to.regress = "percent.mt")
Seurat_scaled  <- RunPCA(Seurat_scaled, features = VariableFeatures(object = Seurat_scaled))

VizDimLoadings(Seurat_scaled, dims = 1:2, reduction = "pca")
DimHeatmap(Seurat_scaled, dims = 1, cells = 500, balanced = TRUE)
RidgePlot(Seurat_scaled, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

Seurat_scaled <- RunPCA(Seurat_scaled, features = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))
DimPlot(Seurat_scaled)

Seurat_scaled_all <- ScaleData(Seurat, features = all.genes, assay = "RNA", vars.to.regress = c("percent.mt", "nCount_RNA", "CC.Difference"), verbose = FALSE)
Seurat_scaled_all  <- RunPCA(Seurat_scaled_all, features = VariableFeatures(object = Seurat_scaled_all))

VizDimLoadings(Seurat_scaled_all, dims = 1:2, reduction = "pca")
DimHeatmap(Seurat_scaled_all, dims = 1, cells = 500, balanced = TRUE)
RidgePlot(Seurat_scaled_all, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

Seurat_scaled_all <- RunPCA(Seurat_scaled_all, features = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))
DimPlot(Seurat_scaled_all)

#run sctransform
BiocManager::install("glmGamPoi")
Seurat_SCT <- SCTransform(Seurat_scaled_all, method = "glmGamPoi", vars.to.regress = c("percent.mt", "nCount_RNA", "CC.Difference"), verbose = FALSE)

save(Seurat_SCT)

Seurat_SCT <- RunPCA(Seurat_SCT, verbose = FALSE)
Seurat_SCT <- RunUMAP(Seurat_SCT, dims = 1:50, verbose = FALSE)
Seurat_SCT <- FindNeighbors(Seurat_SCT, dims = 1:50, verbose = FALSE)
Seurat_SCT <- FindClusters(Seurat_SCT, resolution = 0.4, verbose = FALSE)

DimPlot(Seurat_SCT , group.by = "seurat_clusters", label=TRUE)

Seurat_SCT$initial_seurat_clustering_0.4 <- Seurat_SCT$seurat_clusters
Seurat_SCT@reductions$initial_seurat_clustering_umap <- Seurat_SCT@reductions$umap

VlnPlot(Seurat_SCT, features = c("nFeature_RNA"), pt.size = 0.5)
VlnPlot(Seurat_SCT, features = c("nCount_RNA"), pt.size = 0.5)
VlnPlot(Seurat_SCT, features = c("percent.mt"), pt.size = 0.5)

save(Seurat_SCT)

#Remove Ambient RNA: SoupX & DecontX --> USE RAW COUNTS!
#preserve original counts prior to ambient RNA removal
Seurat_SCT[["original_counts"]] <- CreateAssayObject(counts = Seurat_SCT@assays$RNA@counts) #raw counts


#SoupX 
#install.packages('SoupX')
library(SoupX)

#raw background
raw_cell_ranger_out <- paste(source, "compbio/cellranger/outs/raw_feature_bc_matrix/", sep = '')
raw <- Read10X(data.dir = raw_cell_ranger_out, gene.column = 2)

#modify genes in raw
genes <- rownames(Seurat_SCT@assays$RNA@counts)
raw <- raw[genes,]

#run SoupX algo (this package is what I ended up using instead of DecontX)
sc <- SoupChannel(raw, Seurat_SCT@assays$RNA@counts)
sc <- setClusters(sc, Seurat_SCT $initial_seurat_clustering_0.5)
sc <- autoEstCont(sc, doPlot = FALSE)
soup_out <- adjustCounts(sc, roundToInt = TRUE) #adjusted counts matrix

Seurat_SCT[["SoupX_counts"]] <- CreateAssayObject(counts = soup_out)
Seurat_SCT@assays$RNA@counts <- soup_out #USED

#calculate percent of mRNA removed
percent_removed_soup <- (1 - (sum(Seurat_SCT@assays$RNA@counts) / sum(Seurat_SCT@assays$original_counts@counts)))*100

#DecontX (alternative method)
library(celda)
library(singleCellTK)

Seurat_sparse <- Seurat_SCT@assays$RNA@counts

#raw background
raw_cell_ranger_out <- paste(source, "/compbio/cellranger/outs/raw_feature_bc_matrix/", sep = '')
raw <- Read10X(data.dir = raw_cell_ranger_out, gene.column = 2)

#modify genes in raw
genes <- rownames(Seurat_SCT@assays$RNA@counts)
raw <- raw[genes,]
decont_out <- decontX(Seurat_sparse, background = raw)

Seurat_SCT@assays$RNA@counts <- decont_out$decontXcounts
Seurat_SCT $decontX_contamination <- decont_out$contamination
Seurat_SCT[["DecontX_counts"]] <- CreateAssayObject(counts = decont_out$decontXcounts) 

#calculate percent of mRNA removed
percent_removed_decont <- ( 1 - (sum(Seurat_SCT@assays$DecontX_counts@counts) / sum(Seurat_SCT@assays$original_counts@counts)) )*100

#Plot contamination results
FeaturePlot(Seurat_SCT, features = "decontX_contamination")

#Clustering after ambient RNA Removal
Seurat_SCT <- RunPCA(Seurat_SCT, verbose = FALSE)
Seurat_SCT <- RunUMAP(Seurat_SCT, dims = 1:50, verbose = FALSE)
Seurat_SCT <- FindNeighbors(Seurat_SCT, dims = 1:50, verbose = FALSE)
Seurat_SCT <- FindClusters(Seurat_SCT, resolution = 0.6, verbose = FALSE)

DimPlot(Seurat_SCT, reduction = "umap", group.by = "SoupX_seurat_clustering_0.6", label = TRUE)

Seurat_SCT$SoupX_seurat_clustering_0.6 <- Seurat_SCT$seurat_clusters

Idents(Seurat_SCT) <- "SoupX_seurat_clustering_0.6"

VlnPlot(Seurat_SCT, features = c("nFeature_RNA"), pt.size = 0.5)
VlnPlot(Seurat_SCT, features = c("nCount_RNA"), pt.size = 0.5)
VlnPlot(Seurat_SCT, features = c("percent.mt"), pt.size = 0.5)

#Homotypic Doublet & Ambient RNA Removal (must be performed AFTER standard QC)
#pipeline: https://www.nature.com/articles/s41467-022-29212-9 (SingleCellTK Pipeline)

#Homotypic Doublet Removal: scDblFinder & DoubletFinder (intersection)

#scDbl
BiocManager::install("scDblFinder")
library(scDblFinder)

Seurat_scbdl <- scDblFinder(GetAssayData(Seurat_sparse, slot="counts"), clusters= Seurat_sparse$SoupX_seurat_clustering_0.6) 

Seurat_SCT$scDblFinder.score <- Seurat_scbdl$scDblFinder.score 
Seurat_SCT$scDblFinder.class <- Seurat_scbdl$scDblFinder.class
scDbl_doublets <- subset(Seurat_sparse, scDblFinder.class %in% "doublet") #168

#DoubletFinder
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

#pK identification
sweep.res.list <- paramSweep_v3(Seurat_SCT, PCs = 1:50, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]])) #optimal pK = 0.23

#pN = 0.25 (default)
#exp = 12% (high throughput) or 24% (standard) #USED 12

#Homotypic Doublet Estimation
annotations <- Seurat_SCT@meta.data$SoupX_seurat_clustering_0.6
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.12*nrow(Seurat_sparse@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#run doublet finder
Seurat_doubletfinder <- doubletFinder_v3(Seurat_SCT, 
                                         PCs = 1:50, 
                                         pN = 0.25, 
                                         pK = pK, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = FALSE, sct = FALSE)

Seurat_SCT$DF.classifications_0.25_0.23_841 <- Seurat_doubletfinder$DF.classifications_0.25_0.23_841

DF_doublets <- subset(Seurat_doubletfinder, DF.classifications_0.25_0.23_841 %in% "Doublet") #841

#Find intersecting doublets
select_doublets <- intersect(colnames(scDbl_doublets), colnames(DF_doublets)) #76 cells

#Remove the intersect cells
Seurat_SCT <- Seurat_SCT[, !(colnames(Seurat_SCT) %in% select_doublets)]

#Re-do prelim clustering
Seurat_SCT <- RunPCA(Seurat_SCT, verbose = FALSE)
Seurat_SCT <- RunUMAP(Seurat_SCT, dims = 1:50, verbose = FALSE)
Seurat_SCT <- FindNeighbors(Seurat_SCT, dims = 1:50, verbose = FALSE)
Seurat_SCT <- FindClusters(Seurat_SCT, resolution = 0.6, verbose = FALSE)

DimPlot(Seurat_SCT, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

Seurat_SCT$initial_post_qc_clustering_0.6 <- Seurat_SCT$seurat_clusters

#At this point, the cells that remain are all singlets and free from ambient RNA
save(Seurat_SCT, file = paste(dir,"/Fully_Filtered_Singlets.R", sep =""))

################################################################################
#Harmony
install.packages("ggthemes")
install.packages("harmony")
#devtools::install_github("immunogenomics/harmony", build_vignettes=TRUE)
library(harmony)

#perform Harmony integration: remove bias by individual
Seurat_SCT <- RunHarmony(Seurat_SCT, group.by.vars = "Individual", assay.use = "SCT", reduction.save = "harmony")
#converged after 2 iterations

nn_list = c(10, 20, 30, 40, 50)
res_list = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

for (x in nn_list) {
  for (y in res_list) {
    
    #cluster
    Seurat_SCT <- Seurat_SCT %>%
      RunUMAP(reduction = "harmony", n.neighbors = x, dims = 1:50) %>% 
      FindNeighbors(reduction = "harmony", dims = 1:50) %>%
      FindClusters(resolution = y) 
    
    #umap visualizations
    after <- DimPlot(Seurat_SCT, reduction = 'umap', group.by ='APOE_Genotype', cols = APOE_colors)
  
    pdf(file = paste(dir, "/Harmony_Clustering/UMAP_by_APOE_genotype_after_harmony_nn_",x,"_res_",y,"_.pdf", sep = ""))
    print(after)
    dev.off()
    
    after <- DimPlot(Seurat_SCT, reduction = 'umap', group.by ='Individual', cols = Individual_colors)
    
    pdf(file = paste(dir, "/Harmony_Clustering/UMAP_by_Individual_after_harmony_nn_",x,"_res_",y,"_.pdf", sep = ""))
    print(after)
    dev.off()
    
    after <- DimPlot(Seurat_SCT, reduction = 'umap', group.by ='seurat_clusters')
    
    pdf(file = paste(dir, "/Harmony_Clustering/UMAP_by_clusters_after_harmony_nn_",x,"_res_",y,"_.pdf", sep = ""))
    print(after)
    dev.off()
    
  }
}

#USED: res = 0.5, nn = 30.

Seurat_SCT <- Seurat_SCT %>%
  RunUMAP(reduction = "harmony", n.neighbors = 30, dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>%
  FindClusters(resolution = 0.5) 

save(Seurat_SCT, file = paste(dir,"/harmonized_clustered.R", sep =""))

Seurat_SCT@reductions$harmony_res_0.5_umap <- Seurat_SCT@reductions$umap

#Cell type annotation, DEA, fgsea, & network analysis will be conducted next
