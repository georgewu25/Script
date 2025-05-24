library(Seurat)

#HTO Matrix
cite_seq_out <- paste(source, "/citeseq_outs/umi_count/", sep = "")
HTO_matrix <- Read10X(data.dir = cite_seq_out, gene.column = 1)

#UMI Matrix: filtered by CellRanger
cell_ranger_out <- paste(source, "/compbio/cellranger/outs/filtered_feature_bc_matrix/", sep = '')
UMI_matrix <- Read10X(data.dir = cell_ranger_out, gene.column = 1)


#Generate TAG list
tag_table <- data.frame()

col_1 <- sapply(rownames(HTO_matrix)[-nrow(HTO_matrix)], function(x)strsplit(x, "-")[[1]][1])
col_2 <- sapply(rownames(HTO_matrix)[-nrow(HTO_matrix)], function(x)strsplit(x, "-")[[1]][2])
col_1 <- as.data.frame(as.vector(col_1))
col_2 <- as.data.frame(as.vector(col_2))

tag_table <- rbindlist(list(tag_table, col_2)) 
tag_table$HTO <- col_1 


#Generate white list: cells present in the UMI Matrix

whitelist <- sapply(colnames(UMI_matrix)[-ncol(UMI_matrix)], function(x)strsplit(x, "-")[[1]][1])
whitelist <- as.data.frame(as.vector(whitelist))

##############################################
#After re-doing Cite-Seq-Count

#UMI Matrix: filtered by CellRanger
cell_ranger_out <- paste(source, "/compbio/cellranger/outs/filtered_feature_bc_matrix/", sep = '')
UMI_matrix <- Read10X(data.dir = cell_ranger_out, gene.column = 1)
length(colnames(UMI_matrix)) #30697

#HTO Matrix
cite_seq_out <- paste(source, "/compbio/citeseq_outs/umi_count/", sep = "")
HTO_matrix <- Read10X(data.dir = cite_seq_out, gene.column = 1)


#Intersecting Barcodes: Cells Detected by Both
joint.bcs <- intersect(colnames(UMI_matrix), colnames(HTO_matrix)) #19,204

# Subset RNA and HTO counts by joint cell barcodes
UMI_matrix <- UMI_matrix[, joint.bcs]
HTO_matrix <- HTO_matrix[, joint.bcs]


Seurat <- CreateSeuratObject(counts = UMI_matrix)

#Log Normalization
Seurat <- NormalizeData(Seurat)


# Find and scale variable features
Seurat <- FindVariableFeatures(Seurat, selection.method = "mean.var.plot")
Seurat <- ScaleData(Seurat, features = VariableFeatures(TCW))

#Add hastag data
Seurat[["HTO"]] <- CreateAssayObject(counts = HTO_matrix)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Seurat <- NormalizeData(Seurat, assay = "HTO", normalization.method = "CLR")

#Demultiplex
Seurat <- HTODemux(Seurat, assay = "HTO", qstat -u = 0.995, seed = 100)

# Global classification results
table(Seurat$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(Seurat) <- "HTO_maxID"

RidgePlot(Seurat, assay = "HTO", features = rownames(Seurat[["HTO"]]), ncol = 9)

FeatureScatter(Seurat, feature1 = "hto_HTO-A", feature2 = "hto_HTO-B")

Idents(Seurat) <- "HTO_classification.global"

VlnPlot(Seurat, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

HTO_ID_list <- unique(Seurat $HTO_maxID)

Idents(TCW) <- Seurat$HTO_maxID

done_list = list()
for (x in HTO_ID_list){
  for (y in HTO_ID_list){
    if ((y != x) & !(y %in% done_list)){
      pdf(paste(dir, "/Demultiplexing/feat_scatter_",x,"_and_",y,".pdf", sep = ""), height = 10, width = 15)
      print(FeatureScatter(TCW, feature1 = x, feature2 = y))
      dev.off()
    }
  }
  done_list <- append(done_list, x)
}


#Remove Negative Cells
Seurat_filtered <- subset(Seurat, idents = "Negative", invert = TRUE)


#Calculate a tSNE embedding of the HTO data: Singlets vs. Doublets
DefaultAssay(Seurat_filtered) <- "HTO"
Seurat_filtered <- ScaleData(Seurat_filtered, features = rownames(Seurat_filtered), verbose = FALSE)
Seurat_filtered <- RunPCA(Seurat_filtered, features = rownames(Seurat_filtered), approx = FALSE)
Seurat_filtered <- RunTSNE(Seurat_filtered, dims = 1:8, perplexity = 100)

DimPlot(Seurat_filtered)

HTOHeatmap(Seurat_filtered, assay = "HTO"))

# Extract the singlets
Seurat_filtered <- subset(Seurat_filtered, idents = "Singlet")

Further QC must be performed to remove homotypic doublets (doublets where the cells come from the same sample) and remove ambient RNA existing within droplets containing a single cell.

#LOWLY EXPRESSED GENE filtration: keep only the genes expressed in at least 3 cells.
raw_counts_mat <- as.data.frame(Seurat_filtered@assays$RNA@counts)
raw_counts_mat <- raw_counts_mat %>% mutate(num_cells_expressing =rowSums(.!=0)) #count non-zero cells

ggplot(raw_counts_mat, aes(num_cells_expressing)) + geom_density()

raw_counts_mat <- subset(raw_counts_mat, num_cells_expressing >= 3)
geneNames <- rownames(raw_counts_mat) #genes in >= 3 cells

Seurat_subset <- subset(Seurat_filtered, features = geneNames) 

raw_counts_mat <- as.data.frame(Seurat_subset@assays$RNA@counts)
raw_counts_mat <- raw_counts_mat %>% mutate(num_cells_expressing =rowSums(.!=0)) #count non-zero cells

ggplot(raw_counts_mat, aes(num_cells_expressing)) + geom_density()

#LOW QUALITY CELL filtration: remove cells with very little or very high expression 
#total number of reads (UMIs) detected per cell = library size

umi_counts <- as.data.frame(Seurat_subset$nCount_RNA)
umi_counts %>%
  dplyr::rename("nCount_RNA" = "TCW$nCount_RNA") %>%
  ggplot(aes(x = nCount_RNA)) + 
  geom_density() +
  theme_bw()

ggplot(umi_counts, aes(x = cell_index, y = nCount_RNA)) + geom_bar(stat='identity')+
  theme_bw()


#total number of genes detected per cell
gene_counts <- as.data.frame(Seurat_subset$nFeature_RNA)
gene_counts %>%
  rename("nFeature_RNA" = "Seurat_subset$nFeature_RNA") %>%
  ggplot(aes(x = nFeature_RNA)) + 
  geom_density() +
  theme_bw()

ggplot(gene_counts, aes(x = cell_index, y = nFeature_RNA)) + geom_bar(stat='identity')+
  theme_bw()

#MITO GENES
Idents(Seurat_subset) <- Seurat_subset$orig.ident
DefaultAssay(Seurat_subset) <- "RNA"

#label MT percent in each cell
Seurat_subset[["percent.mt"]] <- PercentageFeatureSet(Seurat_subset, pattern = "^MT-")

mt_counts <- as.data.frame(TCW$percent.mt)
mt_counts %>%
  dplyr::rename("percent.mt" = "TCW$percent.mt") %>%
  ggplot(aes(x = percent.mt)) + 
  geom_density() +
  theme_bw()


#Ribosomal percent per cell
Seurat_subset[["percent.ribo"]] <- PercentageFeatureSet(Seurat_subset, pattern =  "^RP[SL]")

ribo_counts <- as.data.frame(Seurat_subset$percent.ribo)
ribo_counts %>%
  dplyr::rename("percent.ribo" = "Seurat_subset$percent.ribo") %>%
  ggplot(aes(x = percent.ribo)) + 
  geom_density() +
  theme_bw()

#QC metric violin plots
vps<-VlnPlot(object = Seurat_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = "counts", pt.size = 0)
vp1<-vps[[1]] +geom_hline(yintercept = 200,colour="red")
vp2<-vps[[2]]+geom_hline(yintercept = 75000,colour="red") +geom_hline(yintercept = 500,colour="red")
vp3<-vps[[3]]+geom_hline(yintercept = 14,colour="red")

p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank()) +vp4 + theme(axis.text.x = element_blank())
plot(p)


#FeatureScatter is used to visualize feature-feature relationships
#Ensure linear relationship between nFeature_RNA and nCount_RNA
plot1 <- FeatureScatter(Seurat_subset, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Seurat_subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(Seurat_subset, feature1 = "nCount_RNA", feature2 = "percent.ribo")


#Remove Cells with High MT Percentage
#MiQC FILTERING: REMOVE low quality cells and set MT threshold
#vingette: https://bioconductor.org/packages/devel/bioc/vignettes/miQC/inst/doc/miQC.html
#BiocManager::install("flexmix")
library(miQC)
#remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
library(flexmix)

#Run MiQC Algorithm
#note: posterior cutoff = the posterior probability of a cell being part of the compromised distribution, a number between 0 and 1.
#Any cells below the appointed cutoff will be marked to keep. Defaults to 0.75.
Seurat_subset <- RunMiQC(Seurat_subset, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.97, model.slot = "flexmix_model") 

PlotMiQC(Seurat_subset, color.by = "miQC.probability") + ggplot2::scale_color_gradient(low = "grey", high = "purple")

PlotMiQC(Seurat_subset, color.by = "miQC.keep")

#Run MiQC Filtering
Seurat_subset2 <- subset(Seurat_subset, miQC.keep == "keep")

#Check new mito.pct vs. nCount_RNA scatter
plot1 <- FeatureScatter(Seurat_subset2, feature1 = "nCount_RNA", feature2 = "percent.mt")

#Remove Additional Outliers
#Only include cells with at least 200 genes expressed
Seurat_subset2 <- subset(Seurat_subset2, subset = nFeature_RNA > 200)

#Only include cells with at least 500 molecules expressed
Seurat_subset2 <- subset(Seurat_subset2, subset = nCount_RNA > 500) #number of genes detected in each cell

#Remove the upper tail observed on the "nCount_RNA" violin plot (99.5% percentile)
ub <- quantile(Seurat_subset2[["nCount_RNA"]]$nCount_RNA, probs = 0.99)
Seurat_subset2 <- Seurat_subset2[, Seurat_subset2[["nCount_RNA"]] < ub]

#FeatureScatter is used to visualize feature-feature relationships
#Ensure linear relationship between nFeature_RNA and nCount_RNA
plot1 <- FeatureScatter(Seurat_subset2, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(Seurat_subset2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(Seurat_subset2, feature1 = "nCount_RNA", feature2 = "percent.ribo")


#QC metric violin plots
vps<-VlnPlot(object = Seurat_subset2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = "counts", pt.size = 0) 
vp1<-vps[[1]] +labs(x = "TCW_14311_run1")
vp2<-vps[[2]]+geom_hline(yintercept = 75000,colour="red")
vp3<-vps[[3]]+geom_hline(yintercept = 14,colour="red")
vp4<-vps[[4]]+geom_hline(yintercept = 30,colour="red")

p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank())
plot(p)


#QC metric violin plots with ribo
vps<-VlnPlot(object = Seurat_subset2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, slot = "counts", pt.size = 0)
vp1<-vps[[1]] +geom_hline(yintercept = 200,colour="red")
vp2<-vps[[2]]+geom_hline(yintercept = 75000,colour="red") +geom_hline(yintercept = 500,colour="red")
vp3<-vps[[3]]+geom_hline(yintercept = 14,colour="red")
vp4<-vps[[4]]+geom_hline(yintercept = 30,colour="red")

p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank()) +vp4 + theme(axis.text.x = element_blank())
plot(p)

save(Seurat_subset2, file = "Seurat_Final")

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

#run SoupX algo
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

load(file = paste(dir_data, "Prelim_Analysis/Cell_Type_Annotation/FCT", sep = ""))

#Annotation with Gene Symbol
t2g <- read_tsv(paste(dir_data, "SingleR_Reference_Data/gencode.v26.annotation_1.gtf", sep = ""), skip = 5, col_names = FALSE)

t2g %<>%
  dplyr::rename(feature = X3, attribute = X9) %>%
  filter(feature == "gene") %>%
  mutate(gene_id = str_match(attribute, "gene_id \"(\\S+)\";")[, 2],
         gene_name = str_match(attribute, "gene_name \"(\\S+)\";")[, 2],) %>%
  select(gene_id, gene_name)

t2g %>% filter(is.na(gene_id)) %>% nrow

#Get metadata table
metadata <- Percent_Expressing(Seurat_SCT,
  features = feats,
  threshold = 0,
  group_by = "Final_Cell_Type",
  split_by = NULL,
  entire_object = FALSE,
  slot = "data",
  assay = "SCT")


#Get cell types
harmony_cell_type_list <- as.vector(unique(Seurat$harmony_FCT))

#Exclude cell types that do not have enough cells in each group for comparison
harmony_cell_type_list <- harmony_cell_type_list[harmony_cell_type_list %in% c("NPC 3", "VLMC", "OPC + Oligodendrocyte", "Astrocyte 2") == FALSE]

#Run DEA (MAST & Wilcoxon)
for (set in sets) { #set = All, individual_1, or individual_2
  
  #create data table for DEA method comparison
  DT = data.table()
  wilcox_stat_meth_DT = data.table()
  seurat_mast_meth_DT = data.table()
  
  for (x in harmony_cell_type_list){
    
    Idents(set) <- set$harmony_FCT
    sub <- subset(set, idents = x)
    
    Idents(set) <- set$Individual
    
    Idents(sub) <- sub$APOE_Genotype 
    
    #Seurat Wilcoxon via FindMarkers function (only.pos = FALSE) to find differentially expressed genes (both up & down).
    Wilcoxon <- FindMarkers(sub, ident.1 = "APOE 44", ident.2 = "APOE 33", logfc.threshold = 0, min.pct = 0.1, assay="SCT", test.use = "wilcox", only.pos = FALSE)
   
      genes_to_test <- rownames(Wilcoxon)
      
      Wilcoxon_stat <- wilcoxauc(sub, 'APOE_Genotype', seurat_assay = "SCT",  groups_use = c('APOE 44', 'APOE 33')) 
      Wilcoxon_stat <- as.data.table(Wilcoxon_stat)
      
      #get correct statistic & auc
      for(i in 1:nrow(Wilcoxon_stat)){ #each row
        
        gene <- Wilcoxon_stat$feature[i] 
        
        stat_44 <- Wilcoxon_stat[Wilcoxon_stat$feature == gene & Wilcoxon_stat$group == 'APOE 44']$statistic 
        stat_33 <- Wilcoxon_stat[Wilcoxon_stat$feature == gene & Wilcoxon_stat$group == 'APOE 33']$statistic 
        
        auc_44 <- Wilcoxon_stat[Wilcoxon_stat$feature == gene & Wilcoxon_stat$group == 'APOE 44']$auc 
        auc_33 <- Wilcoxon_stat[Wilcoxon_stat$feature == gene & Wilcoxon_stat$group == 'APOE 33']$auc
        
        statistic <- min(stat_33, stat_44) #minimum of u-stat_1 & u-stat_2
        auc <- min(auc_33, auc_44)
        
        Wilcoxon_stat[Wilcoxon_stat$feature == gene & Wilcoxon_stat$group == 'APOE 44']$statistic <- statistic
        Wilcoxon_stat[Wilcoxon_stat$feature == gene & Wilcoxon_stat$group == 'APOE 44']$auc <- auc
        
      }
      
      Wilcoxon_stat <- Wilcoxon_stat[Wilcoxon_stat$feature %in% genes_to_test]
      Wilcoxon_stat <- Wilcoxon_stat[Wilcoxon_stat$group %in% 'APOE 44']
      
      Wilcoxon_stat <- Wilcoxon_stat %>% 
        dplyr::rename('gene' = 'feature') %>%  #left is what is new #rename(new_column_name = old_column_name)
        dplyr::rename('p_val' = 'pval') %>% 
        dplyr::rename('BH_p_val_adj' = 'padj') #benjamini hochberg
      Wilcoxon_stat <- Wilcoxon_stat[, ('group'):=NULL] 
      
      #BC p value adjustment
      Wilcoxon_stat$BC_p_val_adj <- p.adjust(Wilcoxon_stat$p_val, method = "bonferroni") 
      
      #get avg log 2 fold change measurement from seurat wilcoxon --> logFC to avg_log2FC
      Wilcoxon_stat <- Wilcoxon_stat %>% slice_max(n = 1000000000, order_by = p_val)
      Wilcoxon <- Wilcoxon %>% slice_max(n = 1000000000, order_by = p_val)
      
      #Wilcoxon <- subset(Wilcoxon, rownames(Wilcoxon) %in% Wilcoxon_stat$gene)
      Wilcoxon_stat$avg_log2FC <- Wilcoxon$avg_log2FC
      Wilcoxon_stat <- as.data.frame(Wilcoxon_stat) #only df have rownames
      rownames(Wilcoxon_stat) <- Wilcoxon_stat$gene
      Wilcoxon_stat$gene <- NULL
      
      # MAST DEA (applies hurdle model & removes variation due to individual)
      #latent variable = individual
      if (identical(set, Seurat)){
        MAST <- FindMarkers(sub, ident.1 = "APOE 44", ident.2 = "APOE 33", logfc.threshold = 0, min.pct = 0.1, assay="RNA", test.use = "MAST", latent.vars = "Individual", only.pos = FALSE)
      } 
      else{
        MAST <- FindMarkers(sub, ident.1 = "APOE 44", ident.2 = "APOE 33", logfc.threshold = 0, min.pct = 0.1, assay="RNA", test.use = "MAST", only.pos = FALSE)
      }
      
      MAST <- MAST %>% 
        dplyr::rename('BC_p_val_adj' = 'p_val_adj') #bonferroni correction
        
        #BH p value adjustment
        MAST$BH_p_val_adj <- p.adjust(MAST$p_val, method = "hochberg") 
      
      print("done MAST")
      
    if (identical(set, Seurat)){
      y = "All"
    } 
    
    if (identical(set, Individual_1)){
      y = "Individual 1"
    } 
    
    if (identical(set, Individual_2)){
      y = "Individual 2"
    } 
          
    write.csv(Wilcoxon_stat %>% slice_max(n = 10000000, order_by = avg_log2FC), paste("Harmony/Specific/",y,"/",x,"/wilcox_stat_dea/",x,"_wilcoxon_stat_dea.csv", sep = ''))

    write.csv(MAST %>% slice_max(n = 10000000, order_by = avg_log2FC), paste("Harmony/Specific/",y,"/",x,"/seurat_mast_dea/",x,"seurat_mast_dea.csv", sep = ''))

    #Annotation gene name
    
    for (d in list(Wilcoxon_stat, MAST)){
      
      res_ordered <- d

      annotation <- t2g[match(rownames(res_ordered), t2g$gene_name),]
      resorderedAnnotated <- cbind(res_ordered, annotation)
      
      if (identical(d, Wilcoxon_stat)){
        w = "wilcox_stat_dea"
      }
      if (identical(d, MAST)){
        w = "seurat_mast_dea"
      }
      
      write.csv(as.data.frame(resorderedAnnotated), file= paste("Harmony/Specific/",y,"/",x,"/",w,"/harmony_",x,"_",w,"_Results_GTFAnnotated.csv", sep = ''))
      
      nodup <- resorderedAnnotated
      nodup$absvalFC <- abs(nodup$avg_log2FC)
     
      nodup <- nodup[order(nodup$gene_name,-nodup$absvalFC),]
      nodup <- nodup[!duplicated(nodup$gene_name),]
      write.csv(as.data.frame(nodup), file= paste("Harmony/Specific/",y,"/",x,"/",w,"/harmony_",x,"_",w,"_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
       
      res <- nodup #table of genes (exp > 10% cells, ranked by logFC E44 vs. E33) #this has been annotated
      
      res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
      
      #FDR of 0.05 
      res_table_pval_thres <- res_tbl %>% 
        mutate(threshold = p_val < 0.05) #significant DEGs

      res_table_pval_thres_only <- res_table_pval_thres[res_table_pval_thres$p_val < 0.05,] #only look at sig DEGs based on p-value < 0.05.
      res_table_pval_thres_only <- res_table_pval_thres_only[!is.na(res_table_pval_thres_only$p_val),] 
      res_table_pval_thres_only <- as.data.frame(res_table_pval_thres_only)
      write.csv(as.data.frame(res_table_pval_thres_only), file= paste("Harmony/Specific/",y,"/",x,"/",w,"/harmony_",x,"_",w,"_significant_0.05_degs_by_p_val_only_results.csv", sep = ''))
      
      #BH: Benjamini Hochberg
      res_table_BH_padj_thres <- res_tbl %>%
        mutate(threshold = BH_p_val_adj < 0.05) #significant DEG
     
      #Select Sig DEGS only BH p val adj 
      res_table_BH_padj_thres_only <- res_table_BH_padj_thres[res_table_BH_padj_thres$BH_p_val_adj < 0.05,] #only look at significant DEGs based on BH padj-value < 0.05
      res_table_BH_padj_thres_only <- res_table_BH_padj_thres_only[!is.na(res_table_BH_padj_thres_only$BH_p_val_adj),] 
      res_table_BH_padj_thres_only <- as.data.frame(res_table_BH_padj_thres_only)
      write.csv(as.data.frame(res_table_BH_padj_thres_only), file= paste("Harmony/Specific/",y,"/",x,"/",w,"/harmony_",x,"_",w,"_significant_0.05_degs_by_BH_padj_val_only_results.csv", sep = ''))
      
      #BC: Bonferroni Correction
      res_table_BC_padj_thres <- res_tbl %>%
        mutate(threshold = BC_p_val_adj < 0.05)
      
      #Select Sig DEGS only BC p val adj 
      res_table_BC_padj_thres_only <- res_table_BC_padj_thres[res_table_BC_padj_thres$BC_p_val_adj < 0.05,] #only look at significant DEGs based on BC padj-value < 0.05
      res_table_BC_padj_thres_only <- res_table_BC_padj_thres_only[!is.na(res_table_BC_padj_thres_only$BC_p_val_adj),] 
      res_table_BC_padj_thres_only <- as.data.frame(res_table_BC_padj_thres_only)
      write.csv(as.data.frame(res_table_BC_padj_thres_only), file= paste("Harmony/Specific/",y,"/",x,"/",w,"/harmony_",x,"_",w,"_significant_0.05_degs_by_BC_padj_val_only_results.csv", sep = ''))
   
    } 
  }


#FGSEA Prep
harmony_cell_type_list <- as.vector(unique(Seurat$harmony_FCT)) #explore APOE DEGs within the specific cell types
harmony_cell_type_list <- harmony_cell_type_list[harmony_cell_type_list %in% c("NPC 3", "VLMC", "OPC + Oligodendrocyte", "Astrocyte 2") == FALSE]

gen_harmony_cell_type_list <- as.vector(unique(Seurat$general_harmony_FCT))
gen_harmony_cell_type_list <- gen_harmony_cell_type_list[gen_harmony_cell_type_list %in% c("VLMC", "OPC + Oligodendrocyte") == FALSE]

for (t in c("Specific", "General")){
  
  if (identical(t, "Specific")){
    #cell_list <- harmony_cell_type_list
    cell_list <- c("Astrocyte 1", "Astrocyte 3", "Astrocyte 5", "Astrocyte 6")
  }
  
  if (identical(t, "General")){
    #cell_list <- gen_harmony_cell_type_list
    cell_list <- "Astrocyte"
  }
  
  for (x in cell_list){ #x = cell type
    
    #annotated no dup files, no p val adj threshold, no logfc threshold: all degs found via preso wilcox method
    
    wilcox_1 <- read_csv(paste("Harmony/",t,"/individual_1/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
    wilcox_2 <- read_csv(paste("Harmony/",t,"/individual_2/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
    
    wilcox_1 <- as.data.frame(wilcox_1)
    wilcox_2 <- as.data.frame(wilcox_2)
    
    #get only those genes in both individual 1 and individual 2
    com_degs <- c(unique(wilcox_1$gene_name), unique(wilcox_2$gene_name))
    com_degs <- as.data.table(com_degs)
    com_degs <- com_degs[duplicated(com_degs)] #gene names only present in both
    com_degs <- com_degs[!(is.na(com_degs))] #remove NA
    
    #build table
    colnames(com_degs)<- "gene_name"
    com_degs[, statistic := 0] #wilcoxon u-stat
    com_degs[, auc := 0] #wilcoxon auc --> for magnitude of fgsea stat
    com_degs[, avg_log2FC := 0] #log2fc --> for sign of fgsea stat
    com_degs[, fgsea_stat := 0] #fgsea_stat 
    com_degs[, gene_id := ""] #gene annotation 
    

    if (length(com_degs$gene_name) != 0){
      
      for (gene_test in com_degs$gene_name){ #extract stats for each deg from each ind
        
        stat_1 <- (subset(wilcox_1, gene_name == gene_test))$statistic #individual 1 u-stat
        stat_2 <- (subset(wilcox_2, gene_name == gene_test))$statistic #individual 2 u-stat
        auc_1 <- (subset(wilcox_1, gene_name == gene_test))$auc #individual 1 auc
        auc_2 <- (subset(wilcox_2, gene_name == gene_test))$auc #individual 2 auc
        logFC_1 <- (subset(wilcox_1, gene_name == gene_test))$avg_log2FC #individual 1 avg_log2FC
        logFC_2 <- (subset(wilcox_2, gene_name == gene_test))$avg_log2FC #individual 2 avg_log2FC
        gene_anno <- (subset(wilcox_1, gene_name == gene_test))$gene_id #gene annotation
        
        x_1 = -log10(auc_1) * sign(logFC_1)
        x_2 = -log10(auc_2) * sign(logFC_2)
          
        #if (sign(logFC_1) == sign(logFC_2)){ #the gene should be regulated in the same direction in both individuals
          
        com_degs[com_degs$gene_name == gene_test]$statistic <- mean(stat_1, stat_2) #mean(stat_1, stat_2) #take the avg of the stats --> for fgsea
        com_degs[com_degs$gene_name == gene_test]$fgsea_stat <- mean(x_1, x_2) #mean(stat_1, stat_2) #take the avg of the stats --> for fgsea
        com_degs[com_degs$gene_name == gene_test]$auc <- mean(auc_1, auc_2) #take the avg of the auc
        com_degs[com_degs$gene_name == gene_test]$avg_log2FC <- mean(logFC_1, logFC_2) #take the avg of the log fold change
        com_degs[com_degs$gene_name == gene_test]$gene_id <- gene_anno #gene annotation
        #}
        
        #else{ #else, remove the gene from the combined list
          #com_degs <- com_degs[gene_name != gene]
        #}
        
      } #for each deg
      
      write.csv(as.data.frame(com_degs), file= paste("Harmony/",t,"/combined/",x,"/fgsea/harmony_",x,"wilcox_combined_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
      
    } #if there are shared degs, build the table

  } #cell type
} #cell annotation type


#BH vs. BC Decision

for (t in c("Specific", "General")){
  
  if (identical(t, "Specific")){
    cell_list <- harmony_cell_type_list
  }
  
  if (identical(t, "General")){
    cell_list <- gen_harmony_cell_type_list
  }
  
  DT = data.table()
  wilcox_meth_DT = data.table()
  mast_meth_DT = data.table()
  bar_DT = data.table()
  
  for (x in cell_list){ #x = cell type
  
    for (d in list("wilcox_stat_dea", "seurat_mast_dea")){
          
          #WILCOXON
          if(identical(d, "wilcox_stat_dea")){ #wilcoxon: combine
            
            #all files have been gene annotated
            
            #Individual 1
            y = "individual_1"
      
              #BH: Benjamini Hochberg
          
              wilcox_res_table_BH_padj_thres_only_1 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BH_padj_val_only_results.csv", sep = ''))
              
              wilcox_res_table_BH_padj_thres_0.15_1 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BH_padj_val_and_log2fc_0.15_results.csv", sep = ''))
              
              wilcox_res_table_BH_padj_thres_0.25_1 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BH_padj_val_and_log2fc_0.25_results.csv", sep = ''))
              
              #BC: Bonferroni Correction
  
              wilcox_res_table_BC_padj_thres_only_1 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BC_padj_val_only_results.csv", sep = ''))
              
              wilcox_res_table_BC_padj_thres_0.15_1 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BC_padj_val_and_log2fc_0.15_results.csv", sep = ''))
              
              wilcox_res_table_BC_padj_thres_0.25_1 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BC_padj_val_and_log2fc_0.25_results.csv", sep = ''))
            
              
            #Individual 2
            y = "individual_2"
            
              #BH: Benjamini Hochberg
              wilcox_res_table_BH_padj_ <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BH_padj_val_only_results.csv", sep = ''))
             
              
              #BC: Bonferroni Correction
              wilcox_res_table_BC_padj<- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BC_padj_val_only_results.csv", sep = ''))
             
              
             
            #COMBINE WILCOXON individual 1 & 2 overlapping degs
              for (method in c("BH", "BC")){ #DEA method
                for (thresh in c("only")){ 
                    
                    wilcox_1 <- get(paste0("wilcox_res_table_",method,"_padj_thres_",thresh,"_1", sep='')) #individual 1
                    wilcox_2 <- get(paste0("wilcox_res_table_",method,"_padj_thres_",thresh,"_2", sep='')) #individual 2
          
                    wilcox_1 <- as.data.frame(wilcox_1)
                    wilcox_2 <- as.data.frame(wilcox_2)
                    
                    #get only those genes present in both lists
                    com_degs <- c(unique(wilcox_1$gene_name), unique(wilcox_2$gene_name)) #all gene names 
                    com_degs <- as.data.table(com_degs)
                    com_degs <- com_degs[duplicated(com_degs)] #gene names only present in both
                    com_degs <- com_degs[!(is.na(com_degs))] #remove NA
                    
                    #build table
                    colnames(com_degs)<- "gene_name"
                    com_degs[, statistic := 0] #wilcox u-stat
                    com_degs[, avg_log2FC := 0] #log2fc
                    com_degs[, ind_1_p_val := 0] #non-adjusted p-value
                    com_degs[, ind_2_p_val := 0] #non-adjusted p-value
                    com_degs[, ind_1_p_val_adj := 0] #adjusted p-value
                    com_degs[, ind_2_p_val_adj := 0] #adjusted p-value
                    com_degs[, gene_id := ""] #gene annotation
                    
                    
                    #fill table
                    if (length(com_degs != 0)){ #must be >0 shared degs
                      
                      for (gene_test in com_degs$gene_name){ #extract stats for each deg from each ind
 
                        stat_1 <- (subset(wilcox_1, gene_name == gene_test))$statistic #individual 1 u-stat
                        stat_2 <- (subset(wilcox_2, gene_name == gene_test))$statistic #individual 2 u-stat
                        logFC_1 <- (subset(wilcox_1, gene_name == gene_test))$avg_log2FC #individual 1 avg_log2FC
                        logFC_2 <- (subset(wilcox_2, gene_name == gene_test))$avg_log2FC #individual 2 avg_log2FC
                        gene_anno <- (subset(wilcox_1, gene_name == gene_test))$gene_id #gene annotation
                        
                          if (sign(logFC_1) == sign(logFC_2)){ #the gene should be regulated in the same direction in both individuals
                            
                            com_degs[com_degs$gene_name == gene_test]$statistic <- mean(stat_1, stat_2) #mean(stat_1, stat_2) #take the avg of the stats --> for fgsea
                            com_degs[com_degs$gene_name == gene_test]$avg_log2FC <- mean(logFC_1, logFC_2) #take the avg of the log fold change
                            com_degs[com_degs$gene_name == gene_test]$ind_1_p_val <- (subset(wilcox_1, gene_name == gene_test))$p_val #ind 1 p_val
                            com_degs[com_degs$gene_name == gene_test]$ind_2_p_val <- (subset(wilcox_2, gene_name == gene_test))$p_val #ind 2 p_val
                            
                            p_val_adj = paste(method, "_p_val_adj", sep = '')
                            
                            com_degs[com_degs$gene_name == gene_test]$ind_1_p_val_adj <- (subset(wilcox_1, gene_name == gene_test))[[paste0(method,"_p_val_adj")]] #ind 1 p_val_adj 
                            com_degs[com_degs$gene_name == gene_test]$ind_2_p_val_adj <- (subset(wilcox_2, gene_name == gene_test))[[paste0(method,"_p_val_adj")]] #ind 2 p_val_adj
                            
                            com_degs[com_degs$gene_name == gene_test]$gene_id <- gene_anno #gene annotation
                          }
                          
                          else{ #else, remove the gene from the combined list
                            com_degs <- com_degs[gene_name != gene_test]
                          }
                 
                      } #for each deg
     
                    } #if there are shared degs, build the table
                    
                    #re-assign table to variable name
                    
                    assign(paste0("wilcox_res_table_",method,"_padj_thres_",thresh,sep=''), com_degs)
                    
                } #avg_log2fc threshold
              } #for method
            } #if Wilcoxon
      
          #MAST 
          if(identical(d, "seurat_mast_dea")){
           
            y = "All"
            
            #BH: Benjamini Hochberg
            
            mast_res_table_BH_padj_thres_only <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_significant_0.05_degs_by_BH_padj_val_only_results.csv", sep = ''))
            
            mast_res_table_BH_padj_thres_0.15 <- read_csv(file= paste("Harmony/",t,"/",y,"/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_significant_0.05_degs_by_BH_padj_val_and_log2fc_0.15_results.csv", sep = ''))
            
            mast_res_table_BH_padj_thres_0.25 <- read_csv(file= paste("Harmony/",t,"/",y,"/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_significant_0.05_degs_by_BH_padj_val_and_log2fc_0.25_results.csv", sep = ''))
            
            #BC: Bonferroni Correction
            
            mast_res_table_BC_padj_thres_only <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_significant_0.05_degs_by_BC_padj_val_only_results.csv", sep = ''))
            
            mast_res_table_BC_padj_thres_0.15 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_significant_0.05_degs_by_BC_padj_val_and_log2fc_0.15_results.csv", sep = ''))
            
            mast_res_table_BC_padj_thres_0.25 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_significant_0.05_degs_by_BC_padj_val_and_log2fc_0.25_results.csv", sep = ''))
          
          } #if MAST 
        
        ###############
        #For each combination of dea method and logfc threshold, plot volcano plot & bp for each ell type
      
        for (m in c("BH", "BC")){ #p adj method
          for (thresh in c("only", "0.15", "0.25")){ #avg_log2fc
            
            if (identical(thresh, "only")){
              t_num = 0
            }
            
            if (identical(thresh, "0.15")){
              t_num = 0.15
            }
            
            if (identical(thresh, "0.25")){
              t_num = 0.25
            }
            
            if (identical(d, "wilcox_stat_dea")){
              j = "wilcox"
            }
            
            if (identical(d, "seurat_mast_dea")){
              j = "mast"
            }
            
            #VOLCANO PLOT: scatterplot of p-value vs. log fold change of DEGs --> visualize significance of DEGs with large FC
            #colors, downregulated = blue, upregulared = red, no change = gray
            
            #Volcano for MAST, wilcox individual 1, wilcox individual 2
        
            #Volcano plot, still display all degs --> don't use edited tables
            volcano <- get(paste(j,"_res_table_",m,"_padj_thres_only", sep = '')) %>% 
                mutate(threshold = case_when(
                avg_log2FC > t_num ~ "Up",
                avg_log2FC < -t_num  ~ "Down",
                ((avg_log2FC < t_num) & (avg_log2FC > -t_num)) ~ "No Change"))
                
            if (j == "mast"){
                #Volcano
                #p-val (not adjusted) with log2fc > 0.15
                pdf(paste("Harmony/",t,"/combined/",x,"/",d,"/",x,"_volcano_plot_p_val_adj_",m,"_log2fc_",thresh,".pdf", sep = ""), width = 10, height = 10) 
                p<- ggplot() + 
                  geom_point(volcano, mapping = aes(x = avg_log2FC, y = -log10(.data[[paste0(m,"_p_val_adj",sep='')]]), colour = threshold)) +
                  ggtitle(paste("Volcano Plot of ",x," Significant APOE E44 vs. E33 DEGs with ",expression("Log"[2]*"(FC)")," > ",thresh," & ",m," P Val Adj < 0.05", sep = "")) +
                  xlab(expression("Log"[2]*"(FC)")) + 
                  ylab(expression("-Log"[10]*"(Adjusted P-value)")) +
                  scale_y_continuous(limits = c(0,50)) +
                  scale_colour_manual(values = c("Up"="red", "Down"="blue", "No Change"="lightgray")) +
                  theme(legend.position = "none",
                        plot.title = element_text(size = rel(1.5), hjust = 0.5),
                        axis.title = element_text(size = rel(1.25))) + 
                  theme_bw()
                print(p)
                dev.off()
            }
            
            if (j == "wilcox"){ #here
              #Volcano
              
              #individual 1
              pdf(paste("Harmony/",t,"/combined/",x,"/",d,"/",x,"individual_1_volcano_plot_p_val_adj_",m,"_log2fc_",thresh,".pdf", sep = ""), width = 10, height = 10) 
              p<- ggplot() + 
                geom_point(volcano, mapping = aes(x = avg_log2FC, y = -log10(ind_2_p_val_adj), colour = threshold)) +
                ggtitle(paste("Volcano Plot of Individual 1 ",x," Significant APOE E44 vs. E33 DEGs with ",expression("Log"[2]*"(FC)")," > ",thresh," & ",m," P Val Adj < 0.05", sep = "")) +
                xlab(expression("Log"[2]*"(FC)")) + 
                ylab(expression("-Log"[10]*"(Adjusted P-value)")) +
                scale_y_continuous(limits = c(0,50)) +
                scale_colour_manual(values = c("Up"="red", "Down"="blue", "No Change"="lightgray")) +
                theme(legend.position = "none",
                      plot.title = element_text(size = rel(1.5), hjust = 0.5),
                      axis.title = element_text(size = rel(1.25))) + 
                theme_bw()
              print(p)
              dev.off()
              
              #individual 2
              pdf(paste("Harmony/",t,"/combined/",x,"/",d,"/",x,"individual_2_volcano_plot_p_val_adj_",m,"_log2fc_",thresh,".pdf", sep = ""), width = 10, height = 10) 
              p<- ggplot() + 
                geom_point(volcano, mapping = aes(x = avg_log2FC, y = -log10(ind_1_p_val_adj), colour = threshold)) +
                ggtitle(paste("Volcano Plot of Individual 2 ",x," Significant APOE E44 vs. E33 DEGs with ",expression("Log"[2]*"(FC)")," > ",thresh," & ",m," P Val Adj < 0.05", sep = "")) +
                xlab(expression("Log"[2]*"(FC)")) + 
                ylab(expression("-Log"[10]*"(Adjusted P-value)")) +
                scale_y_continuous(limits = c(0,50)) +
                scale_colour_manual(values = c("Up"="red", "Down"="blue", "No Change"="lightgray")) +
                theme(legend.position = "none",
                      plot.title = element_text(size = rel(1.5), hjust = 0.5),
                      axis.title = element_text(size = rel(1.25))) + 
                theme_bw()
              print(p)
              dev.off()
              
            }
            
            print('done Volcano')
        
            #BARPLOT: number of up & down regulated genes in this cell type for this dea test
        
              #Get subsetted data table by p adj value threshold for this p val corr method
              deg_summary <- get(paste(j,"_res_table_",m,"_padj_thres_",thresh, sep = ''))  
              
              deg_summary <- deg_summary %>%
                mutate(deg_type = case_when(
                  avg_log2FC > 0 ~ "Up",
                  avg_log2FC < 0  ~ "Down"))
              
              deg_summary <- deg_summary %>%
                mutate(deg_count = case_when(
                  deg_type == "Up" ~ 1,
                  deg_type == 'Down' ~ 1)) 
              
              #save barplot
              pdf(paste("Harmony/",t,"/combined/",x,"/",d,"/",x,"_DEG_barplot_p_val_adj_",m,"_log2fc_",thresh,".pdf", sep = ""), width = 10, height = 6) 
              p<-ggplot(deg_summary, aes(x = factor(deg_type), y = deg_count)) + 
                geom_bar(stat = "identity") +ggtitle(paste("Number of Significant Up & Down Regulated DEGs in", x))+ 
                xlab('Direction of Expression in APOE 44') + ylab('Number of APOE DEGs')+ 
                theme_bw()
              print(p)
              dev.off()
              
              print('done bar')
              
          #####################
          #large barplot: number of total degs found by each method for each cell type
          
          DT_new_row <- data.table("dea_method" = d, 
                                   "p_adj_method" = method, 
                                   "log2FC_threshold" = thresh,
                                   "cell_type" = x, 
                                   "Up" = sum(deg_summary$deg_type == 'Up'), 
                                   "Down"= sum(deg_summary$deg_type == 'Down'),
                                   "degs" = sum(deg_summary$deg_count))
          
          DT <- rbindlist(list(DT, DT_new_row))
          
          ###
          #up and down degs found for each cell type by specific method and each logfc threshold
          
          meth_DT_new_row_1 <- data.table("cell_type" = x, 
                                          "p_adj_method" = method, 
                                          "log2FC_threshold" = thresh,
                                          "direction" = "Up",
                                          "degs" = sum(deg_summary$deg_type == 'Up'))
          
          meth_DT_new_row_2 <- data.table("cell_type" = x, 
                                          "p_adj_method" = method,
                                          "log2FC_threshold" = thresh,
                                          "direction" = "Down",
                                          "degs" = sum(deg_summary$deg_type == 'Down'))
          
          
          if (identical(d, "seurat_mast_dea")){ 
            mast_meth_DT <- rbindlist(list(mast_meth_DT, meth_DT_new_row_1))
            mast_meth_DT <- rbindlist(list(mast_meth_DT, meth_DT_new_row_2))
          }
          
          if (identical(d, "wilcoxon_stat_dea")){ 
            wilcox_meth_DT <- rbindlist(list(wilcox_meth_DT, meth_DT_new_row_1))
            wilcox_meth_DT <- rbindlist(list(wilcox_meth_DT, meth_DT_new_row_2))
          }
          
          bar_DT_new_row <- data.table("cell_type" = x, 
                                       "p_adj_method" = method,
                                       "log2FC_threshold" = thresh,
                                       "dea_method" = d,
                                       "Up" = sum(deg_summary$deg_type == 'Up'),
                                       "Down" = sum(deg_summary$deg_type == 'Down'))
          
          bar_DT <- rbindlist(list(bar_DT, bar_DT_new_row))
          
        } #vol & bar log2fc thresh
      } #vol & bar dea method'
    } #dea method
  }#cell type 
  
  #after data tables have been filled for each cell type for both methods:
  for (m in c("BH", "BC")){ #p adj method
    for (thresh in c("only", "0.15", "0.25")){ #avg_log2fc

      print(length(DT[0]))
      
      corr_name <- m
      log2FC_thres <- thresh
      
      DT <- as.data.frame(DT)
      wilcox_stat_meth_DT <- as.data.frame(wilcox_stat_meth_DT)
      seurat_mast_meth_DT <- as.data.frame(seurat_mast_meth_DT)
      
      print(length(DT[0]))
      
      #subset only for this correction method
      DT_x <- subset(DT, p_adj_method == corr_name)
      DT_x <- subset(DT, log2FC_threshold == log2FC_thres)
      
      print(length(DT[0]))
      
      wilcox_stat_meth_DT_x <- subset(wilcox_stat_meth_DT, p_adj_method == corr_name)
      wilcox_stat_meth_DT_x <- subset(wilcox_stat_meth_DT, log2FC_threshold == log2FC_thres)
      
      seurat_mast_meth_DT_x <- subset(seurat_mast_meth_DT, p_adj_method == corr_name)
      seurat_mast_meth_DT_x <- subset(seurat_mast_meth_DT, log2FC_threshold == log2FC_thres)
      
      #Barplot of all DEA Methods by cell type
      identities <- levels(Seurat$harmony_FCT)
      identities <- identities[identities %in% DT$cell_type]
      harmony_FCT_palette <- hue_pal()(length(identities))

      DT_x$cell_type <- factor(DT_x$cell_type, levels = identities)
      
      pdf(paste("Harmony/",t,"/combined/DEA_Method_Comparison_Barplot_",corr_name,"_",log2FC_thres,"_Log2FC.pdf", sep = ""), width = 12, height = 8)
      p <- ggplot(DT_x, aes(x= dea_method, y= degs, fill= cell_type)) +
        geom_bar(stat="identity", position="dodge") +
        scale_fill_manual(values = harmony_FCT_palette)+
        theme(axis.text.x = element_text(angle = 90))+ 
        theme_bw()
      print(p)
      dev.off()
      
      print("done DT bar")

      wilcox_stat_meth_DT_x$direction <- factor(wilcox_stat_meth_DT_x$direction, levels = c("Up", "Down"))
      
      pdf(paste("Harmony/",t,"/combined/Wilcox_Stat_DEA_Results_by_Cell_Type_Barplot_",corr_name,"_",log2FC_thres,"_Log2FC.pdf", sep = ""), width = 12, height = 6)
      p <- ggplot(wilcox_meth_DT_x, aes(x= cell_type, y= degs, fill= direction)) +
        geom_bar(stat="identity", position="dodge") +
        scale_fill_manual(values = c("Up" = "red","Down" = "blue")) +
        theme(axis.text.x = element_text(angle = 90))+ 
        theme_bw()
      print(p)
      dev.off()
      
      seurat_mast_meth_DT$direction <- factor(seurat_mast_meth_DT$direction, levels = c("Up", "Down"))
  
      pdf(paste("Harmony/",t,"/combined/Seurat_Mast_DEA_Results_by_Cell_Type_Barplot_",corr_name,"_",log2FC_thres,"_Log2FC.pdf", sep = ""), width = 12, height = 6)
      p <- ggplot(mast_meth_DT_x, aes(x= cell_type, y= degs, fill= direction)) +
        geom_bar(stat="identity", position="dodge") +
        scale_fill_manual(values = c("Up" = "red","Down" = "blue")) +
        theme(axis.text.x = element_text(angle = 90))+ 
        theme_bw()
      print(p)
      dev.off()
      
      print("done barplots")
    
      
      DT_x <- as.data.table(DT_x)
      wilcox_stat_meth_DT_x <- as.data.table(wilcox_stat_meth_DT_x)
      seurat_mast_meth_DT_x <- as.data.table(seurat_mast_meth_DT_x)
      
      library(plotrix)
      
      cell_name_list <- list()
      for (name in rownames(bar_DT)){
        if (name == 'Astrocyte'){
          cell_name_list <- append(cell_name_list, 'Ast')
        }
        if (name == 'Astrocyte 1'){
          cell_name_list <- append(cell_name_list, 'Ast 1')
        }
        if (name == 'Astrocyte 2'){
          cell_name_list <- append(cell_name_list, 'Ast 2')
        }
        if (name == 'Astrocyte 3'){
          cell_name_list <- append(cell_name_list, 'Ast 3')
        }
        if (name == 'Astrocyte 4'){
          cell_name_list <- append(cell_name_list, 'Ast 4')
        }
        if (name == 'Astrocyte 5'){
          cell_name_list <- append(cell_name_list, 'Ast 5')
        }
        if (name == 'Astrocyte 6'){
          cell_name_list <- append(cell_name_list, 'Ast 6')
        }
        if (name == 'NPC'){
          cell_name_list <- append(cell_name_list, 'NPC')
        }
        if (name == 'NPC 1 (Inhibitory Neuron)'){
          cell_name_list <- append(cell_name_list, 'NPC 1')
        }
        if (name == 'NPC 2 (Cycling NPC + GPC)'){
          cell_name_list <- append(cell_name_list, 'NPC 2')
        }
        if (name == 'NPC 3'){
          cell_name_list <- append(cell_name_list, 'NPC 3')
        }
        if (name == 'Excitatory Neuron'){
          cell_name_list <- append(cell_name_list, 'Ex')
        }
        if (name == 'Inhibitory Neuron'){
          cell_name_list <- append(cell_name_list, 'In')
        }
        if (name == 'Mature excitatory Neuron'){
          cell_name_list <- append(cell_name_list, 'Mat ex')
        }
        if (name == 'OPC + Oligodendrocyte'){
          cell_name_list <- append(cell_name_list, 'OPC/Oli')
        }
        if (name == 'VLMC'){
          cell_name_list <- append(cell_name_list, 'VLMC')
        }
        
      }
      
      #edit bar DT
      mast_bar_DT <- bar_DT[bar_DT$dea_method] == 'seurat_mast_dea'
      mast_bar_DT <- mast_bar_DT[(mast_bar_DT$p_adj_method == corr_name) & (mast_bar_DT$log2FC_threshold == log2FC_thres)]
      mast_bar_DT$p_adj_method <- NULL
      mast_bar_DT$log2FC_threshold <- NULL
      mast_bar_DT <- as.data.frame(mast_bar_DT)
      rownames(mast_bar_DT) <- mast_bar_DT$cell_type
      mast_bar_DT$cell_type <- NULL
      
      wilcox_bar_DT <- bar_DT[bar_DT$dea_method] == 'wilcox_stat_dea'
      wilcox_bar_DT <- wilcox_bar_DT[(wilcox_bar_DT$p_adj_method == corr_name) & (wilcox_bar_DT$log2FC_threshold == log2FC_thres)]
      wilcox_bar_DT$p_adj_method <- NULL
      wilcox_bar_DT$log2FC_threshold <- NULL
      wilcox_bar_DT <- as.data.frame(wilcox_bar_DT)
      rownames(wilcox_bar_DT) <- wilcox_bar_DT$cell_type
      wilcox_bar_DT$cell_type <- NULL
      
      pdf(paste("Harmony/",t,"/combined/Seurat_Mast_DEA_Results_by_Cell_Type_Barplot_",corr_name,"_",log2FC_thres,"_Log2FC.pdf", sep = ""), width = 10, height = 6)
      p <- color2D.matplot(mast_bar_DT, 
                           show.values = 0.5,
                           axes = FALSE,
                           main = "DEGs",
                           cex.main = 2,
                           xlab = "",
                           ylab = "",
                           vcex = 2,
                           vcol = "black",
                           extremes = c("white", "#009999"),
                           Hinton = TRUE)
      axis(1, at = seq_len(ncol(bar_DT)) - 0.5,
           labels = colnames(bar_DT), tick = FALSE, cex.axis = 2)
      axis(2, at = seq_len(nrow(bar_DT)) -0.5,
           labels = rev(cell_name_list), tick = FALSE, las = 1, cex.axis = 2)
      print(p)
      
      dev.off()
      
      pdf(paste("Harmony/",t,"/combined/Wilcox_Stat_DEA_Results_by_Cell_Type_Barplot_",corr_name,"_",log2FC_thres,"_Log2FC.pdf", sep = ""), width = 10, height = 6)
      p <- color2D.matplot(wilcox_bar_DT, 
                           show.values = 0.5,
                           axes = FALSE,
                           main = "DEGs",
                           cex.main = 2,
                           xlab = "",
                           ylab = "",
                           vcex = 2,
                           vcol = "black",
                           extremes = c("white", "#009999"),
                           Hinton = TRUE)
      axis(1, at = seq_len(ncol(bar_DT)) - 0.5,
           labels = colnames(bar_DT), tick = FALSE, cex.axis = 2)
      axis(2, at = seq_len(nrow(bar_DT)) -0.5,
           labels = rev(cell_name_list), tick = FALSE, las = 1, cex.axis = 2)
      print(p)
      
      dev.off()
      
      print("done all")
     
     }
   } #corr method

} #annotation type: Specific, General

library(VennDiagram)

harmony_cell_type_list <- as.vector(unique(Seurat$harmony_FCT))

p_val_DT <- data.table()
p_val_DT[, cell_type := ""]
p_val_DT[, gene := ""]
p_val_DT[, mast_p_value := 0]
p_val_DT[, individual_1_p_value := 0]
p_val_DT[, individual_2_p_value := 0]

for (x in harmony_cell_type_list){
  
  #check conditions for DEA to have been run
  
  Idents(Seurat) <- Seurat$individual
  individual_1 <-  subset(Seurat, idents = "individual_1")
  individual_2 <-  subset(Seurat, idents = "individual_2")
  
  Idents(individual_1) <- individual_1$harmony_FCT
  sub_1 <- subset(individual_1, idents = x) 
  
  Idents(individual_2) <- individual_2$harmony_FCT
  sub_2 <- subset(individual_2, idents = x) 
  
  apoe33_1 <- subset(sub_1, APOE_Genotype %in% "APOE 33")
  apoe44_1 <- subset(sub_1, APOE_Genotype %in% "APOE 44")
  
  apoe33_2 <- subset(sub_2, APOE_Genotype %in% "APOE 33")
  apoe44_2 <- subset(sub_2, APOE_Genotype %in% "APOE 44")
  
  if ((length(colnames(apoe33_1)) >= 5) & (length(colnames(apoe44_1)) >= 5) & (length(colnames(apoe44_2)) >= 5) & (length(colnames(apoe44_2)) >= 5)){
  
    for (m in list('BH', 'BC')){
      
      #all degs, not just significant
      # wilcox_1 <- read.csv(paste("Harmony/Specific/individual_1/",x,"/wilcox_stat_dea/harmony_",x,"_significant_0.05_degs_by_",m,"_padj_value_results.csv", sep = ""))
      # wilcox_2 <- read.csv(paste("Harmony/Specific/individual_2/",x,"/wilcox_stat_dea/harmony_",x,"_significant_0.05_degs_by_",m,"_padj_value_results.csv", sep = ""))
      # 
      wilcox_1 <- read.csv(paste("Harmony/Specific/individual_1/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
      wilcox_2 <- read.csv(paste("Harmony/Specific/individual_2/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
      
      #get only those genes present in both lists
      com_degs <- c(unique(wilcox_1$gene_name), unique(wilcox_2$gene_name))
      com_degs <- as.data.table(com_degs)
      com_degs <- com_degs[duplicated(com_degs)]
      
      colnames(com_degs)<- "gene_name"
      com_degs[, statistic := 0]
      com_degs[, logFC := 0]
      com_degs[, ind_1_p_val := 0]
      com_degs[, ind_2_p_val := 0]
      com_degs[, ind_1_p_val_adj := 0]
      com_degs[, ind_2_p_val_adj := 0]
      
      if (length(com_degs != 0)){
      
      for (gene_test in com_degs$gene_name){
        
        stat_1 <- (subset(wilcox_1, gene_name == gene_test))$statistic #u-stat
        stat_2 <- (subset(wilcox_2, gene_name == gene_test))$statistic #u-stat
        logFC_1 <- (subset(wilcox_1, gene_name == gene_test))$logFC #avg_log2FC
        logFC_2 <- (subset(wilcox_2, gene_name == gene_test))$logFC #avg_log2FC
        
        #the gene should be regulated in the same direction in both individuals
        if (sign(logFC_1) == sign(logFC_2)){
          com_degs[com_degs$gene_name == gene_test]$statistic <- mean(stat_1, stat_2) #mean(stat_1, stat_2) #take the avg of the stats
          com_degs[com_degs$gene_name == gene_test]$logFC <- mean(logFC_1, logFC_2) #take the avg of the log fold change
          com_degs[com_degs$gene_name == gene_test]$ind_1_p_val <- (subset(wilcox_1, gene_name == gene_test))$p_val #ind 1 p_val
          com_degs[com_degs$gene_name == gene_test]$ind_2_p_val <- (subset(wilcox_2, gene_name == gene_test))$p_val #ind 2 p_val
          
          p_val_adj = paste(m, "_p_val_adj", sep = '')
          
          com_degs[com_degs$gene_name == gene_test]$ind_1_p_val_adj <- (subset(wilcox_1, gene_name == gene_test))$p_val_adj #ind 1 p_val_adj
          com_degs[com_degs$gene_name == gene_test]$ind_2_p_val_adj <- (subset(wilcox_2, gene_name == gene_test))$p_val_adj #ind 2 p_val_adj
        }
        else{ #else, remove the gene from the combined list
          com_degs <- com_degs[gene_name != gene_test]
        }
      }
      
      write.csv(as.data.frame(com_degs), file= paste("Harmony/Specific/combined/",x,"/harmony_",x,"_combined_wilcox_stat_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
      
      print("done csv 1")
      
      #FDR < 0.05
      res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
      
      
      #p-val (not adjusted)
      
      res_table_pval_thres <- res_tbl %>% 
        mutate(threshold = p_val < 0.05) #significant DEGs
      
      res_table_pval_thres <- res_table_pval_thres[res_table_pval_thres$p_val < 0.05,] #only look at significant DEGs based on p-value < 0.05
      res_table_pval_thres <- res_table_pval_thres[!is.na(res_table_pval_thres$p_val),] 
      res_table_pval_thres <- as.data.frame(res_table_pval_thres)
      
      
      #Venn diagram of overlapping DEGs between individual_1 & individual_2
      
      ind_1_ex <- wilcox_1$gene_name
      ind_2_ex <- wilcox_2$gene_name
      
      ind_1_ex_up <- subset(wilcox_1, logFC > 0)$gene_name
      ind_2_ex_up <- subset(wilcox_2,logFC > 0)$gene_name
  
      ind_1_ex_down <- subset(wilcox_1, logFC < 0)$gene_name
      ind_2_ex_down <- subset(wilcox_2,logFC < 0)$gene_name
      
      venn.diagram(
        x = list(ind_1_ex, ind_2_ex),
        category.names = c("Individual 1", "Individual 2"),
        print.mode=c("raw","percent"),
        lwd = 1,
        cex = 2,
        cat.cex = 2,
        main.cex = 2,
        cat.pos = 0,
        main = paste("Overlapping ",x, " DEGs Between Individuals", sep = ''),
        filename = paste('Harmony/Specific/combined/',x,'/individual_shared_DEG_venn.png', sep = ''),
        height = 5000, 
        width = 5000) 
      
      venn.diagram(
        x = list(ind_1_ex_up, ind_2_ex_up),
        category.names = c("Individual 1" , "Individual 2"),
        print.mode=c("raw","percent"),
        lwd = 1,
        cex = 2,
        cat.cex = 2,
        main.cex = 2,
        cat.pos = 0,
        main = paste("Overlapping ",x, " Upregulated DEGs Between Individuals", sep = ''),
        filename = paste('Harmony/Specific/combined/',x,'/individual_shared_up_DEG_venn.png', sep = ''),
        height = 5000, 
        width = 5000) 
      
      venn.diagram(
        x = list(ind_1_ex_down, ind_2_ex_down),
        category.names = c("Individual 1" , "Individual 2"),
        print.mode=c("raw","percent"),
        lwd = 1,
        cex = 2,
        cat.cex = 2,
        main.cex = 2,
        cat.pos = 0,
        main = paste("Overlapping ",x," Downregulated DEGs Between Individuals", sep = ''),
        filename = paste('Harmony/Specific/combined/',x,'/individual_shared_down_DEG_venn.png',  sep = ''),
        height = 5000, 
        width = 5000) 
      
      print("done venn 1")

      #Venn diagram of overlapping DEGs between combined wilcox & MAST, split by up & down regulated
      mast <- read.csv(file = paste("Harmony/Specific/All/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
      
      mast <- subset(mast, p_val_adj < 0.05)
      
      #get common genes
      meth_com_degs <- c(unique(mast$gene_name), unique(com_degs$gene_name))
      meth_com_degs <- meth_com_degs[duplicated(meth_com_degs)]
      meth_com_degs <- as.data.table(meth_com_degs)
      
      colnames(meth_com_degs)<- "X"
      meth_com_degs[, statistic := 0]
      meth_com_degs[, p_val_adj_1 := 0]
      meth_com_degs[, p_val_adj_2 := 0]
      meth_com_degs[, p_val_1 := 0]
      meth_com_degs[, p_val_2 := 0]
      meth_com_degs[, mast_p_val_adj := 0]
      meth_com_degs[, wilcox_logFC := 0] #logFC
      meth_com_degs[, mast_logFC := 0] #avg_log2FC
      
      for (gene_test in meth_com_degs$gene_name){
        
        statistic <- subset(com_degs, gene_name == gene_test)$statistic #avg stat btw ind 1 & ind 2
        p_val_adj_1 <- subset(com_degs, gene_name == gene_test)$ind_1_p_val_adj #ind 1 p_val_adj
        p_val_adj_2 <- subset(com_degs, gene_name == gene_test)$ind_2_p_val_adj #ind 2 p_val_adj
        p_val_1 <- subset(com_degs, gene_name == gene_test)$ind_1_p_val_adj #ind 1 p_val
        p_val_2 <- subset(com_degs, gene_name == gene_test)$ind_2_p_val_adj #ind 2 p_val
        mast_p_val_adj <- subset(mast, gene_name == gene_test)$p_val_adj #mast p_val_adj
        wilcox_logFC <- subset(com_degs, gene_name == gene_test)$logFC #wilcox logFC
        mast_logFC <- subset(mast, gene_name == gene_test)$avg_log2FC #mast logFC
        
        if (sign(stat_1) == sign(stat_2)){
          meth_com_degs[meth_com_degs$gene_name == gene_test]$statistic <- statistic #take the avg of the stats
          meth_com_degs[meth_com_degs$gene_name == gene_test]$p_val_adj_1 <- p_val_adj_1 #take the avg of the log fold change
          meth_com_degs[meth_com_degs$gene_name == gene_test]$p_val_adj_2 <- p_val_adj_2 #ind 1 p_val
          meth_com_degs[meth_com_degs$gene_name == gene_test]$p_val_1 <- p_val_1 #ind 1 p_val
          meth_com_degs[meth_com_degs$gene_name == gene_test]$p_val_2 <- p_val_2 #ind 1 p_val
          meth_com_degs[meth_com_degs$gene_name == gene_test]$mast_p_val_adj <- mast_p_val_adj  #ind 2 p_val
          meth_com_degs[meth_com_degs$gene_name == gene_test]$wilcox_logFC <- wilcox_logFC #wilcox logFC
          meth_com_degs[meth_com_degs$gene_name == gene_test]$mast_logFC <- mast_logFC #mast logFC
        }
        else{ #else, remove the gene from the combined list
          meth_com_degs$gene_test <- NULL
        }
        
      }
      
      print("check 1")
      
      mast_ex <- mast$gene_name
      wilcox_ex <- com_degs$gene_name
  
      
      mast_ex_up <- subset(mast, avg_log2FC > 0)$gene_name
      wilcox_ex_up <- subset(com_degs, logFC > 0)$gene_name
      
      mast_ex_down <- subset(mast, avg_log2FC < 0)$gene_name
      wilcox_ex_down <- subset(com_degs, logFC < 0)$gene_name
      
      venn.diagram(
        x = list(wilcox_ex, mast_ex),
        category.names = c("Wilcoxon" , "MAST"),
        print.mode=c("raw","percent"),
        lwd = 1,
        cex = 2,
        cat.cex = 2,
        main.cex = 2,
        cat.pos = 0,
        main = paste("Overlapping ",x, " DEGs Across DEA Methods", sep = ''),
        filename = paste('Harmony/Specific/combined/',x,'/dea_method_shared_DEG_venn.png', sep = ''),
        height = 5000, 
        width = 5000) 
      
      venn.diagram(
        x = list(wilcox_ex_up, mast_ex_up),
        category.names = c("Wilcoxon" , "MAST"),
        print.mode=c("raw","percent"),
        lwd = 1,
        cex = 2,
        cat.cex = 2,
        main.cex = 2,
        cat.pos = 0,
        main = paste("Overlapping ",x, " Upregulated DEGs Across DEA Methods", sep = ''),
        filename = paste('Harmony/Specific/combined/',x,'/dea_method_shared_up_DEG_venn.png', sep = ''),
        height = 5000, 
        width = 5000) 
      
      venn.diagram(
        x = list(wilcox_ex_down, mast_ex_down),
        category.names = c("Wilcoxon" , "MAST"),
        print.mode=c("raw","percent"),
        lwd = 1,
        cex = 2,
        cat.cex = 2,
        main.cex = 2,
        cat.pos = 0,
        main = paste("Overlapping ",x, " Downregulated DEGs Across DEA Methods", sep = ''),
        filename = paste('Harmony/Specific/combined/',x,'/dea_method_shared_down_DEG_venn.png', sep = ''),
        height = 5000, 
        width = 5000)  
      
      print("done venn 2")
      
      meth_com_degs <- meth_com_degs %>% slice_max(n= 10000000, order_by = abs(wilcox_logFC))
      
      mast_p =  meth_com_degs$mast_p_val
      ind_1_p =  meth_com_degs$p_val_1
      ind_2_p =   meth_com_degs$p_val_2
      
      if (length(meth_com_degs$gene_name) > 0){
        pdf(paste("Harmony/Specific/combined/",x,"/ind1_v_ind2_p_value_scatterplot.pdf", sep = "")) 
        
        # creating scatterplot 
        data <- data.table()
        data$individual_1_p_value <- -log10(ind_1_p)
        data$individual_2_p_value <- -log10(ind_2_p)
        data <- as.data.frame(data)
        
        p <- ggplot(data) + 
          geom_point(aes(x= individual_1_p_value, y= individual_2_p_value))+
          xlab("individual 1 -log10(p value)")+
          ylab("individual 2 -log10(p value)")+
          xlim(0, 50)+
          ylim(0, 50)
        
        print(p)
        
        dev.off()
        
        print("done scatter 1")
        
        #Barplot for top 50 DEGs in celltype: wilcox avg log fold change, MAST log fold change
        
        pdf(paste("Harmony/Specific/combined/",x,"/ind1_v_mast_p_value_scatterplot.pdf", sep = "")) 
      
        data <- data.table()
        data$individual_1_p_value <- -log10(ind_1_p)
        data$mast_p_value <- -log10(mast_p)
        data <- as.data.frame(data)
        
        p <- ggplot(data) + 
          geom_point(aes(x= individual_1_p_value, y= mast_p_value))+
          xlab("individual 1 -log10(p value)")+
          ylab("mast -log10(p value)")+
          xlim(0, 50)+
          ylim(0, 50)
        
        print(p)
        
        dev.off()
        print("done scatter 2")
        
        pdf(paste("Harmony/Specific/combined/",x,"/ind2_v_mast_p_value_scatterplot.pdf", sep = ""))
        
        data <- data.table()
        data$individual_2_p_value <- -log10(ind_2_p)
        data$mast_p_value <- -log10(mast_p)
        data <- as.data.frame(data)
        
        p <- ggplot(data) + 
          geom_point(aes(x= individual_2_p_value, y= mast_p_value))+
          xlab("individual 2 -log10(p value)")+
          ylab("mast -log10(p value)")+
          xlim(0, 50)+
          ylim(0, 50)
        
        plot(p)
        
        dev.off()
        print("done scatter 3")
        
          for (gene_test in meth_com_degs$gene_name){
              a<- subset(meth_com_degs, gene_name == gene_test)$mast_p_val
              b<- subset(meth_com_degs, gene_name == gene_test)$p_val_1
              c<- subset(meth_com_degs, gene_name == gene_test)$p_val_2
            
              new_row <- data.table("cell_type" = x, "gene" = gene_test, "mast_p_value" = -log10(a), "individual_1_p_value" = -log10(b), "individual_2_p_value" = -log10(c))
              p_val_DT <- rbindlist(list(p_val_DT, new_row))
            } #for 
          } #if
      
        } #if deg list length
      
      } #p adj method: BH or BC
  } #condition for dea
} #cell type

#p value correlation table

p_val_DT <- as.data.frame(p_val_DT)[-1,]

pdf(paste("Harmony/Specific/combined/giant_p_value_scatterplot.pdf", sep = ""))
  
  p1 <- ggplot(p_val_DT)+
        geom_point(aes(x = individual_1_p_value, y = individual_2_p_value))+ 
        facet_wrap(vars(cell_type))
  
  p2 <- ggplot(p_val_DT, aes(x = individual_1_p_value, y = mast_p_value))+
    geom_point()
  
  p2 <- p2 + facet_wrap(vars(cell_type))
  
  p3 <- ggplot(p_val_DT, aes(x = individual_2_p_value, y = mast_p_value))+
    geom_point()
  
  p3 <- p3 + facet_wrap(vars(cell_type))
   
  giant_p <- p1 + p2 + p3 + plot_layout(widths = c(50, 50, 50)) #+ plot_layout(heights = c(50, 50, 50))
  
  #dev.new(width=50, height=50)
  giant_p

dev.off()



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


pfc_tom_h <- readRDS("/output/Post_mortem_brain_validation/PFC/PFC_TOM.rds")

pfc_tom_h <- ScaleData(pfc_tom_h, features=VariableFeatures(pfc_tom_h))

pfc_tom_h <- ModuleEigengenes(
  pfc_tom_h,
  group.by.vars="projid",
  vars.to.regress = c("Study"), #Variables to regress
  scale.model.use = "linear",
  assay = "SCT",
  pc_dim = 1#First PC
)


pfc_ME <- GetMEs(pfc_tom_h, harmonized = T)

pfc_module_df <- pfc_tom_h@misc$PFC_Ast_WGCNA$wgcna_modules

pfc_module_df <- pfc_module_df %>%
  dplyr::filter(module != "grey")

#Hypergeometric Test
pathways<-fread('/projectnb/tcwlab/MSigDB/all_CPandGOs_gene_and_genesets.csv.gz')

pathways_infos<-fread('/projectnb/tcwlab/MSigDB/all_CPandGOs_genesets_metadata.csv.gz')
setnames(pathways_infos,old = 'pathway','term')

pathwaysf<-pathways[pathway.size>5&pathway.size<2000]
length(unique(pathwaysf$pathway))#12k

module_gene_df_final <- merge(pfc_module_df, pathwaysf, by.x = "gene_name", by.y = "gene")

#rm non annotated or unassigned genes
module_gene_df_final <- module_gene_df_final[!(is.na(module_gene_df_final$gene_name) | module_gene_df_final$gene_name==''), ]

res_enr<-rbindlist(lapply(split(pathwaysf, by='subcat'),function(msigdbf)OR3(split(module_gene_df_final$gene_name,module_gene_df_final$module),
                                                                             terms_list = split(msigdbf$gene,msigdbf$pathway),
                                                                             background =module_gene_df_final$gene_name)))

#add subcategory and pathway size info
res_enr_final<-merge(res_enr[padj<0.1&n.overlap>5],unique(pathways_infos,by='term'),by='term')[order(query,term,pval)]


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

outdir <- "/output/Post_mortem_brain_validation/PFC/Module_GO/"

for (i in 1:length(unique(res_enr_copy$query))) {
  df <- res_enr_copy[res_enr_copy$query. == unique(res_enr_copy$query)[i], ]
  
  df_filtered <- df[df$padj <= 0.1, ]
  
  df_sorted <- df_filtered[order(abs(df_filtered$NES), decreasing = T), ]
  
  df_sorted <- df_sorted[!duplicated(df_sorted$pathway), ]
  
  if (nrow(df_sorted) > 0) {
    
    #Top 50 pathways
    n.pathway <- min(nrow(df_sorted), 50)
    
    emmaplot(df_sorted[1:n.pathway, ])
    ggsave(paste0(outdir, unique(res_enr_copy$query)[i], "_Emma.pdf"), height = 10, width = 20)
  }
}


pfc_tom$Study <- ifelse(pfc_tom$Study == "ROS", 0, 1)
#pfc_tom$APOE <- ifelse(pfc_tom$APOE == 33, 0, ifelse(pfc_tom$APOE == 34, NA, 1))
pfc_tom$APOE <- ifelse(pfc_tom$APOE == 33, 0, 1)

traits <- c("Study", "msex", "APOE", "braaksc", "ceradsc", "cogdx", "dcfdx_lv")

res_enr_copy <- read.csv("/output/Post_mortem_brain_validation/PFC/res_enr.csv")

pfc_tom <- ModuleTraitCorrelation(
  pfc_tom,
  traits = traits,
  group.by = NULL,
  features = 'hMEs', 
  cor_method = 'spearman',
)

trait_res <- GetModuleTraitCorrelation(pfc_tom)

trait_df <- trait_res$cor$all_cells
trait_df <- trait_df[, colnames(trait_df) %in% res_enr_copy$query.]

trait_p_df <- trait_res$fdr$all_cells
trait_p_df <- trait_p_df[, colnames(trait_p_df) %in% res_enr_copy$query.]

plot_trait_heatmap <- function(cor, p_val, filename) {
  
  text_p <- paste(signif(cor, 2), "\n(", round(signif(p_val, 2), 3), ")", sep="")
  text_value <- signif(cor, 2)
  dim(text_value) <- dim(cor)
  
  outdir <- "/output/Post_mortem_brain_validation/PFC/"
  
  # Create figure with * instead of p-values and no correlation values:
  module_trait_p_mark <- sapply(as.data.frame(p_val), function(col) {
    sapply(col, function(x) {
      if (x < 0.0001 & x > 0) {
        return("\u2731\u2731\u2731")
      } else if (x < 0.001 & x >= 0.0001) {
        return("\u2731\u2731")
      } else if (x < 0.01 & x >= 0.001) {
        return("\u2731")
      } else if (x == 0) {
        return("")
      } else {
        return(" ")
      }
    })
  })
  
  color=colorRampPalette(c("dodgerblue", "white", "firebrick1"))(50)
  
  cor_t <- t(cor)
  module_trait_p_mark_t <- t(module_trait_p_mark)
  
  print(pheatmap::pheatmap(cor_t, cluster_cols = FALSE, cluster_rows = TRUE, 
                           color = color, angle_col = "45", main = filename, 
                           fontsize_row = 7, fontsize_col = 10, na_col = "white", 
                           cellwidth = 25, display_numbers = as.matrix(module_trait_p_mark_t), 
                           number_color = "black", breaks = seq(-0.2, 0.2, length.out = 51)))
}

plot_trait_heatmap(trait_df, trait_p_df, "PFC 44_vs_33(+34) Module Trait Analysis")


res_enr <- read.csv( "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Post_mortem_brain_validation/PFC_427/PFC_427_scdemon_enr_res.csv")
ME_of_interest <- unique(res_enr$query.)
ME_data <- eigengene_df[paste0("X",ME_of_interest)]

metadata <- read.csv("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Post_mortem_brain_validation/PFC_427/Metadata.csv")

ME_data$Genotype <- metadata$apoe_genotype[match(rownames(ME_data), metadata$read)]
ME_data$Genotype <- ifelse(ME_data$Genotype == 33, "APOE 33", "APOE 44")

ME_data_Genotype <- ME_data %>%
  pivot_longer(cols = -Genotype, # Select columns to pivot
               names_to = "Module",     # New column for former column names
               values_to = "Eigengene")  

ggplot(ME_data_Genotype, aes(x = Genotype, y = Eigengene, fill = Genotype)) +
  geom_boxplot() +
  facet_wrap(~Module, scales = "free_x", strip.position = "top") +
  stat_compare_means(
      aes(label = paste0("wilcox p = ", ..p.format..)), 
      method = "wilcox.test",
      label.y = max(ME_data_Genotype$Count) + 10 # Adjust position of labels
    ) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.text = element_text(size = 14)) +
  labs(title = "Post-mortem PFC Module Eigengene Expression", x = "", y = "Expression Level")


# Use TFBSTools to get the motif position weight matrices (JASPAR 2020 database)
pfm_core <- TFBSTools::getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# run the motif scan
seurat_obj <- MotifScan(
  seurat_obj,
  species_genome = 'hg38',
  pfm = pfm_core,
  EnsDb = EnsDb.Hsapiens.v86
)

# get the motif df:
motif_df <- GetMotifs(seurat_obj)

# keep all TFs, and then remove all genes from the grey module
tf_genes <- unique(motif_df$gene_name)
modules <- GetModules(seurat_obj)
nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name
genes_use <- c(tf_genes, nongrey_genes)

# update the gene list and re-run SetDatExpr
seurat_obj <- SetWGCNAGenes(seurat_obj, genes_use)
seurat_obj <- SetDatExpr(seurat_obj, group_name = "Ast", group.by= "major_cell_type", assay = 'SCT', slot = 'data')

saveRDS(seurat_obj, "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Post_mortem_brain_validation/PFC_427/PFC_TOM.rds")

# define model params:
model_params <- list(
    objective = 'reg:squarederror',
    max_depth = 1,
    eta = 0.1,
    nthread=16,
    alpha=0.5
)

# construct the TF network
seurat_obj <- ConstructTFNetwork(seurat_obj, model_params)

results <- GetTFNetwork(seurat_obj)
head(results)

#Top TFs for each gene.
# seurat_obj <- AssignTFRegulons(
#     seurat_obj,
#     strategy = "A",
#     reg_thresh = 0.01,
#     n_tfs = 10
# )

#Top genes for each TF
# seurat_obj <- AssignTFRegulons(
#     seurat_obj,
#     strategy = "B",
#     reg_thresh = 0.01,
#     n_genes = 50
# )

#All TF-gene pair above reg_thresh
seurat_obj <- AssignTFRegulons(
    seurat_obj,
    strategy = "C",
    reg_thresh = 0.1
)

#Plot gene targets of selected TF
RegulonBarPlot(seurat_obj, selected_tf='ZNF263')


# positive regulons
seurat_obj <- RegulonScores(
    seurat_obj,
    target_type = 'positive',
    ncores=8
)

# negative regulons
seurat_obj <- RegulonScores(
    seurat_obj,
    target_type = 'negative',
    cor_thresh = -0.05,
    ncores=8
)

# access the results:
pos_regulon_scores <- GetRegulonScores(seurat_obj, target_type='positive')
neg_regulon_scores <- GetRegulonScores(seurat_obj, target_type='negative')

seurat_obj@meta.data
Idents(seurat_obj) <- "apoe_genotype"
VlnPlot(seurat_obj, features = "pos_regulon_scores", pt.size = 0) + NoLegend() +
  labs(x = "Cluster") +
  ggtitle("PFC nFeature_RNA")


cur_tfs <- c('MAF', 'CTCF', 'ZNF263')

#Plot TF Network
TFNetworkPlot(
    seurat_obj, selected_tfs=cur_tfs, 
    target_type = 'both', 
    label_TFs=0, depth=2
) + ggtitle("pos & neg targets")

#Plot Correlation among TFs
ModuleRegulatoryHeatmap(
  seurat_obj, feature='delta', dendrogram=FALSE
) + ggtitle('TFs only')
