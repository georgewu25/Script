#Demultiplexing

#Load Libraries 
library(cowplot)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(patchwork)
library(RColorBrewer)
library(rhdf5)
library(Seurat)
library(scales)
library(utils)

#############################################
#Redone: Cite Seq Count

#HTO Matrix
cite_seq_out <- paste(source, "Sample_TCW-Hash-01-30-2020-P1-HTO/compbio/citeseq_outs/umi_count/", sep = "")
HTO_matrix <- Read10X(data.dir = cite_seq_out, gene.column = 1)

#UMI Matrix: filtered by CellRanger
cell_ranger_out <- paste(source, "Sample_TCW-Hash-01-30-2020-P1-cDNA/compbio/cellranger/outs/filtered_feature_bc_matrix/", sep = '')
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
cell_ranger_out <- paste(source, "Sample_TCW-Hash-01-30-2020-P1-cDNA/compbio/cellranger/outs/filtered_feature_bc_matrix/", sep = '')
UMI_matrix <- Read10X(data.dir = cell_ranger_out, gene.column = 1)
head(UMI_matrix)
length(colnames(UMI_matrix)) #30697

#HTO Matrix
cite_seq_out <- paste(source, "Sample_TCW-Hash-01-30-2020-P1-HTO/compbio/citeseq_outs/umi_count/", sep = "")
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
