library(Seurat)
library(WGCNA)

#######iPSC
#Transpose Count matrix so genes as rows
ipsc_mtx <- read.csv("/output/WGCNA/normalized_expression.csv")

#Get a deseq data frame fro gene name annotation
deseq_df <- read.csv("/results/DESEQ2/APOE33_d24_vs_d0/APOE33_d24_vs_d0/Results_GTFAnnotated_NoGeneIDDuplicates.csv")

length(intersect(deseq_df$X, colnames(ipsc_mtx, -1))) #10574 out of 10576 genes are annotated
length(colnames(ipsc_mtx, -1))

#Annotation
gene_names <- deseq_df$gene_name[match(colnames(ipsc_mtx)[-1], deseq_df$X)]

colnames(ipsc_mtx)[-1] <- ifelse(is.na(gene_names), colnames(ipsc_mtx)[-1], gene_names)

ipsc_mtx_t <- t(as.matrix(ipsc_mtx[,- 1]))
colnames(ipsc_mtx_t) <- ipsc_mtx[,1]
dim(ipsc_mtx_t) #10575 X 32

ipsc_tom <- readRDS("/output/WGCNA/deepsplit2/TOM_deepsplit2.rds")

ipsc_module_color <- labels2colors(ipsc_tom$colors, zeroIsGrey=TRUE)
names(ipsc_module_color)<-names(ipsc_tom$colors)


#######Post_mortem Brain Processed
seurat_obj_normalized <- readRDS("/output/Post_mortem_brain_validation/PFC_processed_Seurat.rds")

PFC_mtx <- seurat_obj_normalized@assays$RNA$counts
dim(PFC_mtx) #33538 X 11434

#######Module Preservation
#Keep only common genes
common_genes <- intersect(rownames(PFC_mtx), rownames(ipsc_mtx_t))

ipsc_mtx_t_subset <- ipsc_mtx_t[rownames(ipsc_mtx_t) %in% common_genes,]
dim(ipsc_mtx_t_subset) #9973 X 32

PFC_mtx_subset <- PFC_mtx[rownames(PFC_mtx) %in% common_genes,]
dim(PFC_mtx_subset) #9973 X 11434

#Change colorlist
names(ipsc_module_color) <- deseq_df$gene_name[match(names(ipsc_module_color), deseq_df$X)]

#Module Preservation
multiExpr <- list(ipsc = list(data = ipsc_mtx_t_subset), brain = list(data = PFC_mtx_subset))

multiExpr <- list(ipsc = ipsc_mtx_t_subset, brain = PFC_mtx_subset)
#Convert to same data type
multiExpr$brain$data <- as.matrix(multiExpr$brain$data)

colorList <- list(ipsc = ipsc_module_color[names(ipsc_module_color) %in% common_genes])

length(colorList[[1]]) #9973

# Calculate module preservation statistics
preservation <- WGCNA::modulePreservation(multiExpr, colorList, 
                                          dataIsExpr = TRUE,
                                          networkType = "signed",
                                          checkData = TRUE,
                                          greyName = "grey", goldName = NULL,
                                          referenceNetworks = 1, 
                                          nPermutations = 1, 
                                          verbose = 3)


saveRDS(preservation, "output/hdWGCNA/preservation.rds")