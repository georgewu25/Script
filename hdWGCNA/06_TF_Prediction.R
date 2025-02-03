library(JASPAR2020)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(xgboost)
library(hdWGCNA)


seurat_obj <- readRDS("/output/Post_mortem_brain_validation/PFC_427/PFC_TOM.rds")

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