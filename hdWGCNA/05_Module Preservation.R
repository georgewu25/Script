library(Seurat)
library(WGCNA)

#######iPSC Count Matrix
#Transpose Count matrix so genes as rows
ipsc_mtx <- read.csv("/output/WGCNA/normalized_expression.csv")

#Get a deseq data frame for gene name annotation
deseq_df <- read.csv("/results/DESEQ2/APOE33_d24_vs_d0/APOE33_d24_vs_d0/Results_GTFAnnotated_NoGeneIDDuplicates.csv")

length(intersect(deseq_df$X, colnames(ipsc_mtx, -1))) #10574 out of 10576 genes are annotated
length(colnames(ipsc_mtx, -1))

#Annotate gene names
gene_names <- deseq_df$gene_name[match(colnames(ipsc_mtx)[-1], deseq_df$X)]
#Use gene ID if gene name is absent
colnames(ipsc_mtx)[-1] <- ifelse(is.na(gene_names), colnames(ipsc_mtx)[-1], gene_names)
rownames(ipsc_mtx) <- ipsc_mtx$X
ipsc_mtx <- ipsc_mtx[,-1]

#######Post_mortem Brain
pfc_tom <- readRDS("/output/Post_mortem_brain_validation/PFC/PFC_TOM.rds")
pfc_mtx <- pfc_tom[["SCT"]]$data
dim(pfc_mtx) #25625 X 20021

#Transpose count matrix, genes as columns
pfc_mtx_t <- t(pfc_mtx)

#######Module Preservation
#Filter common genes
common_genes <- intersect(colnames(pfc_mtx_t), colnames(ipsc_mtx))
length(common_genes) #9722

#iPSC filter
ipsc_mtx_subset <- ipsc_mtx[, colnames(ipsc_mtx) %in% common_genes]
dim(ipsc_mtx_subset) #32 X 9722

#Brain filter
pfc_mtx_t_subset <- pfc_mtx_t[, colnames(pfc_mtx_t) %in% common_genes]
dim(pfc_mtx_t_subset) # 20021 X 9722

######iPSC module colorlist
ipsc_tom <- readRDS("/output/WGCNA/deepsplit2/TOM_deepsplit2.rds")

ipsc_module_color <- labels2colors(ipsc_tom$colors, zeroIsGrey=TRUE)
names(ipsc_module_color)<-names(ipsc_tom$colors)

#Change colorlist names to gene names
names(ipsc_module_color) <- deseq_df$gene_name[match(names(ipsc_module_color), deseq_df$X)]
#Colorlist Filter
ipsc_module_color <- ipsc_module_color[names(ipsc_module_color) %in% common_genes]

######Brain module colorlist
pfc_module_gene_df <- read.csv("/output/Post_mortem_brain_validation/PFC/Module_gene.csv")

pfc_module_color <- setNames(as.vector(pfc_module_gene_df$module), pfc_module_gene_df$gene_name)
pfc_module_color <- pfc_module_color[names(pfc_module_color) %in% common_genes]

#######WGCNA::ModulePreservation
#iPSC Preservation in PFC
ipsc_multiExpr <- list(ipsc = list(data = as.matrix(ipsc_mtx_subset)), brain = list(data = as.matrix(pfc_mtx_t_subset)))

ipsc_colorList <- list(ipsc = ipsc_module_color)

#QC
ipsc_qc<-goodSamplesGenes(ipsc_multiExpr$ipsc$data)
ipsc_qc$allOK#TRUE

brain_qc<-goodSamplesGenes(ipsc_multiExpr$brain$data)
brain_qc$allOK#False

#Remove bad samples and genes from brain dataset
goodSamples <- which(brain_qc$goodSamples)
goodGenes <- which(brain_qc$goodGenes)

ipsc_multiExpr$brain$data <- ipsc_multiExpr$brain$data[goodSamples, goodGenes]

# Calculate module preservation statistics
ipsc_preservation <- WGCNA::modulePreservation(ipsc_multiExpr, ipsc_colorList, 
                                               networkType = "signed",
                                               corFnc='bicor',
                                               quickCor = 0,
                                               greyName = "grey", 
                                               goldName = NULL,
                                               referenceNetworks = 1, 
                                               nPermutations = 50, 
                                               parallelCalculation=FALSE,
                                               verbose = 3
)

saveRDS(ipsc_preservation, "/output/Post_mortem_brain_validation/PFC/ipsc_preservation_in_PFC.rds")

mp = readRDS("/output/Post_mortem_brain_validation/PFC/ipsc_preservation_in_PFC.rds")

res_enr_copy <- read.csv("/output/Post_mortem_brain_validation/PFC/res_enr.csv")

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
            signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#preservation of rna modules in proteins
# Module labels and module sizes are also contained in the results
pfc_preservation_res <- mp$preservation$observed[[ref]][[test]]
pfc_preservation_Z <- mp$preservation$Z[[ref]][[test]]

modColors = rownames(pfc_preservation_res)
moduleSizes = pfc_preservation_Z$moduleSize

# leave grey and gold modules out
plotMods = modColors %in% unique(res_enr_copy$query.)

# Text labels for points
text = modColors[plotMods]

# Auxiliary convenience variable
plotData = cbind(pfc_preservation_res$medianRank.pres, pfc_preservation_Z$Zsummary.pres)

# Start the plot
#sizeGrWindow(10, 5);
#pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))

for (p in 1:2){
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  #labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}

plot(moduleSizes[plotMods], plotData[plotMods, 2], col = 1, bg = modColors[plotMods], pch = 21,
     main = "Preservation Zsummary in PFC",
     cex = 2.4,
     ylab = mains[2], xlab = "Module size", log = "x",
     ylim = ylim,
     xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)