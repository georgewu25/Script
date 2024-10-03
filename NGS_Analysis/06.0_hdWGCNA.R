library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(devtools)
library(WGCNA)
library(hdWGCNA)


seurat_obj_normalized <- readRDS("/output/Post_mortem_brain_validation/PFC_processed_Seurat.rds")

hdwgcna_obj <- SetupForWGCNA(
  seurat_obj_normalized,
  gene_select = "fraction",
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "WB_Ast_WGCNA"
)

#Constructing MetaCells
hdwgcna_obj <- MetacellsByGroups(
  hdwgcna_obj,
  group.by = c("major_cell_type", "individualID"),
  reduction = 'pca',
  k = 25,
  max_shared = 10,
  ident.group = 'major_cell_type'
)

hdwgcna_obj <- NormalizeMetacells(hdwgcna_obj)

#Using full Seurat object since Metacell is not performed due to using subset and processed data
hdwgcna_obj_final <- SetDatExpr(
  hdwgcna_obj,
  group_name = "Ast",
  group.by= "major_cell_type",
  assay = 'RNA',
  slot = 'data'
)

hdwgcna_obj_final <- TestSoftPowers(
  hdwgcna_obj_final,
  networkType = 'signed'
)
#Warning: executing %dopar% sequentially: no parallel backend registered
#Warning: bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
#Warning: bicor: zero MAD in variable 'y'. Pearson correlation was used for individual columns with zero (or missing) MAD.


plot_list <- PlotSoftPowers(hdwgcna_obj_final)
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(hdwgcna_obj_final)

#Construct co-expression network:
hdwgcna_obj_tom <- ConstructNetwork(
  hdwgcna_obj_final,
  tom_name = 'PFC_Ast',
  overwrite_tom = TRUE
)

saveRDS(hdwgcna_obj_tom, "/output/Post_mortem_brain_validation/PFC_tom.rds")


PlotDendrogram(hdwgcna_obj)

hdwgcna_obj_tom <- GetTOM(hdwgcna_obj)

metadata <- hdwgcna_obj@meta.data
count_mtx <- hdwgcna_obj@assays$RNA$data

module_gene_df <- hdwgcna_obj@misc$PFC_Ast_WGCNA$wgcna_modules
connectivity <- hdwgcna_obj@misc$PFC_Ast_WGCNA$wgcna_powerTable
eigengene_df <- hdwgcna_obj@misc$PFC_Ast_WGCNA$datExpr

#Transpose eigengene_df
eigengene_df_t <- as.data.frame(t(eigengene_df))
eigengene_df_t$gene_name <- rownames(eigengene_df_t)

module_gene_df_final <- merge(eigengene_df_t, module_gene_df[, c(1,2)], by = "gene_name")
#Remove genes in grey module
module_gene_df_final <- module_gene_df_final %>%
  dplyr::filter(module != "grey")

#Hypergeometric Test
pathways<-fread('/projectnb/tcwlab/MSigDB/all_CPandGOs_gene_and_genesets.csv.gz')

pathways_infos<-fread('/projectnb/tcwlab/MSigDB/all_CPandGOs_genesets_metadata.csv.gz')
setnames(pathways_infos,old = 'pathway','term')

pathwaysf<-pathways[pathway.size>5&pathway.size<2000]
length(unique(pathwaysf$pathway))#12k

#rm non annotated or unassigned genes
module_gene_df_final <- module_gene_df_final[!(is.na(module_gene_df_final$gene_name) | module_gene_df_final$gene_name==''), ]

res_enr<-rbindlist(lapply(split(pathwaysf, by='subcat'),function(msigdbf)OR3(split(module_gene_df_final$gene_name,module_gene_df_final$module),
                                                                             terms_list = split(msigdbf$gene,msigdbf$pathway),
                                                                             background =module_gene_df_final$gene_name)))

#add subcategory and pathway size info
res_enr_final<-merge(res_enr[padj<0.1&n.overlap>5],unique(pathways_infos,by='term'),by='term')[order(query,term,pval)]


write.csv(res_enr_final, "/output/hdWGCNA/res_enr.csv", row.names = FALSE, quote = FALSE)
