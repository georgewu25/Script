library(Seurat)

pfc_tom <- readRDS("/output/Post_mortem_brain_validation/PFC_427/PFC_TOM.rds")

#########Plot Metadata
plot_df <- pfc_tom@meta.data

plot_df$apoe_genotype <- factor(plot_df$apoe_genotype, levels = c(33, 34, 44))

ggplot(plot_df, aes(x = apoe_genotype)) +
  geom_bar() +
  theme_minimal() +
  labs(x = "APOE Genotyoe", y = "Cell Count", title = "PFC MetaCell Count Distribution")

pfc_tom <- ScaleData(pfc_tom, features=VariableFeatures(pfc_tom))

#########Compute module eigengenes
pfc_tom <- ModuleEigengenes(
  pfc_tom,
  group.by.vars="projid",
  vars.to.regress = c("Study"), #Variables to regress
  scale.model.use = "linear",
  assay = "SCT",
  pc_dim = 1#First PC
)

pfc_ME <- GetMEs(pfc_tom, harmonized = T)


########Module Trait Analysis
#Metadata processing
pfc_tom$Study <- ifelse(pfc_tom$Study == "ROS", 0, 1)
pfc_tom$apoe_genotype <- ifelse(pfc_tom$apoe_genotype == 33, 0, ifelse(pfc_tom$apoe_genotype == 34, NA, 1))
#pfc_tom$apoe_genotype <- ifelse(pfc_tom$apoe_genotype == 33, 0, 1)

traits <- c("Study", "msex", "apoe_genotype", "pmi", "braaksc", "ceradsc", "cogdx", "dcfdx_lv")

#Use hub module genes
pfc_tom <- ModuleTraitCorrelation(
  pfc_tom,
  traits = traits,
  group.by = NULL,
  features = 'hMEs', 
  cor_method = 'spearman',
)

#Use module genes
pfc_tom <- ModuleTraitCorrelation(
  pfc_tom,
  traits = traits,
  group.by = NULL,
  features = 'MEs', 
  cor_method = 'spearman',
)

trait_res <- GetModuleTraitCorrelation(pfc_tom)

plot_trait_heatmap <- function(eigen_gene, trait_data, expression_data, filename) {
  
  # Correlate modules with traits using Spearman's correlation function:
  module_trait_cor <- trait_res$cor$all_cells
  
  module_trait_p <- trait_res$fdr$all_cells
  
  text_p <- paste(signif(module_trait_cor, 2), "\n(", round(signif(module_trait_p, 2), 3), ")", sep="")
  text_value <- signif(module_trait_cor, 2)
  dim(text_value) <- dim(module_trait_cor)
  
  outdir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Post_mortem_brain_validation/HP/"
  
  # Create figure with * instead of p-values and no correlation values:
  # p<0.001 '***', p<0.01 '**', p<0.05 '*', p<0.1 '.' 
  module_trait_p_mark <- module_trait_p
  module_trait_p_mark[module_trait_p_mark < 0.001] <- "***"
  module_trait_p_mark[module_trait_p_mark < 0.01 & module_trait_p_mark >= 0.001] <- "**"
  module_trait_p_mark[module_trait_p_mark < 0.1 & module_trait_p_mark >= 0.01] <- "*"
  module_trait_p_mark[module_trait_p_mark > 0.1] <- ""
  
  module_trait_p_mark <- sapply(as.data.frame(module_trait_p), function(col) {
    sapply(col, function(x) {
      if (x < 0.001) {
        return("\u2731\u2731\u2731")
      } else if (x < 0.01 & x >= 0.001) {
        return("\u2731\u2731")
      } else if (x < 0.1 & x >= 0.01) {
        return("\u2731")
      } else {
        return(" ")
      }
    })
  })
  
  color=colorRampPalette(c("dodgerblue", "white", "firebrick1"))(50)
  
  print(pheatmap::pheatmap(module_trait_cor, cluster_cols = FALSE, cluster_rows = T, color =  color, angle_col = "45", main = filename, fontsize_row = 7, fontsize_col = 10, na_col = "white", cellwidth = 25, display_numbers = as.matrix(module_trait_p_mark), number_color = "black", breaks = seq(-0.2, 0.2, length.out = 51)))
}


eigen_gene <- pfc_ME

expression_data <- t(pfc_tom@assays$SCT@data)

traits <- c("msex", "apoe_genotype", "pmi", "braaksc", "ceradsc", "cogdx", "dcfdx_lv")

trait_data <- pfc_tom@meta.data %>%
  dplyr::select(all_of(traits))

plot_trait_heatmap(eigen_gene, expression_data, trait_data, "44_vs_33 Hub Module Gene Trait Analysis")