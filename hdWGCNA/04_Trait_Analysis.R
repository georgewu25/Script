library(Seurat)
library(tidyverse)
library(cowplot)
#library(patchwork)
library(WGCNA)
library(hdWGCNA)

pfc_tom <- readRDS("/output/Post_mortem_brain_validation/PFC/PFC_TOM.rds")

length(unique(pfc_tom@meta.data$projid)) #423 Individuals

#########Plot Metadata
plot_df <- pfc_tom@meta.data

plot_df$apoe_genotype <- factor(plot_df$apoe_genotype, levels = c(33, 34, 44))

ggplot(plot_df, aes(x = apoe_genotype)) +
  geom_bar() +
  theme_minimal() +
  labs(x = "APOE Genotyoe", y = "Cell Count", title = "PFC MetaCell Count Distribution")

pfc_ME <- GetMEs(pfc_tom, harmonized = T)

########Module Trait Analysis
#Metadata processing
colnames(pfc_tom@meta.data)[13] <- "APOE"

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