plot_trait_heatmap <- function(eigen_gene, trait_data, expression_data, filename) {
  
  # Correlate modules with traits using Spearman's correlation function:
  module_trait_cor <- cor(eigen_gene, trait_data, use="p", method="spearman") # Obtain correlation value
  
  module_trait_p <- corPvalueStudent(module_trait_cor, nrow(expression_data)) # Obtain Student asymptotic p-values. 
  
  text_p <- paste(signif(module_trait_cor, 2), "\n(", round(signif(module_trait_p, 2), 3), ")", sep="")
  text_value <- signif(module_trait_cor, 2)
  dim(text_value) <- dim(module_trait_cor)
  
  outdir <- "/output/WGCNA/deepsplit2/trait_analysis/"
  
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
  
  colnames(module_trait_cor) <- c("0hr", "8hr", "24hr/D0hr", "D24hr")
  
  pheatmap::pheatmap(module_trait_cor, cluster_cols = FALSE, cluster_rows = T, color =  color, angle_col = "45", main = "", fontsize_row = 7, fontsize_col = 10, na_col = "white", cellwidth = 25, display_numbers = as.matrix(module_trait_p_mark), number_color = "black", breaks = seq(-1, 1, length.out = 51))
}

outdir <- "/output/WGCNA/deepsplit2/"

expression_df <- read.csv("/output/WGCNA/normalized_expression.csv")

trait_df <- read.csv("/Metadata.csv", header=TRUE)

# At 24 hr, it is assumed to have no more uptake
trait_df$Treatment[17:32] <- "N"

colnames(trait_df)[14] <- "Uptake"
colnames(trait_df)[15] <- "Degradation"

# Reformat the sample name
trait_df <- trait_df[, c(1,2,3,5,6,7,10,11,14,15)]

#Numberical variables
num_vars<-colnames(trait_df)[c(2,6,8,10,11)]

#Categorical variables
cat_vars<-colnames(trait_df)[c(3,4,7,9)]

# Convert categorical variables to dummy variables
for (i in cat_vars) {
  trait_df[, colnames(trait_df) == i] <- as.numeric(factor(trait_df[, colnames(trait_df) == i])) -1
}

#Convert individuals to dummy variables due to four levels
dummy <- model.matrix(~ as.factor(trait_df$Individual) - 1)
colnames(dummy) <- levels(as.factor(trait_df$Individual))

# Remove original individual column from the dataframe
trait_df <- trait_df[, -5]

# Combine the original dataframe with the dummy variables dataframe
trait_data_final <- cbind(trait_df, dummy)

# Change sample names to match
trait_data_final$Sample_Name <- gsub("-", ".", trait_data_final$Sample_Name)

# Match the samples
rows_matched <- match(expression_df[, 1], trait_data_final$Sample_Name)

trait_data_final <- trait_data_final[rows_matched, ]

rownames(trait_data_final) <- trait_data_final$Sample_Name

trait_data_final <- trait_data_final[, -1]

#Correlation with all metadata
plot_trait_heatmap(ME, trait_data_final, expression_df, "Trait_analysis")

#Only modules whose genes significantly overlap with pathway references (See Functional Enrichment Analysis) 
plot_trait_heatmap(ME[paste0("ME", unique(res_enr_final$query.))], trait_data_final, expression_df, "Trait_analysis_Final")

rownames(ME) <- expression_df[,1]

#Interaction between genotype and treatment
trait_data_copy <- trait_df
trait_data_copy$APOE.44_vs_33_0h <- ifelse(trait_data_copy$Degradation == 0, ifelse(trait_data_copy$APOE == 0, 0, 1), NA)
trait_data_copy$APOE.44_vs_33_8h <- ifelse(trait_data_copy$Degradation == 8, ifelse(trait_data_copy$APOE == 0, 0, 1), NA)
trait_data_copy$APOE.44_vs_33_24h <- ifelse(trait_data_copy$Degradation == 24, ifelse(trait_data_copy$APOE == 0, 0, 1), NA)
trait_data_copy$APOE.44_vs_33_48h <- ifelse(trait_data_copy$Degradation == 48, ifelse(trait_data_copy$APOE == 0, 0, 1), NA)

trait_df2 <- trait_data_copy[,c(10:13)]

#Only modules whose genes significantly overlap with pathway references (See Functional Enrichment Analysis) 
pdf("/output/manuscript_data/figure/Trait_Analysis_Heatmap.pdf", height = 10, width = 10)
print(plot_trait_heatmap(ME[paste0("ME", unique(res_enr_final$query.))], trait_df2, expression_df, "APOE 44 vs APOE 33_Module_Trait_Relationship"))
dev.off()