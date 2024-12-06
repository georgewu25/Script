library(Seurat)
library(tidyverse)
library(cowplot)
#library(patchwork)
library(WGCNA)
library(hdWGCNA)

source("/scripts/OR3.R")
source("/scripts/Emmaplot.R")

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

res_enr_copy <- read.csv("/output/Post_mortem_brain_validation/PFC/res_enr.csv")

unique(res_enr_copy$query.)

colors <- c("black", "blue", "brown", "green", "magenta", "pink", "red", "turquoise", "yellow")

plot_df <- list()

for (i in 1:length(unique(res_enr_copy$query.))) {
  
  color <- colors[i]
  
  #Top 3 most significant unique pathways in each module
  df_subset <- res_enr_copy %>%
    dplyr::filter(query. == color) %>%
    arrange(padj) %>%
    slice_head(n=3)
  
  plot_df[[color]] <- df_subset
}

plot_df <- rbindlist(plot_df)

which(duplicated(plot_df$term))

plot_df$term[27] <- "GOCC_NEURON_PROJECTION"

plot_df$term <- factor(plot_df$term, levels = plot_df$term)

ggplot(plot_df, aes(x=-log10(padj), y=term, fill= query.))+
  geom_col() +
  scale_fill_manual(values = c("black"="black", "blue"="blue", "brown"="brown", "green"="green", "magenta"="magenta", "pink"="pink", "red"="red", "turquoise"="turquoise", "yellow"="yellow")) +
  geom_vline(xintercept = -log10(0.1), linetype = "solid", color = "red") + 
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 30),
    axis.title.x = element_text(size = 20),
    panel.margin = unit(.05, "lines"),
    panel.border = element_blank(),
    strip.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2)
  ) +
  guides(fill = FALSE) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 76)) +
  scale_y_discrete(limits = rev) +
  labs(title = "", x = "-Log10(FDR)", y = "", fill = "")
