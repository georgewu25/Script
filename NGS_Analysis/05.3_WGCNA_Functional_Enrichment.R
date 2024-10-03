#Filtered input expression data from Feature Count
df_filtered <- read.csv("/output/WGCNA/normalized_expression.csv")

outdir <- "/output/WGCNA/deepsplit2/"

#Transform the expression data
df_filtered_t <- as.data.frame(t(df_filtered[, -1]))
colnames(df_filtered_t) <- df_filtered[, 1]
df_filtered_t$gene_name <- rownames(df_filtered_t)

library(data.table)

# Save gene id and module color infromation from merged module into a table
modules_df<-data.table(gene_id=names(module_color),module=module_color)

# Merge the ME group with the expression data on gene_id
modules_df<-merge(modules_df, df_filtered_t, by.x = 'gene_id', by.y = "gene_name")

# Order by ME group
modules_df<-modules_df[order(module)]

#Convert to gene ID to EnsembleID
#Save the gene ids to csv for online conversion using Biotools
write.csv(modules_df$gene_id, paste0(outdir, "genes_to_convert.csv"), row.names = FALSE, quote = FALSE)


gene_symbol <- read.csv(paste0(outdir, "gene_names_converted.csv"))

#Keep the gene ids if conversion to gene symbols failed
empty_rows <- gene_symbol$gene_symbol == ""
gene_symbol$gene_symbol[empty_rows] <- gene_symbol$gene_id[empty_rows]

modules_df$gene_symbol <- gene_symbol$gene_symbol

write.csv(modules_df, paste0(outdir, "module_genes.csv"), row.names = FALSE, quote = FALSE)

source("/scripts/OR3.R")

outdir <- "/output/WGCNA/deepsplit2/"

#Gene-module
module_gene_df <- read.csv(paste0(outdir, "module_genes.csv"))

library(data.table)

pathways<-fread('/MSigDB/all_CPandGOs_gene_and_genesets.csv.gz')

pathways_infos<-fread('/MSigDB/all_CPandGOs_genesets_metadata.csv.gz')
setnames(pathways_infos,old = 'pathway','term')

#test pathways <2k et >5k
pathwaysf<-pathways[pathway.size>5&pathway.size<2000]
length(unique(pathwaysf$pathway))#12k

#rm non annotated or unassigned genes
module_gene_df<-module_gene_df[!(is.na(module_gene_df$gene_symbol)|module_gene_df$gene_symbol==''), ]

res_enr<-rbindlist(lapply(split(pathwaysf, by='subcat'),function(msigdbf)OR3(split(module_gene_df$gene_symbol,module_gene_df$module),
                                                                             terms_list = split(msigdbf$gene,msigdbf$pathway),
                                                                             background =module_gene_df$gene_symbol)))

#add subcategory and pathway size info
res_enr_final<-merge(res_enr[padj<0.1&n.overlap>5],unique(pathways_infos,by='term'),by='term')[order(query,term,pval)]


write.csv(res_enr_final,paste0(outdir, 'functional_enrichment_padj0.1_overlap5.csv'), row.names = FALSE, quote = FALSE)

res_enr_final <- read.csv(paste0(outdir, 'functional_enrichment_padj0.1_overlap5.csv'))

#Top 3 Postive and Negative Enrichment for each module
for (i in unique(res_enr_final$query.)) {
  
  df_subset <- res_enr_final[res_enr_final$query. == i, ]
  df_subset_sorted <- df_subset[order(df_subset$padj, decreasing = F), ]
  
  pos_enrich <- head(df_subset_sorted, min(3, nrow(df_subset_sorted)))
  pos_enrich$group <- "Positive Enrichment"
  neg_enrich <- tail(df_subset_sorted, min(3, nrow(df_subset_sorted)))
  neg_enrich$group <- "Negative Enrichment"
  
  plot_df <- rbind(pos_enrich, neg_enrich)
  plot_df <- plot_df[!duplicated(plot_df$term), ]
  
  plot_df$group <- factor(plot_df$group, levels = c("Positive Enrichment", "Negative Enrichment"))
  
  ggplot(plot_df, aes(x = fold.enrichment, y = term, fill = -log10(padj))) +
    geom_col() +
    facet_grid(rows = vars(group), scales = "free_y") + 
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 8, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(plot_df$fold.enrichment) * 1.05)) +
    scale_fill_gradient(name = "-log10(padj)", low = "blue", high = "red") +
    labs(title = paste0("Pathway Enrichment in Module ", i), x = "Enrichment", y = "")
}