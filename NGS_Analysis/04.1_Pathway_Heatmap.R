
top_pathway_df_all <- data.frame()

for (i in 1:length(deseq2_conf_list)) {
  comparison <- deseq2_conf_list[i]
  for (j in 1:nrow(pathway_df)) {
    pathway_res <- paste0(gsea_root, comparison, "/", comparison, "_", pathway_df$gmt[j], "_result.rds")
    df <- readRDS(pathway_res)
    df_sorted <- df[order(df$NES), ]
    top_neg_pathway <- head(df_sorted, 5)
    top_pos_pathway <- tail(df_sorted, 5)
    
    top_pathways <- rbind(top_neg_pathway, top_pos_pathway)
    top_pathways$Pathway <- pathway_df$name[j]
    top_pathways$Comparison <- deseq2_conf_list[i]
    
    #Get a single data frame that contains all top pathways from all comparisons
    top_pathway_df_all <- rbind(top_pathway_df_all, top_pathways)
    
    #Pathway visualization
    outdir <- '/output/Pathway/'
    
    plot_dir <- paste0(outdir, comparison, "/", pathway_df$gmt[j], "/")
    if(!file.exists(plot_dir)) {
      dir.create(plot_dir, mode="0755", recursive=TRUE)
    }
    
    #Bar plot
    bar <- ggplot(top_pathways, aes(x = NES, y = pathway, fill = -log10(padj))) +
      geom_col() +
      scale_fill_gradientn(colors = c("yellow", "red")) +
      theme_light() +
      theme(
        panel.margin = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.text.y = element_text(size = 14)
      ) +
      labs(title = deseq2_name_list[i], x = "NES", y = "Pathway")
    
    #Save as pdf
    pdf(paste0(plot_dir, deseq2_name_list[i], "_", pathway_df$gmt[j], "_bar_plot.pdf"), height = 10, width = 15)
    print(bar)
    dev.off()
  }
}

#Takes index of deseq_list
pathway_bar_plot <- function(x) {
  df <- readRDS(paste0("/output/GSEA/", deseq2_conf_list[x], "/", deseq2_conf_list[x], "_c2.cp.v2023.1.Hs.symbols.gmt_result.rds"))
  
  df_sorted <- df[order(df$NES), ]
  top_neg_pathway <- head(df_sorted, 5)
  top_pos_pathway <- tail(df_sorted, 5)
  
  top_pathways <- rbind(top_neg_pathway, top_pos_pathway)
  
  ggplot(top_pathways, aes(x = NES, y = pathway, fill = -log10(padj))) +
    geom_col() +
    scale_fill_gradientn(colors = c("yellow", "red")) +
    theme_light() +
    theme(
      panel.margin = unit(.05, "lines"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.y = element_text(size = 8)
    ) +
    labs(title = deseq2_name_list[x], x = "NES", y = "")
}


ggplot(top_pathways, aes(x = NES, y = pathway, fill = -log10(padj))) +
  geom_col() +
  scale_fill_gradientn(colors = c("yellow", "red")) +
  theme_light() +
  theme(
    panel.margin = unit(.05, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.y = element_text(size = 10)
  ) +
  labs(title = "APOE 44_8h_vs_Ctrl", x = "NES", y = "")

#Table for Lu
for (i in c(5,6,7,10)) {
  
  top_pathways <- data.frame()
  
  comparison <- deseq2_conf_list[i]
  name <- deseq2_name_list[i]
  
  dir <- "/output/GSEA/"
  
  df <- readRDS(paste0(dir, comparison, "/", comparison, "_c2.cp.v2023.1.Hs.symbols.gmt_result.rds"))
  
  df_sorted <- df[order(df$NES), -8]
  top_neg_pathway <- head(df_sorted, 5)
  top_pos_pathway <- tail(df_sorted, 5)
  
  top_pathways <- rbind(top_neg_pathway, top_pos_pathway)
  
  write.csv(top_pathways, paste0("/output/manuscript_data/table/", name, "_Top_Pathways.csv"),)
}


#Takes index of deseq2_conf_list
plot_overlap <- function(x,y,filename) {
  df <- data.frame()
  
  #Read GSEA results
  df1 <- readRDS(paste0(gsea_root, "/", deseq2_conf_list[x], "/", deseq2_conf_list[x], "_c2.cp.v2023.1.Hs.symbols.gmt_result.rds"))
  df1$condition <- deseq2_name_list[x]
  
  df2 <- readRDS(paste0(gsea_root, "/", deseq2_conf_list[y], "/", deseq2_conf_list[y], "_c2.cp.v2023.1.Hs.symbols.gmt_result.rds"))
  df2$condition <- deseq2_name_list[x]
  
  merged_df <- merge(df1, df2, by = "pathway", all = T)
  
  #Keep only pathways that are significant FDR < 0.1
  merged_df_filtered <- merged_df %>%
    dplyr::filter((padj.x < 0.1 & (padj.y < 0.1 | is.na(padj.y))) | (padj.y < 0.1 & (padj.x < 0.1 | is.na(padj.x))))
  
  #Change NA to 0 for NES
  merged_df_filtered$NES.x[is.na(merged_df_filtered$NES.x)] <- 0
  merged_df_filtered$NES.y[is.na(merged_df_filtered$NES.y)] <- 0
  
  round_axis <- function(x) {
    lower_half <- floor(x * 2) / 2
    upper_half <- ceiling(x * 2) / 2
    
    # Compare the absolute values and return the one with the larger absolute value
    if (abs(lower_half) > abs(upper_half)) {
      return(lower_half)
    } else {
      return(upper_half)
    }
  }
  
  axis_min <- round_axis(min(min(merged_df_filtered$NES.x), min(merged_df_filtered$NES.y)))
  axis_max <- round_axis(max(max(merged_df_filtered$NES.x), max(merged_df_filtered$NES.y)))
  
  ggplot(merged_df_filtered, aes(x = NES.x, y = NES.y)) +
    geom_point() +
    geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 0.8) +
    geom_hline(yintercept = 0, color = "black", linetype = "solid", linewidth = 0.8) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
    ) +
    scale_x_continuous(limits = c(axis_min, axis_max), breaks = seq(axis_min, axis_max, by = 0.5)) +
    scale_y_continuous(limits = c(axis_min, axis_max), breaks = seq(axis_min, axis_max, by = 0.5)) +
    labs(x = "APOE 33 NES", y = "APOE 44 NES", title = filename)
  
  #Print the pathways in the four quadrants for manual annotation
  print(merged_df_filtered %>%
          dplyr::filter(NES.x < 0 & NES.y > 0) %>%
          select(pathway, NES.x, NES.y))
  
  print(merged_df_filtered %>%
          dplyr::filter(NES.x < 0 & NES.y < 0) %>%
          select(pathway, NES.x, NES.y))
  
  print(merged_df_filtered %>%
          dplyr::filter(NES.x > 0 & NES.y < 0) %>%
          select(pathway, NES.x, NES.y))
  
  print(merged_df_filtered %>%
          dplyr::filter(NES.x > 0 & NES.y > 0) %>%
          select(pathway, NES.x, NES.y))
}