deg_df <- data.table(
  Comparison = character(length(deseq2_name_list)),
  DEG_count = integer(length(deseq2_name_list)),
  Up = integer(length(deseq2_name_list)),
  Down = integer(length(deseq2_name_list))
)

for (i in 1:length(deseq2_name_list)) {
  
  dir <- paste0(deseq2_root, deseq2_conf_list[i])
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
  
  #p_value threshold
  df_subset <- df[!is.na(df$padj) & df$padj <= 0.1, ]
  
  deg_df[i, Comparison := deseq2_name_list[i]]
  deg_df[i, DEG_count := nrow(df_subset)]
  deg_df[i, Up := sum(df_subset$log2FoldChange > 0)]
  deg_df[i, Down := sum(df_subset$log2FoldChange < 0)]
}


#Takes the index of the deseq2_name_list
plot_deg <- function(x) {
  deg_df <- data.table(
    Comparison = character(length(x)),
    DEG_count = integer(length(x)),
    Up = integer(length(x)),
    Down = integer(length(x))
  )
  
  for (i in 1:length(x)) {
    
    index <- x[i]
    comparison <- deseq2_conf_list[index]
    name <- deseq2_name_list[index]
    
    dir <- paste0(deseq2_root, comparison)
    subdir <- list.dirs(dir)[2]
    df <- read.csv(paste0(subdir, "/", "Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
    
    #p_value threshold
    df_subset <- df[!is.na(df$padj) & df$padj <= 0.1, ]
    
    deg_df[i, Comparison := name]
    deg_df[i, DEG_count := nrow(df_subset)]
    deg_df[i, Up := sum(df_subset$log2FoldChange > 0)]
    deg_df[i, Down := sum(df_subset$log2FoldChange < 0)]
  }
  
  #Melt data
  plot_df <- tidyr::pivot_longer(deg_df, cols = c(Up, Down), names_to = "DEG", values_to = "value")
  plot_df$DEG <- factor(plot_df$DEG, levels = c("Up", "Down"))
  
  return(plot_df)
}

plot_df <- plot_deg(c(1,2,3,4))

plot_df$Genotype <- ifelse(grepl("APOE 33", plot_df$Comparison), "APOE 33", "APOE 44")
plot_df$Genotype <- factor(plot_df$Genotype, levels = c("APOE 33", "APOE 44"))

plot_df$Time <- sapply(plot_df$Comparison, function(x){
  paste0(unlist(strsplit(x, "_"))[2], " vs 0h")
})
plot_df$Time <- factor(plot_df$Time, levels = c("8h vs 0h", "24h vs 0h"))

ggplot(plot_df, aes(x = Genotype, y = value, fill = DEG)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_grid(~ Time) +
  scale_fill_manual(values = c("red", "blue"), name = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.6, size = 16),
        axis.text.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.2),
        panel.margin = unit(0, "lines"),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = "lightblue", color = "black"),
        plot.margin = margin(1,10,1,1)) +
  scale_x_discrete(expand = c(0.6, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2000 * 1.05)) +
  labs(title = "", x = "", y = "Number of DEGs (FDR<0.1)")

plot_df <- plot_deg(c(5,6,7))

plot_df$Time <- ifelse(grepl("Ctrl", plot_df$Comparison), "0 hr", ifelse(grepl("8h", plot_df$Comparison), "8 hr", "24 hr"))
plot_df$Time <- factor(plot_df$Time, levels = c("0 hr", "8 hr", "24 hr"))

ggplot(plot_df, aes(x = Time, y = value, fill = DEG)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("red", "blue"), name = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.6, size = 16),
        axis.text.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.2),
        panel.margin = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", color = "black")) +
  scale_x_discrete(expand = c(0.5, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 250)) +
  labs(title = "APOE 44 vs APOE 33", x = "", y = "Number of DEGs (FDR<0.1)")

plot_df <- plot_deg(c(8,9))

plot_df$Genotype <- ifelse(grepl("APOE 33", plot_df$Comparison), "APOE 33", "APOE 44")
plot_df$Genotype <- factor(plot_df$Genotype, levels = c("APOE 33", "APOE 44"))

ggplot(plot_df, aes(x = Genotype, y = value, fill = DEG)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("red", "blue"), name = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.6, size = 16),
        axis.text.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.2),
        panel.margin = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", color = "black")) +
  scale_x_discrete(expand = c(1, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 12)) +
  labs(title = "D24h vs D0h", x = "", y = "Number of DEGs (FDR<0.1)")

plot_df <- plot_deg(c(7,10))

plot_df$Time <- c(rep("D0h", 2), rep("D24h", 2))

ggplot(plot_df, aes(x = Time, y = value, fill = DEG)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("red", "blue"), name = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.6, size = 16),
        axis.text.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.2),
        panel.margin = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", color = "black")) +
  scale_x_discrete(expand = c(0.5, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 450)) +
  labs(title = "APOE 44 vs APOE 33", x = "", y = "Number of DEGs (FDR<0.1)")