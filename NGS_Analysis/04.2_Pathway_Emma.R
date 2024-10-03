source("/scripts/Emmaplot.R")

run_emma <- function(x){
  comparison <- deseq2_conf_list[x]
  dir <- "/output/GSEA/"
  subdir <- paste0(dir, comparison, "/")
  gsea_res <- paste0(subdir, comparison, "_c2.cp.v2023.1.Hs.symbols.gmt_result.rds")
  df <- readRDS(gsea_res)
  
  outdir <- '/output/Pathway/'
  
  #Top 40 pathways in either directions
  activated <- head(df[order(df$NES, decreasing = TRUE), ], 20)
  suppressed <- head(df[order(df$NES, decreasing = FALSE), ], 20)
  
  top_pathways <- data.frame()
  top_pathways <- rbind(activated, suppressed)
  
  emma <- emmaplot(top_pathways, label.size = 3)
  
  emma_dir <- paste0(outdir, comparison, "/", "c2.cp.v2023.1.Hs.symbols.gmt/")
  if(!file.exists(emma_dir)) {
    dir.create(emma_dir, mode="0755", recursive=TRUE)
  }
  
  #Save as pdf
  ggsave(paste0(emma_dir, comparison, "_c2_emma_plot.pdf"), height = 10, width = 20)
}

for (i in 1:length(deseq2_conf_list)) {
  run_emma(i)
}