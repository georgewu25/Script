library(ChIPseeker)
library(clusterProfiler)
#BiocManager::install("ReactomePA")
library(ReactomePA)
library(ggplot2)


library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(msigdbr)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

outdir <- "../ATAC_seq/output/DBA/DiffBind/MACS2"

files <- list.files(outdir, recursive = T, full.names = T)

file_name <- list.files(outdir, recursive = T, full.names = F)
comparison <- sub(".rds", "", file_name)
comparison <- sub("dba_", "", comparison)


run_chipseeker <- function(index) {
  
  obj <- readRDS(files[index])
  
  name <- comparison[index]
  
  if (length(obj) == 0) {
    print(paste0("Not processing ", name))
  } else{
    
    print(paste0("Processing: ", name))
    
    outdir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/ATAC_seq/output/peak_anno/macs2"
    subdir <- paste0(outdir,"/", name)
    if (!exists(subdir)) {
      dir.create(subdir, mode="0755", recursive=TRUE)
    }
    
    #3kb upstream and downstream of TSS
    peakAnno <- annotatePeak(obj, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Hs.eg.db")
    
    if (length(peakAnno) != 0) {
      
      df <- as.data.frame(peakAnno)
      write.csv(df, paste0(subdir, "/", name, "_annotated_dba_peaks.csv"))
      
      #Plot pie chart of feature distribution
      ggplot(peakAnno@annoStat, aes(x="", y=Frequency, fill=Feature)) +
        geom_bar(stat = "identity", width = 3) +
        coord_polar(theta = "y") +
        scale_fill_manual(values = feature_colors) +
        theme_void() +
        ggtitle(paste0(name, " Feature Distribution"))
      ggsave(paste0(subdir, "/", name, "_Feature_Pie_chart.pdf"), width = 10, height = 8)
      
      
      #Plot pathways
      pathway <- enrichPathway(as.data.frame(peakAnno)$geneId)
      
      if (!is.null(pathway) & nrow(as.data.frame(pathway)) > 1 ) {
        dotplot(pathway)
        ggsave(paste0(subdir, "/", name, "_pathway_plot.pdf"), width = 8, height = 8)
      } else{
        print(paste0("No pathways for ", name))
      }
      
    } #End peakanno length check
  } #End GRange length check
} #End function

feature_colors <- c(
  "Promoter (<=1kb)" = "#E74C3C",  # Red (Promoter 1)
  "Promoter (1-2kb)" = "#C0392B",  # Darker Red (Promoter 2)
  "Promoter (2-3kb)" = "#F1948A",  # Lighter Red (Promoter 3)
  
  "5' UTR" = "#5DADE2",            # Blue (5' UTR)
  "3' UTR" = "#1F77B4",            # Darker Blue (3' UTR)
  
  "1st Exon" = "#28B463",          # Green (1st Exon)
  "Other Exon" = "#52BE80",        # Lighter Green (Other Exon)
  "1st Intron" = "#58D68D",        # Green (1st Intron)
  
  "Other Intron" = "#1ABC9C",      # Dark Green (Other Intron)
  
  "Downstream (<=300)" = "#BDC3C7",# Light Gray (Downstream)
  "Distal Intergenic" = "#7F8C8D"  # Darker Gray (Distal Intergenic)
)

for (i in seq_along(files)) {
  run_chipseeker(i)
}