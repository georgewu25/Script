# add project rpath to include access to the proj-specific packages:
proj_rpath <- "/projectnb/tcwlab/software/R/library/4.2.1"
# BiocManager::install("WGCNA", lib=proj_rpath)
.libPaths(c(proj_rpath, .libPaths()))

library(RColorBrewer)
library(edgeR)
library(WGCNA) #https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/#cranInstall
library(flashClust)
library(DESeq2)
library(parallel)
ncores <- as.numeric(Sys.getenv("NSLOTS"))
print(paste0("number of cores: ", ncores))
cl <- makeCluster(ncores) # 4)
library(doParallel)
library(org.Hs.eg.db)  #https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
library(AnnotationDbi)  #https://www.rdocumentation.org/packages/AnnotationDbi/versions/1.34.4
#https://github.com/Bioconductor/AnnotationDbi
library(ggplot2) 
library(topGO) # not instrall in R/4.2.1;
library(ReportingTools) 
library(variancePartition)
library(tidyverse)
library(magrittr)
library(stringr)

# the following shall be done with separate R code
# and call once outside the pipeline
#
# prepare annotation file, this shall be done only once:
#Annotate
anno_dir <- "/projectnb/tcwlab/Genome/hg38/"
# # better to store it in anno_dir, but due to permission, put it in meta_dir:
# t2g_file <- paste0(meta_dir, "gencode.v26.annotation.t2g.rds")
t2g_file <- paste0(anno_dir, "gencode.v26.annotation.t2g.rds")

anno_file <- paste0(anno_dir, "gencode.v26.annotation.gtf")

t2g <- read_tsv(anno_file,
                skip = 5,
                col_names = FALSE)

# get all gene info: gene id and gene name:
t2g %<>%
  rename(feature = X3, attribute = X9) %>%
  filter(feature == "gene") %>%
  mutate(gene_id = str_match(attribute, "gene_id \"(\\S+)\";")[, 2],
         gene_name = str_match(attribute, "gene_name \"(\\S+)\";")[, 2],) %>%
  select(gene_id, gene_name)

# save to rds:

saveRDS(t2g, t2g_file)

