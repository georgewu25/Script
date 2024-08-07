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


run_deseq2_analysis <- function () {
  conf <- read.csv(deseq2_conf, header=T)
  for (i in 1:nrow(conf)) {
    refGroup <- str_trim(conf$Reference[i])
    targetGroup <- str_trim(conf$Target[i])
    comparison <- str_trim(conf$Comparison[i])
    namespecifier <- str_trim(conf$NameSpecifier[i])
    cmp_cols <- unlist(strsplit(comparison, '+', fixed=T)) 
    refSubGroup <- unlist(strsplit(refGroup, '+', fixed=T))
    targetSubGroup <- unlist(strsplit(targetGroup, '+', fixed=T))
    ref_str <- paste(refSubGroup, collapse ="_")
    target_str <- paste(targetSubGroup, collapse ="_")  
    if(is.null(namespecifier)) {
      namespecifier <- paste0(target_str, "vs", ref_str)
    }
# use 'namespecifier in deseq conf to create the subdir:
    deseq2_subdir<-paste0(deseq2_dir, "/", namespecifier, "/")
    if(dir.exists(deseq2_subdir)) {
      unlink(deseq2_subdir, recursive = TRUE)
    }
    dir.create(deseq2_subdir, mode="0755", recursive=TRUE)
    
    
    # #drop samples 14 and 15
    # cova <-cova[ which (cova$s1to13=="yes"),]
    # exprs <-exprs[,colnames(exprs)%in%as.character(cova$id)]
    
    #select the data
    #select the ref data
    ref_covaA <- ref_cova %>% 
      filter(SampleGroup %in% refSubGroup) %>%
      select(all_of(select_cols))

    ref_exprA<- ref_exprs[,colnames(ref_exprs)%in%as.character(ref_covaA$sample_id)]

    #select the target data
    target_covaA <- target_cova %>% 
      filter(SampleGroup %in% targetSubGroup)%>%
      select(all_of(select_cols))

    target_exprA<- target_exprs[,colnames(target_exprs)%in%as.character(target_covaA$sample_id)]
    #  head(target_exprA)

    # join two data sets:
    exprA <- cbind(target_exprA, ref_exprA)
    covaA <- rbind(target_covaA, ref_covaA)
    
    ###Remove the low expression genes
    ###Remove the low expression genes
    isexpr <- rowSums(cpm(exprA)>1) >= 0.1 * ncol(exprA)
    nolowexprA <- exprA[isexpr,]

    # preprocess variates, convert to factor, scale and center, etc,  first: 
    for (j in 1:length(cmp_cols))  {
      switch(cmp_cols[j],
             SampleGroup={
               covaA$SampleGroup <- factor(
                 ifelse(covaA$SampleGroup %in% refSubGroup[1], ref_str, target_str), 
                 levels=c(ref_str, target_str))
             },
             SEX={
               covaA$SEX <- as.factor(covaA$SEX)
             }, 
             GRSnoAPOE={
               covaA$GRSnoAPOE <- as.factor(covaA$GRSnoAPOE)
             },
             RIN={
               covaA$RIN <- as.factor(covaA$RIN)
             },
             batch={
               covaA$batch <- as.factor(covaA$batch)
             },
             disease={
               covaA$disease <- as.factor(covaA$disease)
             },
             A_P_O_E={
               covaA$A_P_O_E <- as.factor(covaA$A_P_O_E, 
                                          levels=c(unique(ref_covaA$A_P_O_E), unique(target_covaA$A_P_O_E))
               )
             },
             # CellType={
             #   covaA$CellType <- as.factor(covaA$CellType)
             # },
             # individual={
             #   covaA$individual <- as.factor(covaA$individual)
             # },
             # proportion={
             #   covaA$proportion <- as.factor(covaA$proportion)
             # },
             {
               print(paste0("Warning: Unknown variate, ", cmp_cols[j]))
             }
      )
    }
    
    print (paste("i=", i, "comparison=", comparison, "nameSpecifier=", namespecifier))

    #  comparison <- "SEX+RIN+disease+A_P_O_E"
    designForm <- formula(paste0("~", comparison))
    # # check variate 
    # varPartResidA <- fitExtractVarPartModel(nolowexprA, designForm, covaA)
    # # here we need to check if the variable returns big enough differences, if so, remove that variable and report it.
    # varPartPDF <- paste0(deseq2_subdir, 'VariancePartion_Covariates.pdf')
    # #    varPartPDF <- "/projectnb/tcwlab/LabMember/yshen16/Project/APOE22/00_fastq_100/DESeq2/4V_Tvs44_M/VariancePartion_Covariates.pdf"
    # pdf(file=varPartPDF)
    # main_label <- paste0("Variance Partitioning on raw counts, Design=", comparison)
    # plotVarPart(varPartResidA ,main=main_label)
    # # main_label <- "Variance Partitioning on raw counts"    
    # # plotVarPart(varPartResidA)
    # dev.off()
    # 
    tryCatch(
      {
#      stop("Calling DESeqDataSetFromMatrix"), 
      xddsMatnolowA <- DESeqDataSetFromMatrix(nolowexprA, covaA, design = designForm)  
      xddsMatA <- DESeqDataSetFromMatrix(exprA, covaA, design=designForm)
      
      #  xddsMatnolowA <- DESeqDataSetFromMatrix(nolowexprA, covaA, design = ~ SampleGroup)  
      
      
      # ###define factor levels: levels = c(reference, target)#####
      # xddsMatnolowA$A_P_O_E <- factor(xddsMatnolowA$A_P_O_E, levels = c("3_3","4_4"))
      
      
      ############DeSeq2 differential expression #############
      
      yddsnolowA<-DESeq(xddsMatnolowA, parallel=TRUE)
      resnolowA <- results(yddsnolowA)
      resnolowAordered <- resnolowA[order(resnolowA$padj),]
      #write.csv(as.data.frame(resnolowAordered), file="Astrocyte_APOE_Results.csv")
      
      # annotate DESeq2 result according to gene id and gene name:
      annotation <- t2g[match(rownames(resnolowAordered), t2g$gene_id),]
      
      # bind reslult and save a copy: 
      resnolowAorderedAnnotated <- cbind(resnolowAordered,annotation)
      nolowAorderedAnno_file <- paste0(deseq2_subdir, "Results_GTFAnnotated.csv")
      write.csv(as.data.frame(resnolowAorderedAnnotated), file=nolowAorderedAnno_file)
      
      # order result by gene name and FC Log2foldChange, save the result
      nodup <- resnolowAorderedAnnotated
      nodup$absvalFC <- abs(nodup$log2FoldChange)
      nodup <- nodup[order(nodup$gene_name,-nodup$absvalFC),]
      nodup <- nodup[!duplicated(nodup$gene_name),]
      nodup_file <- paste0(deseq2_subdir, "Results_GTFAnnotated_NoGeneIDDuplicates.csv")
      write.csv(as.data.frame(nodup), file=nodup_file) #final result to deliver
      
      ########### Normalized Counts from DESeq2 and CPM ##############
      xddsMatA_estimateSizeFactors <- estimateSizeFactors(xddsMatA)
      xddsMatA_deseq2_normalized_counts <- counts(xddsMatA_estimateSizeFactors, normalized=TRUE)
      normalized_count_file <- paste0(deseq2_subdir, "deseq2_normalized_counts.csv")
      write.csv(as.data.frame(xddsMatA_deseq2_normalized_counts), file=normalized_count_file) # final result 2 to deliver
      
      
      exprA_count_per_million <- cpm(exprA)
      exprA_cpm_file <- paste0(deseq2_subdir, "count_per_million.csv")
      write.csv(as.data.frame(exprA_count_per_million), file=exprA_cpm_file)
      
      ########### MA (Mean Average) plot #############
      folderChg_pdf <- paste0(deseq2_subdir, "FoldChange.pdf")
      pdf(folderChg_pdf) # final result 4
      label_str <- paste(targetGroup, " vs ", refGroup, ", design=~", comparison)
      plotMA(resnolowA,ylim=c(-5,5),main=label_str)
      dev.off()
      },
      error = function(e) 
      {
        print(e$message)
      }
    )
    
  } # end each deseq conf
  
}  # end of run_deseq2_analysis

# START of Main Program entry:


proj_dir <- "{PROJECT_DIR}"
if(! str_sub(proj_dir, -1, -1) == "/") {
  proj_dir <- paste0(proj_dir, "/")
}
#DESEQ2_CONF_NUM:3
#DESEQ2_CONF_LIST:tcw1,tcw2,tcw1_tcw2
deseq2_conf_num <- {DESEQ2_CONF_NUM}
deseq2_conf_liststr <- "{DESEQ2_CONF_LIST}"
deseq2_conf_list <- unlist(str_split(gsub(" ", "", deseq2_conf_liststr), ","))
conf_root <- "{CONFIG_ROOT}" #"/projectnb/tcwlab/LabMember/yshen16/Project/config/APOE_Mye_Micro_trim_100"
if(! str_sub(conf_root, -1, -1) == "/") {
  conf_root <- paste0(conf_root, "/")
}
deseq2_conf_dir <- "{DESEQ2_CONF_DIR}" #"deseq2_conf"
if(! str_sub(deseq2_conf_dir, -1, -1) == "/") {
  deseq2_conf_dir <- paste0(deseq2_conf_dir, "/")
}

# now let's read in the config files and redo the analysis: 
# deconf_list <- list.files(paste0(conf_root, "/", deseq2_conf_dir))
# here let's create a new metadata to standardize the config info: 
meta_root <- paste0(proj_dir,  "metadata", "/") # make sure all the directory end with "/"
deseq2_root <- paste0(proj_dir, "DESEQ2", "/") 

# load gene anno list: 
anno_dir <- "/projectnb/tcwlab/RawData/Genome2.7.9a/hg38/"
t2g_file <- paste0(anno_dir, "gencode.v26.annotation.gtf")
t2g <- read.table(t2g_file, sep = "\t", header = TRUE)

# loop through the config files:
for (c in 1:deseq2_conf_num) {
  cfg_name <- deseq2_conf_list[c]
  meta_dir <- paste0(meta_root, cfg_name, "/");
  # remove existing metadata/cfg_name folder
  if(dir.exists(meta_dir)) {
   unlink(meta_dir, recursive = TRUE, force=1)
  }
  dir.create(meta_dir, mode="0755", recursive=TRUE)
  orig_cfg_dir <- paste0(conf_root, deseq2_conf_dir, cfg_name, "/")
  cfg_file <- paste0(orig_cfg_dir, "/", cfg_name, ".cfg")
  cfg_data <- readLines(cfg_file)
  cfg_data <- str_replace_all(cfg_data, " ", "")
  self_ref <- as.numeric(str_split(cfg_data[grep("SELF_REF:", cfg_data)], ":")[[1]][2])
  if(self_ref == 1) {
    # META:tcw1_metadata_micro_APOE22_w_group.csv
    # DESEQ2_CONF:tcw1_deseq2_config.csv
    # META_COL:sample_id,RIN,SEX,CellType,A_P_O_E,individual,GRSnoAPOE,proportion,batch,select,SampleGroup
    meta_orig <- str_split(cfg_data[grep("META:", cfg_data)], ":")[[1]][2]
    ref_meta_file <- paste0(meta_dir, "ref_metadata.csv")
    file.symlink(paste0(orig_cfg_dir, meta_orig), ref_meta_file)
    target_meta_file <- paste0(meta_dir, "target_metadata.csv")
    file.symlink(paste0(orig_cfg_dir, meta_orig), target_meta_file)
    ref_fc_file <- paste0(meta_dir, "/ref_fc.txt")
    file.symlink(paste0(proj_dir, "FeatureCount/featureCounts_clean.txt"), ref_fc_file)
    target_fc_file <- paste0(meta_dir, '/target_fc.txt')
    file.symlink(paste0(proj_dir, "FeatureCount/featureCounts_clean.txt"), target_fc_file)
  } else { #self_ref==0
    # TARGET_META:metadata_micro_APOE22_w_group.csv
    # REF_FC:/projectnb/tcwlab/Project/PopulationAPOE/JuliaNote/Alignment_final/featureCounts/featureCounts_clean.txt
    # REF_META:RNAseqPhen_sep6_v2_APOE33.csv
    # DESEQ2_CONF:deseq2_config.csv
    ref_meta_orig <- str_split(cfg_data[grep("REF_META:", cfg_data)], ":")[[1]][2]
    ref_meta_file <- paste0(meta_dir, "ref_metadata.csv")
    file.symlink(paste0(orig_cfg_dir, ref_meta_orig), ref_meta_file)
    target_meta_orig <- str_split(cfg_data[grep("TARGET_META:", cfg_data)], ":")[[1]][2]
    target_meta_file <- paste0(meta_dir, "target_metadata.csv")
    file.symlink(paste0(orig_cfg_dir, target_meta_orig), target_meta_file)
    ref_fc_orig <- str_split(cfg_data[grep("REF_FC:", cfg_data)], ":")[[1]][2]
    ref_fc_file <- paste0(meta_dir, "ref_fc.txt")
    file.symlink(paste0(orig_cfg_dir, ref_fc_orig), ref_fc_file)
    target_fc_file <- paste0(meta_dir, 'target_fc.txt')
    file.symlink(paste0(proj_dir, "FeatureCount/featureCounts_clean.txt"), target_fc_file)
  }
  deseq2_dir <- paste0(deseq2_root, cfg_name, "/")  
  if(!dir.exists(deseq2_dir)) {
    dir.create(deseq2_dir, mode="0755", recursive = TRUE)
  }
  deseq2_orig <- str_split(cfg_data[grep("DESEQ2_CONF:", cfg_data)], ":")[[1]][2]
  deseq2_conf <- paste0(meta_dir, "deseq2_config.csv")
  file.symlink(paste0(orig_cfg_dir, deseq2_orig), deseq2_conf)
  meta_col<- str_split(cfg_data[grep("META_COL:", cfg_data)], ":")[[1]][2]
  select_cols <- unlist(str_split(meta_col,","))
  
  ######read data in matrix format####
  # load feature counts: 
  ref_exprs <- as.matrix(read.table(ref_fc_file, header=TRUE, sep="\t", row.names=1,as.is=TRUE))
  head (ref_exprs)
  target_exprs <- as.matrix(read.table(target_fc_file, header=TRUE, sep="\t", row.names=1,as.is=TRUE))
  head (target_exprs)
  
  # load the entire meta data 
  ref_cova <- read.csv(ref_meta_file, header=TRUE, stringsAsFactors=FALSE, colClasses="character")
  # %>%
  #             filter(SampleGroup=="33_M") %>%
  #             select(id,RIN,SEX,A_P_O_E,batch,disease, SampleGroup)
  # head (ref_cova)
  
  
  target_cova <- read.csv(target_meta_file, header=TRUE, as.is=TRUE, stringsAsFactors=FALSE, colClasses="character")
  # head(target_cova)
  run_deseq2_analysis()
}





