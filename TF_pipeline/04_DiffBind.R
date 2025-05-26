library(DiffBind)


extsize73_dir <- "../ATAC_seq/output/peaks/macs2_extsize73"
broadpeak_files <- list.files(extsize73_dir, pattern = "_peaks\\.broadPeak$", full.names = T)

bam_dir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/ATAC_seq/output/alignment_filtered"
bam_files <- list.files(bam_dir, pattern = "\\.blacklist-filtered.bam$", full.names = T)

sample_name <- sub("_peaks.broadPeak.*", "", basename(broadpeak_files))

mtd <- data.frame(
  SampleID = sample_name,
  Tissue = NA,
  Factor = NA,
  Condition = substr(sample_name, 5,7),
  Treatment = sub(".*-(.*)", "\\1", sample_name),
  Replicate = as.numeric(substr(sample_name, 4,4)),
  bamReads = bam_files,
  ControlID = NA,
  bamControl = NA,
  Peaks = broadpeak_files,
  PeakCaller = "macs"
)


#Create DBA Object to extract overlapping peaks after merging. 
dba <- dba(sampleSheet=mtd)

#Generate read-count based score matrix
dba <- dba.count(dba)

saveRDS(dba, "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/ATAC_seq/output/DBA/DiffBind/MACS2/dba_all.rds")
dba_obj <- readRDS("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/ATAC_seq/output/DBA/DiffBind/MACS2/dba_all.rds")
dba_obj <- dba(dba_obj)

#Extract info on reads per sample and fraction of reads that align to overlaping peaks
info <- dba.show(dba_obj)

libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, PeakReads=round(info$Reads * info$FRiP))

#Normalization
dba_obj <- dba.normalize(dba_obj)
dba_normalized <- dba.normalize(dba_obj, bRetrieve=T)

normlibs <- cbind(FullLibSize=dba_normalized$lib.sizes, NormFacs=dba_normalized$norm.factors, NormLibSize=round(dba_normalized$lib.sizes/dba_normalized$norm.factors))

run_dba_within_genotype_analysis <- function(DBA, genotype, target, control) {
  
  #Define contrast
  dba_obj_subset <- dba(dba_obj, mask=dba_obj$samples$Condition == genotype)
  contrast <- dba.contrast(dba_obj_subset, design="~Treatment", contrast=c("Treatment", target, control), 
                        reorderMeta=list(Treatment=c(control)),
                        minMembers=2, bNot=F, bComplex=F, bGetCoefficients=F)
  
  dba_final <- dba.analyze(contrast)

  dba_report <- dba.report(dba_final)
  
  return(dba_report)
}

for (i in seq_along(dba_output)) {
  
  file <- readRDS(dba_output[i])
  
  write.table(file, paste0(macs2_dba_bed_dir, "/", dba_name[i], ".bed"), row.names = F, quote = F, sep = "\t")
}

