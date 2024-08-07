proj_dir <- "/projectnb/tcwlab/yshen16/Project/APOE22/"
fc_dir <- paste0(proj_dir,  "FeatureCount/")
fc_file <- paste0(fc_dir, "featureCounts.txt")
countdata <- read.table(fc_file, header=TRUE)
countdata <- countdata[-c(2:6)]
colnames(countdata) <- gsub("\\.[sA]ligned$", "", colnames(countdata), perl=T)
head(countdata)
tmp_fc_countdata <- paste0(fc_dir, "temp.txt") 
write.table(countdata, file = tmp_fc_countdata, row.name = F, col.name = T, quote = F, sep = "\t")
