library(apeglm)
library(ComplexHeatmap)
library(cowplot)
library(data.table)
library(dittoSeq)
library(dplyr)
library(flashClust)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(patchwork)
library(pheatmap)
library(plotly)
library(RColorBrewer)
library(Seurat) 
library(scales)
library(scCustomize)
library(stringr)
library(tidyverse)

#LOAD LABELED DATASET SEURAT OBJECT
load(file = paste(dir_data, "Prelim_Analysis/Cell_Type_Annotation/FCT", sep = ""))

#Annotation with Gene Symbol
t2g <- read_tsv(paste(dir_data, "SingleR_Reference_Data/gencode.v26.annotation_1.gtf", sep = ""), skip = 5, col_names = FALSE)

t2g %<>%
  dplyr::rename(feature = X3, attribute = X9) %>%
  filter(feature == "gene") %>%
  mutate(gene_id = str_match(attribute, "gene_id \"(\\S+)\";")[, 2],
         gene_name = str_match(attribute, "gene_name \"(\\S+)\";")[, 2],) %>%
  select(gene_id, gene_name)

t2g %>% filter(is.na(gene_id)) %>% nrow

#Get metadata table
metadata <- Percent_Expressing(Seurat_SCT,
  features = feats,
  threshold = 0,
  group_by = "Final_Cell_Type",
  split_by = NULL,
  entire_object = FALSE,
  slot = "data",
  assay = "SCT")


#Get cell types
harmony_cell_type_list <- as.vector(unique(Seurat$harmony_FCT))

#Exclude cell types that do not have enough cells in each group for comparison
harmony_cell_type_list <- harmony_cell_type_list[harmony_cell_type_list %in% c("NPC 3", "VLMC", "OPC + Oligodendrocyte", "Astrocyte 2") == FALSE]

#Run DEA (MAST & Wilcoxon)
for (set in sets) { #set = All, individual_1, or individual_2
  
  #create data table for DEA method comparison
  DT = data.table()
  wilcox_stat_meth_DT = data.table()
  seurat_mast_meth_DT = data.table()
  
  for (x in harmony_cell_type_list){
    
    Idents(set) <- set$harmony_FCT
    sub <- subset(set, idents = x)
    
    Idents(set) <- set$Individual
    
    Idents(sub) <- sub$APOE_Genotype 
    
    #Seurat Wilcoxon via FindMarkers function (only.pos = FALSE) to find differentially expressed genes (both up & down).
    Wilcoxon <- FindMarkers(sub, ident.1 = "APOE 44", ident.2 = "APOE 33", logfc.threshold = 0, min.pct = 0.1, assay="SCT", test.use = "wilcox", only.pos = FALSE)
   
      genes_to_test <- rownames(Wilcoxon)
      
      Wilcoxon_stat <- wilcoxauc(sub, 'APOE_Genotype', seurat_assay = "SCT",  groups_use = c('APOE 44', 'APOE 33')) 
      Wilcoxon_stat <- as.data.table(Wilcoxon_stat)
      
      #get correct statistic & auc
      for(i in 1:nrow(Wilcoxon_stat)){ #each row
        
        gene <- Wilcoxon_stat$feature[i] 
        
        stat_44 <- Wilcoxon_stat[Wilcoxon_stat$feature == gene & Wilcoxon_stat$group == 'APOE 44']$statistic 
        stat_33 <- Wilcoxon_stat[Wilcoxon_stat$feature == gene & Wilcoxon_stat$group == 'APOE 33']$statistic 
        
        auc_44 <- Wilcoxon_stat[Wilcoxon_stat$feature == gene & Wilcoxon_stat$group == 'APOE 44']$auc 
        auc_33 <- Wilcoxon_stat[Wilcoxon_stat$feature == gene & Wilcoxon_stat$group == 'APOE 33']$auc
        
        statistic <- min(stat_33, stat_44) #minimum of u-stat_1 & u-stat_2
        auc <- min(auc_33, auc_44)
        
        Wilcoxon_stat[Wilcoxon_stat$feature == gene & Wilcoxon_stat$group == 'APOE 44']$statistic <- statistic
        Wilcoxon_stat[Wilcoxon_stat$feature == gene & Wilcoxon_stat$group == 'APOE 44']$auc <- auc
        
      }
      
      Wilcoxon_stat <- Wilcoxon_stat[Wilcoxon_stat$feature %in% genes_to_test]
      Wilcoxon_stat <- Wilcoxon_stat[Wilcoxon_stat$group %in% 'APOE 44']
      
      Wilcoxon_stat <- Wilcoxon_stat %>% 
        dplyr::rename('gene' = 'feature') %>%  #left is what is new #rename(new_column_name = old_column_name)
        dplyr::rename('p_val' = 'pval') %>% 
        dplyr::rename('BH_p_val_adj' = 'padj') #benjamini hochberg
      Wilcoxon_stat <- Wilcoxon_stat[, ('group'):=NULL] 
      
      #BC p value adjustment
      Wilcoxon_stat$BC_p_val_adj <- p.adjust(Wilcoxon_stat$p_val, method = "bonferroni") 
      
      #get avg log 2 fold change measurement from seurat wilcoxon --> logFC to avg_log2FC
      Wilcoxon_stat <- Wilcoxon_stat %>% slice_max(n = 1000000000, order_by = p_val)
      Wilcoxon <- Wilcoxon %>% slice_max(n = 1000000000, order_by = p_val)
      
      #Wilcoxon <- subset(Wilcoxon, rownames(Wilcoxon) %in% Wilcoxon_stat$gene)
      Wilcoxon_stat$avg_log2FC <- Wilcoxon$avg_log2FC
      Wilcoxon_stat <- as.data.frame(Wilcoxon_stat) #only df have rownames
      rownames(Wilcoxon_stat) <- Wilcoxon_stat$gene
      Wilcoxon_stat$gene <- NULL
      
      # MAST DEA (applies hurdle model & removes variation due to individual)
      #latent variable = individual
      if (identical(set, Seurat)){
        MAST <- FindMarkers(sub, ident.1 = "APOE 44", ident.2 = "APOE 33", logfc.threshold = 0, min.pct = 0.1, assay="RNA", test.use = "MAST", latent.vars = "Individual", only.pos = FALSE)
      } 
      else{
        MAST <- FindMarkers(sub, ident.1 = "APOE 44", ident.2 = "APOE 33", logfc.threshold = 0, min.pct = 0.1, assay="RNA", test.use = "MAST", only.pos = FALSE)
      }
      
      MAST <- MAST %>% 
        dplyr::rename('BC_p_val_adj' = 'p_val_adj') #bonferroni correction
        
        #BH p value adjustment
        MAST$BH_p_val_adj <- p.adjust(MAST$p_val, method = "hochberg") 
      
      print("done MAST")
      
    if (identical(set, Seurat)){
      y = "All"
    } 
    
    if (identical(set, Individual_1)){
      y = "Individual 1"
    } 
    
    if (identical(set, Individual_2)){
      y = "Individual 2"
    } 
          
    write.csv(Wilcoxon_stat %>% slice_max(n = 10000000, order_by = avg_log2FC), paste("Harmony/Specific/",y,"/",x,"/wilcox_stat_dea/",x,"_wilcoxon_stat_dea.csv", sep = ''))

    write.csv(MAST %>% slice_max(n = 10000000, order_by = avg_log2FC), paste("Harmony/Specific/",y,"/",x,"/seurat_mast_dea/",x,"seurat_mast_dea.csv", sep = ''))

    #Annotation gene name
    
    for (d in list(Wilcoxon_stat, MAST)){
      
      res_ordered <- d

      annotation <- t2g[match(rownames(res_ordered), t2g$gene_name),]
      resorderedAnnotated <- cbind(res_ordered, annotation)
      
      if (identical(d, Wilcoxon_stat)){
        w = "wilcox_stat_dea"
      }
      if (identical(d, MAST)){
        w = "seurat_mast_dea"
      }
      
      write.csv(as.data.frame(resorderedAnnotated), file= paste("Harmony/Specific/",y,"/",x,"/",w,"/harmony_",x,"_",w,"_Results_GTFAnnotated.csv", sep = ''))
      
      nodup <- resorderedAnnotated
      nodup$absvalFC <- abs(nodup$avg_log2FC)
     
      nodup <- nodup[order(nodup$gene_name,-nodup$absvalFC),]
      nodup <- nodup[!duplicated(nodup$gene_name),]
      write.csv(as.data.frame(nodup), file= paste("Harmony/Specific/",y,"/",x,"/",w,"/harmony_",x,"_",w,"_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
       
      res <- nodup #table of genes (exp > 10% cells, ranked by logFC E44 vs. E33) #this has been annotated
      
      res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
      
      #FDR of 0.05 
      res_table_pval_thres <- res_tbl %>% 
        mutate(threshold = p_val < 0.05) #significant DEGs

      res_table_pval_thres_only <- res_table_pval_thres[res_table_pval_thres$p_val < 0.05,] #only look at sig DEGs based on p-value < 0.05.
      res_table_pval_thres_only <- res_table_pval_thres_only[!is.na(res_table_pval_thres_only$p_val),] 
      res_table_pval_thres_only <- as.data.frame(res_table_pval_thres_only)
      write.csv(as.data.frame(res_table_pval_thres_only), file= paste("Harmony/Specific/",y,"/",x,"/",w,"/harmony_",x,"_",w,"_significant_0.05_degs_by_p_val_only_results.csv", sep = ''))
      
      #BH: Benjamini Hochberg
      res_table_BH_padj_thres <- res_tbl %>%
        mutate(threshold = BH_p_val_adj < 0.05) #significant DEG
     
      #Select Sig DEGS only BH p val adj 
      res_table_BH_padj_thres_only <- res_table_BH_padj_thres[res_table_BH_padj_thres$BH_p_val_adj < 0.05,] #only look at significant DEGs based on BH padj-value < 0.05
      res_table_BH_padj_thres_only <- res_table_BH_padj_thres_only[!is.na(res_table_BH_padj_thres_only$BH_p_val_adj),] 
      res_table_BH_padj_thres_only <- as.data.frame(res_table_BH_padj_thres_only)
      write.csv(as.data.frame(res_table_BH_padj_thres_only), file= paste("Harmony/Specific/",y,"/",x,"/",w,"/harmony_",x,"_",w,"_significant_0.05_degs_by_BH_padj_val_only_results.csv", sep = ''))
      
      #BC: Bonferroni Correction
      res_table_BC_padj_thres <- res_tbl %>%
        mutate(threshold = BC_p_val_adj < 0.05)
      
      #Select Sig DEGS only BC p val adj 
      res_table_BC_padj_thres_only <- res_table_BC_padj_thres[res_table_BC_padj_thres$BC_p_val_adj < 0.05,] #only look at significant DEGs based on BC padj-value < 0.05
      res_table_BC_padj_thres_only <- res_table_BC_padj_thres_only[!is.na(res_table_BC_padj_thres_only$BC_p_val_adj),] 
      res_table_BC_padj_thres_only <- as.data.frame(res_table_BC_padj_thres_only)
      write.csv(as.data.frame(res_table_BC_padj_thres_only), file= paste("Harmony/Specific/",y,"/",x,"/",w,"/harmony_",x,"_",w,"_significant_0.05_degs_by_BC_padj_val_only_results.csv", sep = ''))
   
    } 
  }
} 
