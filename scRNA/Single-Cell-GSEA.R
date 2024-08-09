library(apeglm)
library(ComplexHeatmap)
library(cowplot)
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


#FGSEA Prep
harmony_cell_type_list <- as.vector(unique(Seurat$harmony_FCT)) #explore APOE DEGs within the specific cell types
harmony_cell_type_list <- harmony_cell_type_list[harmony_cell_type_list %in% c("NPC 3", "VLMC", "OPC + Oligodendrocyte", "Astrocyte 2") == FALSE]

gen_harmony_cell_type_list <- as.vector(unique(Seurat$general_harmony_FCT))
gen_harmony_cell_type_list <- gen_harmony_cell_type_list[gen_harmony_cell_type_list %in% c("VLMC", "OPC + Oligodendrocyte") == FALSE]

for (t in c("Specific", "General")){
  
  if (identical(t, "Specific")){
    #cell_list <- harmony_cell_type_list
    cell_list <- c("Astrocyte 1", "Astrocyte 3", "Astrocyte 5", "Astrocyte 6")
  }
  
  if (identical(t, "General")){
    #cell_list <- gen_harmony_cell_type_list
    cell_list <- "Astrocyte"
  }
  
  for (x in cell_list){ #x = cell type
    
    #annotated no dup files, no p val adj threshold, no logfc threshold: all degs found via preso wilcox method
    
    wilcox_1 <- read_csv(paste("Harmony/",t,"/individual_1/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
    wilcox_2 <- read_csv(paste("Harmony/",t,"/individual_2/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
    
    wilcox_1 <- as.data.frame(wilcox_1)
    wilcox_2 <- as.data.frame(wilcox_2)
    
    #get only those genes in both individual 1 and individual 2
    com_degs <- c(unique(wilcox_1$gene_name), unique(wilcox_2$gene_name))
    com_degs <- as.data.table(com_degs)
    com_degs <- com_degs[duplicated(com_degs)] #gene names only present in both
    com_degs <- com_degs[!(is.na(com_degs))] #remove NA
    
    #build table
    colnames(com_degs)<- "gene_name"
    com_degs[, statistic := 0] #wilcoxon u-stat
    com_degs[, auc := 0] #wilcoxon auc --> for magnitude of fgsea stat
    com_degs[, avg_log2FC := 0] #log2fc --> for sign of fgsea stat
    com_degs[, fgsea_stat := 0] #fgsea_stat 
    com_degs[, gene_id := ""] #gene annotation 
    

    if (length(com_degs$gene_name) != 0){
      
      for (gene_test in com_degs$gene_name){ #extract stats for each deg from each ind
        
        stat_1 <- (subset(wilcox_1, gene_name == gene_test))$statistic #individual 1 u-stat
        stat_2 <- (subset(wilcox_2, gene_name == gene_test))$statistic #individual 2 u-stat
        auc_1 <- (subset(wilcox_1, gene_name == gene_test))$auc #individual 1 auc
        auc_2 <- (subset(wilcox_2, gene_name == gene_test))$auc #individual 2 auc
        logFC_1 <- (subset(wilcox_1, gene_name == gene_test))$avg_log2FC #individual 1 avg_log2FC
        logFC_2 <- (subset(wilcox_2, gene_name == gene_test))$avg_log2FC #individual 2 avg_log2FC
        gene_anno <- (subset(wilcox_1, gene_name == gene_test))$gene_id #gene annotation
        
        x_1 = -log10(auc_1) * sign(logFC_1)
        x_2 = -log10(auc_2) * sign(logFC_2)
          
        #if (sign(logFC_1) == sign(logFC_2)){ #the gene should be regulated in the same direction in both individuals
          
        com_degs[com_degs$gene_name == gene_test]$statistic <- mean(stat_1, stat_2) #mean(stat_1, stat_2) #take the avg of the stats --> for fgsea
        com_degs[com_degs$gene_name == gene_test]$fgsea_stat <- mean(x_1, x_2) #mean(stat_1, stat_2) #take the avg of the stats --> for fgsea
        com_degs[com_degs$gene_name == gene_test]$auc <- mean(auc_1, auc_2) #take the avg of the auc
        com_degs[com_degs$gene_name == gene_test]$avg_log2FC <- mean(logFC_1, logFC_2) #take the avg of the log fold change
        com_degs[com_degs$gene_name == gene_test]$gene_id <- gene_anno #gene annotation
        #}
        
        #else{ #else, remove the gene from the combined list
          #com_degs <- com_degs[gene_name != gene]
        #}
        
      } #for each deg
      
      write.csv(as.data.frame(com_degs), file= paste("Harmony/",t,"/combined/",x,"/fgsea/harmony_",x,"wilcox_combined_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
      
    } #if there are shared degs, build the table

  } #cell type
} #cell annotation type


#BH vs. BC Decision

for (t in c("Specific", "General")){
  
  if (identical(t, "Specific")){
    cell_list <- harmony_cell_type_list
  }
  
  if (identical(t, "General")){
    cell_list <- gen_harmony_cell_type_list
  }
  
  DT = data.table()
  wilcox_meth_DT = data.table()
  mast_meth_DT = data.table()
  bar_DT = data.table()
  
  for (x in cell_list){ #x = cell type
  
    for (d in list("wilcox_stat_dea", "seurat_mast_dea")){
          
          #WILCOXON
          if(identical(d, "wilcox_stat_dea")){ #wilcoxon: combine
            
            #all files have been gene annotated
            
            #Individual 1
            y = "individual_1"
      
              #BH: Benjamini Hochberg
          
              wilcox_res_table_BH_padj_thres_only_1 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BH_padj_val_only_results.csv", sep = ''))
              
              wilcox_res_table_BH_padj_thres_0.15_1 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BH_padj_val_and_log2fc_0.15_results.csv", sep = ''))
              
              wilcox_res_table_BH_padj_thres_0.25_1 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BH_padj_val_and_log2fc_0.25_results.csv", sep = ''))
              
              #BC: Bonferroni Correction
  
              wilcox_res_table_BC_padj_thres_only_1 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BC_padj_val_only_results.csv", sep = ''))
              
              wilcox_res_table_BC_padj_thres_0.15_1 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BC_padj_val_and_log2fc_0.15_results.csv", sep = ''))
              
              wilcox_res_table_BC_padj_thres_0.25_1 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BC_padj_val_and_log2fc_0.25_results.csv", sep = ''))
            
              
            #Individual 2
            y = "individual_2"
            
              #BH: Benjamini Hochberg
              wilcox_res_table_BH_padj_ <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BH_padj_val_only_results.csv", sep = ''))
             
              
              #BC: Bonferroni Correction
              wilcox_res_table_BC_padj<- read_csv(paste("Harmony/",t,"/",y,"/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_significant_0.05_degs_by_BC_padj_val_only_results.csv", sep = ''))
             
              
             
            #COMBINE WILCOXON individual 1 & 2 overlapping degs
              for (method in c("BH", "BC")){ #DEA method
                for (thresh in c("only")){ 
                    
                    wilcox_1 <- get(paste0("wilcox_res_table_",method,"_padj_thres_",thresh,"_1", sep='')) #individual 1
                    wilcox_2 <- get(paste0("wilcox_res_table_",method,"_padj_thres_",thresh,"_2", sep='')) #individual 2
          
                    wilcox_1 <- as.data.frame(wilcox_1)
                    wilcox_2 <- as.data.frame(wilcox_2)
                    
                    #get only those genes present in both lists
                    com_degs <- c(unique(wilcox_1$gene_name), unique(wilcox_2$gene_name)) #all gene names 
                    com_degs <- as.data.table(com_degs)
                    com_degs <- com_degs[duplicated(com_degs)] #gene names only present in both
                    com_degs <- com_degs[!(is.na(com_degs))] #remove NA
                    
                    #build table
                    colnames(com_degs)<- "gene_name"
                    com_degs[, statistic := 0] #wilcox u-stat
                    com_degs[, avg_log2FC := 0] #log2fc
                    com_degs[, ind_1_p_val := 0] #non-adjusted p-value
                    com_degs[, ind_2_p_val := 0] #non-adjusted p-value
                    com_degs[, ind_1_p_val_adj := 0] #adjusted p-value
                    com_degs[, ind_2_p_val_adj := 0] #adjusted p-value
                    com_degs[, gene_id := ""] #gene annotation
                    
                    
                    #fill table
                    if (length(com_degs != 0)){ #must be >0 shared degs
                      
                      for (gene_test in com_degs$gene_name){ #extract stats for each deg from each ind
 
                        stat_1 <- (subset(wilcox_1, gene_name == gene_test))$statistic #individual 1 u-stat
                        stat_2 <- (subset(wilcox_2, gene_name == gene_test))$statistic #individual 2 u-stat
                        logFC_1 <- (subset(wilcox_1, gene_name == gene_test))$avg_log2FC #individual 1 avg_log2FC
                        logFC_2 <- (subset(wilcox_2, gene_name == gene_test))$avg_log2FC #individual 2 avg_log2FC
                        gene_anno <- (subset(wilcox_1, gene_name == gene_test))$gene_id #gene annotation
                        
                          if (sign(logFC_1) == sign(logFC_2)){ #the gene should be regulated in the same direction in both individuals
                            
                            com_degs[com_degs$gene_name == gene_test]$statistic <- mean(stat_1, stat_2) #mean(stat_1, stat_2) #take the avg of the stats --> for fgsea
                            com_degs[com_degs$gene_name == gene_test]$avg_log2FC <- mean(logFC_1, logFC_2) #take the avg of the log fold change
                            com_degs[com_degs$gene_name == gene_test]$ind_1_p_val <- (subset(wilcox_1, gene_name == gene_test))$p_val #ind 1 p_val
                            com_degs[com_degs$gene_name == gene_test]$ind_2_p_val <- (subset(wilcox_2, gene_name == gene_test))$p_val #ind 2 p_val
                            
                            p_val_adj = paste(method, "_p_val_adj", sep = '')
                            
                            com_degs[com_degs$gene_name == gene_test]$ind_1_p_val_adj <- (subset(wilcox_1, gene_name == gene_test))[[paste0(method,"_p_val_adj")]] #ind 1 p_val_adj 
                            com_degs[com_degs$gene_name == gene_test]$ind_2_p_val_adj <- (subset(wilcox_2, gene_name == gene_test))[[paste0(method,"_p_val_adj")]] #ind 2 p_val_adj
                            
                            com_degs[com_degs$gene_name == gene_test]$gene_id <- gene_anno #gene annotation
                          }
                          
                          else{ #else, remove the gene from the combined list
                            com_degs <- com_degs[gene_name != gene_test]
                          }
                 
                      } #for each deg
     
                    } #if there are shared degs, build the table
                    
                    #re-assign table to variable name
                    
                    assign(paste0("wilcox_res_table_",method,"_padj_thres_",thresh,sep=''), com_degs)
                    
                } #avg_log2fc threshold
              } #for method
            } #if Wilcoxon
      
          #MAST 
          if(identical(d, "seurat_mast_dea")){
           
            y = "All"
            
            #BH: Benjamini Hochberg
            
            mast_res_table_BH_padj_thres_only <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_significant_0.05_degs_by_BH_padj_val_only_results.csv", sep = ''))
            
            mast_res_table_BH_padj_thres_0.15 <- read_csv(file= paste("Harmony/",t,"/",y,"/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_significant_0.05_degs_by_BH_padj_val_and_log2fc_0.15_results.csv", sep = ''))
            
            mast_res_table_BH_padj_thres_0.25 <- read_csv(file= paste("Harmony/",t,"/",y,"/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_significant_0.05_degs_by_BH_padj_val_and_log2fc_0.25_results.csv", sep = ''))
            
            #BC: Bonferroni Correction
            
            mast_res_table_BC_padj_thres_only <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_significant_0.05_degs_by_BC_padj_val_only_results.csv", sep = ''))
            
            mast_res_table_BC_padj_thres_0.15 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_significant_0.05_degs_by_BC_padj_val_and_log2fc_0.15_results.csv", sep = ''))
            
            mast_res_table_BC_padj_thres_0.25 <- read_csv(paste("Harmony/",t,"/",y,"/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_significant_0.05_degs_by_BC_padj_val_and_log2fc_0.25_results.csv", sep = ''))
          
          } #if MAST 
        
        ###############
        #For each combination of dea method and logfc threshold, plot volcano plot & bp for each ell type
      
        for (m in c("BH", "BC")){ #p adj method
          for (thresh in c("only", "0.15", "0.25")){ #avg_log2fc
            
            if (identical(thresh, "only")){
              t_num = 0
            }
            
            if (identical(thresh, "0.15")){
              t_num = 0.15
            }
            
            if (identical(thresh, "0.25")){
              t_num = 0.25
            }
            
            if (identical(d, "wilcox_stat_dea")){
              j = "wilcox"
            }
            
            if (identical(d, "seurat_mast_dea")){
              j = "mast"
            }
            
            #VOLCANO PLOT: scatterplot of p-value vs. log fold change of DEGs --> visualize significance of DEGs with large FC
            #colors, downregulated = blue, upregulared = red, no change = gray
            
            #Volcano for MAST, wilcox individual 1, wilcox individual 2
        
            #Volcano plot, still display all degs --> don't use edited tables
            volcano <- get(paste(j,"_res_table_",m,"_padj_thres_only", sep = '')) %>% 
                mutate(threshold = case_when(
                avg_log2FC > t_num ~ "Up",
                avg_log2FC < -t_num  ~ "Down",
                ((avg_log2FC < t_num) & (avg_log2FC > -t_num)) ~ "No Change"))
                
            if (j == "mast"){
                #Volcano
                #p-val (not adjusted) with log2fc > 0.15
                pdf(paste("Harmony/",t,"/combined/",x,"/",d,"/",x,"_volcano_plot_p_val_adj_",m,"_log2fc_",thresh,".pdf", sep = ""), width = 10, height = 10) 
                p<- ggplot() + 
                  geom_point(volcano, mapping = aes(x = avg_log2FC, y = -log10(.data[[paste0(m,"_p_val_adj",sep='')]]), colour = threshold)) +
                  ggtitle(paste("Volcano Plot of ",x," Significant APOE E44 vs. E33 DEGs with ",expression("Log"[2]*"(FC)")," > ",thresh," & ",m," P Val Adj < 0.05", sep = "")) +
                  xlab(expression("Log"[2]*"(FC)")) + 
                  ylab(expression("-Log"[10]*"(Adjusted P-value)")) +
                  scale_y_continuous(limits = c(0,50)) +
                  scale_colour_manual(values = c("Up"="red", "Down"="blue", "No Change"="lightgray")) +
                  theme(legend.position = "none",
                        plot.title = element_text(size = rel(1.5), hjust = 0.5),
                        axis.title = element_text(size = rel(1.25))) + 
                  theme_bw()
                print(p)
                dev.off()
            }
            
            if (j == "wilcox"){ #here
              #Volcano
              
              #individual 1
              pdf(paste("Harmony/",t,"/combined/",x,"/",d,"/",x,"individual_1_volcano_plot_p_val_adj_",m,"_log2fc_",thresh,".pdf", sep = ""), width = 10, height = 10) 
              p<- ggplot() + 
                geom_point(volcano, mapping = aes(x = avg_log2FC, y = -log10(ind_2_p_val_adj), colour = threshold)) +
                ggtitle(paste("Volcano Plot of Individual 1 ",x," Significant APOE E44 vs. E33 DEGs with ",expression("Log"[2]*"(FC)")," > ",thresh," & ",m," P Val Adj < 0.05", sep = "")) +
                xlab(expression("Log"[2]*"(FC)")) + 
                ylab(expression("-Log"[10]*"(Adjusted P-value)")) +
                scale_y_continuous(limits = c(0,50)) +
                scale_colour_manual(values = c("Up"="red", "Down"="blue", "No Change"="lightgray")) +
                theme(legend.position = "none",
                      plot.title = element_text(size = rel(1.5), hjust = 0.5),
                      axis.title = element_text(size = rel(1.25))) + 
                theme_bw()
              print(p)
              dev.off()
              
              #individual 2
              pdf(paste("Harmony/",t,"/combined/",x,"/",d,"/",x,"individual_2_volcano_plot_p_val_adj_",m,"_log2fc_",thresh,".pdf", sep = ""), width = 10, height = 10) 
              p<- ggplot() + 
                geom_point(volcano, mapping = aes(x = avg_log2FC, y = -log10(ind_1_p_val_adj), colour = threshold)) +
                ggtitle(paste("Volcano Plot of Individual 2 ",x," Significant APOE E44 vs. E33 DEGs with ",expression("Log"[2]*"(FC)")," > ",thresh," & ",m," P Val Adj < 0.05", sep = "")) +
                xlab(expression("Log"[2]*"(FC)")) + 
                ylab(expression("-Log"[10]*"(Adjusted P-value)")) +
                scale_y_continuous(limits = c(0,50)) +
                scale_colour_manual(values = c("Up"="red", "Down"="blue", "No Change"="lightgray")) +
                theme(legend.position = "none",
                      plot.title = element_text(size = rel(1.5), hjust = 0.5),
                      axis.title = element_text(size = rel(1.25))) + 
                theme_bw()
              print(p)
              dev.off()
              
            }
            
            print('done Volcano')
        
            #BARPLOT: number of up & down regulated genes in this cell type for this dea test
        
              #Get subsetted data table by p adj value threshold for this p val corr method
              deg_summary <- get(paste(j,"_res_table_",m,"_padj_thres_",thresh, sep = ''))  
              
              deg_summary <- deg_summary %>%
                mutate(deg_type = case_when(
                  avg_log2FC > 0 ~ "Up",
                  avg_log2FC < 0  ~ "Down"))
              
              deg_summary <- deg_summary %>%
                mutate(deg_count = case_when(
                  deg_type == "Up" ~ 1,
                  deg_type == 'Down' ~ 1)) 
              
              #save barplot
              pdf(paste("Harmony/",t,"/combined/",x,"/",d,"/",x,"_DEG_barplot_p_val_adj_",m,"_log2fc_",thresh,".pdf", sep = ""), width = 10, height = 6) 
              p<-ggplot(deg_summary, aes(x = factor(deg_type), y = deg_count)) + 
                geom_bar(stat = "identity") +ggtitle(paste("Number of Significant Up & Down Regulated DEGs in", x))+ 
                xlab('Direction of Expression in APOE 44') + ylab('Number of APOE DEGs')+ 
                theme_bw()
              print(p)
              dev.off()
              
              print('done bar')
              
          #####################
          #large barplot: number of total degs found by each method for each cell type
          
          DT_new_row <- data.table("dea_method" = d, 
                                   "p_adj_method" = method, 
                                   "log2FC_threshold" = thresh,
                                   "cell_type" = x, 
                                   "Up" = sum(deg_summary$deg_type == 'Up'), 
                                   "Down"= sum(deg_summary$deg_type == 'Down'),
                                   "degs" = sum(deg_summary$deg_count))
          
          DT <- rbindlist(list(DT, DT_new_row))
          
          ###
          #up and down degs found for each cell type by specific method and each logfc threshold
          
          meth_DT_new_row_1 <- data.table("cell_type" = x, 
                                          "p_adj_method" = method, 
                                          "log2FC_threshold" = thresh,
                                          "direction" = "Up",
                                          "degs" = sum(deg_summary$deg_type == 'Up'))
          
          meth_DT_new_row_2 <- data.table("cell_type" = x, 
                                          "p_adj_method" = method,
                                          "log2FC_threshold" = thresh,
                                          "direction" = "Down",
                                          "degs" = sum(deg_summary$deg_type == 'Down'))
          
          
          if (identical(d, "seurat_mast_dea")){ 
            mast_meth_DT <- rbindlist(list(mast_meth_DT, meth_DT_new_row_1))
            mast_meth_DT <- rbindlist(list(mast_meth_DT, meth_DT_new_row_2))
          }
          
          if (identical(d, "wilcoxon_stat_dea")){ 
            wilcox_meth_DT <- rbindlist(list(wilcox_meth_DT, meth_DT_new_row_1))
            wilcox_meth_DT <- rbindlist(list(wilcox_meth_DT, meth_DT_new_row_2))
          }
          
          bar_DT_new_row <- data.table("cell_type" = x, 
                                       "p_adj_method" = method,
                                       "log2FC_threshold" = thresh,
                                       "dea_method" = d,
                                       "Up" = sum(deg_summary$deg_type == 'Up'),
                                       "Down" = sum(deg_summary$deg_type == 'Down'))
          
          bar_DT <- rbindlist(list(bar_DT, bar_DT_new_row))
          
        } #vol & bar log2fc thresh
      } #vol & bar dea method'
    } #dea method
  }#cell type 
  
  #after data tables have been filled for each cell type for both methods:
  for (m in c("BH", "BC")){ #p adj method
    for (thresh in c("only", "0.15", "0.25")){ #avg_log2fc

      print(length(DT[0]))
      
      corr_name <- m
      log2FC_thres <- thresh
      
      DT <- as.data.frame(DT)
      wilcox_stat_meth_DT <- as.data.frame(wilcox_stat_meth_DT)
      seurat_mast_meth_DT <- as.data.frame(seurat_mast_meth_DT)
      
      print(length(DT[0]))
      
      #subset only for this correction method
      DT_x <- subset(DT, p_adj_method == corr_name)
      DT_x <- subset(DT, log2FC_threshold == log2FC_thres)
      
      print(length(DT[0]))
      
      wilcox_stat_meth_DT_x <- subset(wilcox_stat_meth_DT, p_adj_method == corr_name)
      wilcox_stat_meth_DT_x <- subset(wilcox_stat_meth_DT, log2FC_threshold == log2FC_thres)
      
      seurat_mast_meth_DT_x <- subset(seurat_mast_meth_DT, p_adj_method == corr_name)
      seurat_mast_meth_DT_x <- subset(seurat_mast_meth_DT, log2FC_threshold == log2FC_thres)
      
      #Barplot of all DEA Methods by cell type
      identities <- levels(Seurat$harmony_FCT)
      identities <- identities[identities %in% DT$cell_type]
      harmony_FCT_palette <- hue_pal()(length(identities))

      DT_x$cell_type <- factor(DT_x$cell_type, levels = identities)
      
      pdf(paste("Harmony/",t,"/combined/DEA_Method_Comparison_Barplot_",corr_name,"_",log2FC_thres,"_Log2FC.pdf", sep = ""), width = 12, height = 8)
      p <- ggplot(DT_x, aes(x= dea_method, y= degs, fill= cell_type)) +
        geom_bar(stat="identity", position="dodge") +
        scale_fill_manual(values = harmony_FCT_palette)+
        theme(axis.text.x = element_text(angle = 90))+ 
        theme_bw()
      print(p)
      dev.off()
      
      print("done DT bar")

      wilcox_stat_meth_DT_x$direction <- factor(wilcox_stat_meth_DT_x$direction, levels = c("Up", "Down"))
      
      pdf(paste("Harmony/",t,"/combined/Wilcox_Stat_DEA_Results_by_Cell_Type_Barplot_",corr_name,"_",log2FC_thres,"_Log2FC.pdf", sep = ""), width = 12, height = 6)
      p <- ggplot(wilcox_meth_DT_x, aes(x= cell_type, y= degs, fill= direction)) +
        geom_bar(stat="identity", position="dodge") +
        scale_fill_manual(values = c("Up" = "red","Down" = "blue")) +
        theme(axis.text.x = element_text(angle = 90))+ 
        theme_bw()
      print(p)
      dev.off()
      
      seurat_mast_meth_DT$direction <- factor(seurat_mast_meth_DT$direction, levels = c("Up", "Down"))
  
      pdf(paste("Harmony/",t,"/combined/Seurat_Mast_DEA_Results_by_Cell_Type_Barplot_",corr_name,"_",log2FC_thres,"_Log2FC.pdf", sep = ""), width = 12, height = 6)
      p <- ggplot(mast_meth_DT_x, aes(x= cell_type, y= degs, fill= direction)) +
        geom_bar(stat="identity", position="dodge") +
        scale_fill_manual(values = c("Up" = "red","Down" = "blue")) +
        theme(axis.text.x = element_text(angle = 90))+ 
        theme_bw()
      print(p)
      dev.off()
      
      print("done barplots")
    
      
      DT_x <- as.data.table(DT_x)
      wilcox_stat_meth_DT_x <- as.data.table(wilcox_stat_meth_DT_x)
      seurat_mast_meth_DT_x <- as.data.table(seurat_mast_meth_DT_x)
      
      library(plotrix)
      
      cell_name_list <- list()
      for (name in rownames(bar_DT)){
        if (name == 'Astrocyte'){
          cell_name_list <- append(cell_name_list, 'Ast')
        }
        if (name == 'Astrocyte 1'){
          cell_name_list <- append(cell_name_list, 'Ast 1')
        }
        if (name == 'Astrocyte 2'){
          cell_name_list <- append(cell_name_list, 'Ast 2')
        }
        if (name == 'Astrocyte 3'){
          cell_name_list <- append(cell_name_list, 'Ast 3')
        }
        if (name == 'Astrocyte 4'){
          cell_name_list <- append(cell_name_list, 'Ast 4')
        }
        if (name == 'Astrocyte 5'){
          cell_name_list <- append(cell_name_list, 'Ast 5')
        }
        if (name == 'Astrocyte 6'){
          cell_name_list <- append(cell_name_list, 'Ast 6')
        }
        if (name == 'NPC'){
          cell_name_list <- append(cell_name_list, 'NPC')
        }
        if (name == 'NPC 1 (Inhibitory Neuron)'){
          cell_name_list <- append(cell_name_list, 'NPC 1')
        }
        if (name == 'NPC 2 (Cycling NPC + GPC)'){
          cell_name_list <- append(cell_name_list, 'NPC 2')
        }
        if (name == 'NPC 3'){
          cell_name_list <- append(cell_name_list, 'NPC 3')
        }
        if (name == 'Excitatory Neuron'){
          cell_name_list <- append(cell_name_list, 'Ex')
        }
        if (name == 'Inhibitory Neuron'){
          cell_name_list <- append(cell_name_list, 'In')
        }
        if (name == 'Mature excitatory Neuron'){
          cell_name_list <- append(cell_name_list, 'Mat ex')
        }
        if (name == 'OPC + Oligodendrocyte'){
          cell_name_list <- append(cell_name_list, 'OPC/Oli')
        }
        if (name == 'VLMC'){
          cell_name_list <- append(cell_name_list, 'VLMC')
        }
        
      }
      
      #edit bar DT
      mast_bar_DT <- bar_DT[bar_DT$dea_method] == 'seurat_mast_dea'
      mast_bar_DT <- mast_bar_DT[(mast_bar_DT$p_adj_method == corr_name) & (mast_bar_DT$log2FC_threshold == log2FC_thres)]
      mast_bar_DT$p_adj_method <- NULL
      mast_bar_DT$log2FC_threshold <- NULL
      mast_bar_DT <- as.data.frame(mast_bar_DT)
      rownames(mast_bar_DT) <- mast_bar_DT$cell_type
      mast_bar_DT$cell_type <- NULL
      
      wilcox_bar_DT <- bar_DT[bar_DT$dea_method] == 'wilcox_stat_dea'
      wilcox_bar_DT <- wilcox_bar_DT[(wilcox_bar_DT$p_adj_method == corr_name) & (wilcox_bar_DT$log2FC_threshold == log2FC_thres)]
      wilcox_bar_DT$p_adj_method <- NULL
      wilcox_bar_DT$log2FC_threshold <- NULL
      wilcox_bar_DT <- as.data.frame(wilcox_bar_DT)
      rownames(wilcox_bar_DT) <- wilcox_bar_DT$cell_type
      wilcox_bar_DT$cell_type <- NULL
      
      pdf(paste("Harmony/",t,"/combined/Seurat_Mast_DEA_Results_by_Cell_Type_Barplot_",corr_name,"_",log2FC_thres,"_Log2FC.pdf", sep = ""), width = 10, height = 6)
      p <- color2D.matplot(mast_bar_DT, 
                           show.values = 0.5,
                           axes = FALSE,
                           main = "DEGs",
                           cex.main = 2,
                           xlab = "",
                           ylab = "",
                           vcex = 2,
                           vcol = "black",
                           extremes = c("white", "#009999"),
                           Hinton = TRUE)
      axis(1, at = seq_len(ncol(bar_DT)) - 0.5,
           labels = colnames(bar_DT), tick = FALSE, cex.axis = 2)
      axis(2, at = seq_len(nrow(bar_DT)) -0.5,
           labels = rev(cell_name_list), tick = FALSE, las = 1, cex.axis = 2)
      print(p)
      
      dev.off()
      
      pdf(paste("Harmony/",t,"/combined/Wilcox_Stat_DEA_Results_by_Cell_Type_Barplot_",corr_name,"_",log2FC_thres,"_Log2FC.pdf", sep = ""), width = 10, height = 6)
      p <- color2D.matplot(wilcox_bar_DT, 
                           show.values = 0.5,
                           axes = FALSE,
                           main = "DEGs",
                           cex.main = 2,
                           xlab = "",
                           ylab = "",
                           vcex = 2,
                           vcol = "black",
                           extremes = c("white", "#009999"),
                           Hinton = TRUE)
      axis(1, at = seq_len(ncol(bar_DT)) - 0.5,
           labels = colnames(bar_DT), tick = FALSE, cex.axis = 2)
      axis(2, at = seq_len(nrow(bar_DT)) -0.5,
           labels = rev(cell_name_list), tick = FALSE, las = 1, cex.axis = 2)
      print(p)
      
      dev.off()
      
      print("done all")
     
     }
   } #corr method

} #annotation type: Specific, General

library(VennDiagram)

harmony_cell_type_list <- as.vector(unique(Seurat$harmony_FCT))

p_val_DT <- data.table()
p_val_DT[, cell_type := ""]
p_val_DT[, gene := ""]
p_val_DT[, mast_p_value := 0]
p_val_DT[, individual_1_p_value := 0]
p_val_DT[, individual_2_p_value := 0]

for (x in harmony_cell_type_list){
  
  #check conditions for DEA to have been run
  
  Idents(Seurat) <- Seurat$individual
  individual_1 <-  subset(Seurat, idents = "individual_1")
  individual_2 <-  subset(Seurat, idents = "individual_2")
  
  Idents(individual_1) <- individual_1$harmony_FCT
  sub_1 <- subset(individual_1, idents = x) 
  
  Idents(individual_2) <- individual_2$harmony_FCT
  sub_2 <- subset(individual_2, idents = x) 
  
  apoe33_1 <- subset(sub_1, APOE_Genotype %in% "APOE 33")
  apoe44_1 <- subset(sub_1, APOE_Genotype %in% "APOE 44")
  
  apoe33_2 <- subset(sub_2, APOE_Genotype %in% "APOE 33")
  apoe44_2 <- subset(sub_2, APOE_Genotype %in% "APOE 44")
  
  if ((length(colnames(apoe33_1)) >= 5) & (length(colnames(apoe44_1)) >= 5) & (length(colnames(apoe44_2)) >= 5) & (length(colnames(apoe44_2)) >= 5)){
  
    for (m in list('BH', 'BC')){
      
      #all degs, not just significant
      # wilcox_1 <- read.csv(paste("Harmony/Specific/individual_1/",x,"/wilcox_stat_dea/harmony_",x,"_significant_0.05_degs_by_",m,"_padj_value_results.csv", sep = ""))
      # wilcox_2 <- read.csv(paste("Harmony/Specific/individual_2/",x,"/wilcox_stat_dea/harmony_",x,"_significant_0.05_degs_by_",m,"_padj_value_results.csv", sep = ""))
      # 
      wilcox_1 <- read.csv(paste("Harmony/Specific/individual_1/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
      wilcox_2 <- read.csv(paste("Harmony/Specific/individual_2/",x,"/wilcox_stat_dea/harmony_",x,"_wilcox_stat_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
      
      #get only those genes present in both lists
      com_degs <- c(unique(wilcox_1$gene_name), unique(wilcox_2$gene_name))
      com_degs <- as.data.table(com_degs)
      com_degs <- com_degs[duplicated(com_degs)]
      
      colnames(com_degs)<- "gene_name"
      com_degs[, statistic := 0]
      com_degs[, logFC := 0]
      com_degs[, ind_1_p_val := 0]
      com_degs[, ind_2_p_val := 0]
      com_degs[, ind_1_p_val_adj := 0]
      com_degs[, ind_2_p_val_adj := 0]
      
      if (length(com_degs != 0)){
      
      for (gene_test in com_degs$gene_name){
        
        stat_1 <- (subset(wilcox_1, gene_name == gene_test))$statistic #u-stat
        stat_2 <- (subset(wilcox_2, gene_name == gene_test))$statistic #u-stat
        logFC_1 <- (subset(wilcox_1, gene_name == gene_test))$logFC #avg_log2FC
        logFC_2 <- (subset(wilcox_2, gene_name == gene_test))$logFC #avg_log2FC
        
        #the gene should be regulated in the same direction in both individuals
        if (sign(logFC_1) == sign(logFC_2)){
          com_degs[com_degs$gene_name == gene_test]$statistic <- mean(stat_1, stat_2) #mean(stat_1, stat_2) #take the avg of the stats
          com_degs[com_degs$gene_name == gene_test]$logFC <- mean(logFC_1, logFC_2) #take the avg of the log fold change
          com_degs[com_degs$gene_name == gene_test]$ind_1_p_val <- (subset(wilcox_1, gene_name == gene_test))$p_val #ind 1 p_val
          com_degs[com_degs$gene_name == gene_test]$ind_2_p_val <- (subset(wilcox_2, gene_name == gene_test))$p_val #ind 2 p_val
          
          p_val_adj = paste(m, "_p_val_adj", sep = '')
          
          com_degs[com_degs$gene_name == gene_test]$ind_1_p_val_adj <- (subset(wilcox_1, gene_name == gene_test))$p_val_adj #ind 1 p_val_adj
          com_degs[com_degs$gene_name == gene_test]$ind_2_p_val_adj <- (subset(wilcox_2, gene_name == gene_test))$p_val_adj #ind 2 p_val_adj
        }
        else{ #else, remove the gene from the combined list
          com_degs <- com_degs[gene_name != gene_test]
        }
      }
      
      write.csv(as.data.frame(com_degs), file= paste("Harmony/Specific/combined/",x,"/harmony_",x,"_combined_wilcox_stat_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
      
      print("done csv 1")
      
      #FDR < 0.05
      res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
      
      
      #p-val (not adjusted)
      
      res_table_pval_thres <- res_tbl %>% 
        mutate(threshold = p_val < 0.05) #significant DEGs
      
      res_table_pval_thres <- res_table_pval_thres[res_table_pval_thres$p_val < 0.05,] #only look at significant DEGs based on p-value < 0.05
      res_table_pval_thres <- res_table_pval_thres[!is.na(res_table_pval_thres$p_val),] 
      res_table_pval_thres <- as.data.frame(res_table_pval_thres)
      
      
      #Venn diagram of overlapping DEGs between individual_1 & individual_2
      
      ind_1_ex <- wilcox_1$gene_name
      ind_2_ex <- wilcox_2$gene_name
      
      ind_1_ex_up <- subset(wilcox_1, logFC > 0)$gene_name
      ind_2_ex_up <- subset(wilcox_2,logFC > 0)$gene_name
  
      ind_1_ex_down <- subset(wilcox_1, logFC < 0)$gene_name
      ind_2_ex_down <- subset(wilcox_2,logFC < 0)$gene_name
      
      venn.diagram(
        x = list(ind_1_ex, ind_2_ex),
        category.names = c("Individual 1", "Individual 2"),
        print.mode=c("raw","percent"),
        lwd = 1,
        cex = 2,
        cat.cex = 2,
        main.cex = 2,
        cat.pos = 0,
        main = paste("Overlapping ",x, " DEGs Between Individuals", sep = ''),
        filename = paste('Harmony/Specific/combined/',x,'/individual_shared_DEG_venn.png', sep = ''),
        height = 5000, 
        width = 5000) 
      
      venn.diagram(
        x = list(ind_1_ex_up, ind_2_ex_up),
        category.names = c("Individual 1" , "Individual 2"),
        print.mode=c("raw","percent"),
        lwd = 1,
        cex = 2,
        cat.cex = 2,
        main.cex = 2,
        cat.pos = 0,
        main = paste("Overlapping ",x, " Upregulated DEGs Between Individuals", sep = ''),
        filename = paste('Harmony/Specific/combined/',x,'/individual_shared_up_DEG_venn.png', sep = ''),
        height = 5000, 
        width = 5000) 
      
      venn.diagram(
        x = list(ind_1_ex_down, ind_2_ex_down),
        category.names = c("Individual 1" , "Individual 2"),
        print.mode=c("raw","percent"),
        lwd = 1,
        cex = 2,
        cat.cex = 2,
        main.cex = 2,
        cat.pos = 0,
        main = paste("Overlapping ",x," Downregulated DEGs Between Individuals", sep = ''),
        filename = paste('Harmony/Specific/combined/',x,'/individual_shared_down_DEG_venn.png',  sep = ''),
        height = 5000, 
        width = 5000) 
      
      print("done venn 1")

      #Venn diagram of overlapping DEGs between combined wilcox & MAST, split by up & down regulated
      mast <- read.csv(file = paste("Harmony/Specific/All/",x,"/seurat_mast_dea/harmony_",x,"_seurat_mast_dea_Results_GTFAnnotated_NoGeneIDDuplicates.csv", sep = ''))
      
      mast <- subset(mast, p_val_adj < 0.05)
      
      #get common genes
      meth_com_degs <- c(unique(mast$gene_name), unique(com_degs$gene_name))
      meth_com_degs <- meth_com_degs[duplicated(meth_com_degs)]
      meth_com_degs <- as.data.table(meth_com_degs)
      
      colnames(meth_com_degs)<- "X"
      meth_com_degs[, statistic := 0]
      meth_com_degs[, p_val_adj_1 := 0]
      meth_com_degs[, p_val_adj_2 := 0]
      meth_com_degs[, p_val_1 := 0]
      meth_com_degs[, p_val_2 := 0]
      meth_com_degs[, mast_p_val_adj := 0]
      meth_com_degs[, wilcox_logFC := 0] #logFC
      meth_com_degs[, mast_logFC := 0] #avg_log2FC
      
      for (gene_test in meth_com_degs$gene_name){
        
        statistic <- subset(com_degs, gene_name == gene_test)$statistic #avg stat btw ind 1 & ind 2
        p_val_adj_1 <- subset(com_degs, gene_name == gene_test)$ind_1_p_val_adj #ind 1 p_val_adj
        p_val_adj_2 <- subset(com_degs, gene_name == gene_test)$ind_2_p_val_adj #ind 2 p_val_adj
        p_val_1 <- subset(com_degs, gene_name == gene_test)$ind_1_p_val_adj #ind 1 p_val
        p_val_2 <- subset(com_degs, gene_name == gene_test)$ind_2_p_val_adj #ind 2 p_val
        mast_p_val_adj <- subset(mast, gene_name == gene_test)$p_val_adj #mast p_val_adj
        wilcox_logFC <- subset(com_degs, gene_name == gene_test)$logFC #wilcox logFC
        mast_logFC <- subset(mast, gene_name == gene_test)$avg_log2FC #mast logFC
        
        if (sign(stat_1) == sign(stat_2)){
          meth_com_degs[meth_com_degs$gene_name == gene_test]$statistic <- statistic #take the avg of the stats
          meth_com_degs[meth_com_degs$gene_name == gene_test]$p_val_adj_1 <- p_val_adj_1 #take the avg of the log fold change
          meth_com_degs[meth_com_degs$gene_name == gene_test]$p_val_adj_2 <- p_val_adj_2 #ind 1 p_val
          meth_com_degs[meth_com_degs$gene_name == gene_test]$p_val_1 <- p_val_1 #ind 1 p_val
          meth_com_degs[meth_com_degs$gene_name == gene_test]$p_val_2 <- p_val_2 #ind 1 p_val
          meth_com_degs[meth_com_degs$gene_name == gene_test]$mast_p_val_adj <- mast_p_val_adj  #ind 2 p_val
          meth_com_degs[meth_com_degs$gene_name == gene_test]$wilcox_logFC <- wilcox_logFC #wilcox logFC
          meth_com_degs[meth_com_degs$gene_name == gene_test]$mast_logFC <- mast_logFC #mast logFC
        }
        else{ #else, remove the gene from the combined list
          meth_com_degs$gene_test <- NULL
        }
        
      }
      
      print("check 1")
      
      mast_ex <- mast$gene_name
      wilcox_ex <- com_degs$gene_name
  
      
      mast_ex_up <- subset(mast, avg_log2FC > 0)$gene_name
      wilcox_ex_up <- subset(com_degs, logFC > 0)$gene_name
      
      mast_ex_down <- subset(mast, avg_log2FC < 0)$gene_name
      wilcox_ex_down <- subset(com_degs, logFC < 0)$gene_name
      
      venn.diagram(
        x = list(wilcox_ex, mast_ex),
        category.names = c("Wilcoxon" , "MAST"),
        print.mode=c("raw","percent"),
        lwd = 1,
        cex = 2,
        cat.cex = 2,
        main.cex = 2,
        cat.pos = 0,
        main = paste("Overlapping ",x, " DEGs Across DEA Methods", sep = ''),
        filename = paste('Harmony/Specific/combined/',x,'/dea_method_shared_DEG_venn.png', sep = ''),
        height = 5000, 
        width = 5000) 
      
      venn.diagram(
        x = list(wilcox_ex_up, mast_ex_up),
        category.names = c("Wilcoxon" , "MAST"),
        print.mode=c("raw","percent"),
        lwd = 1,
        cex = 2,
        cat.cex = 2,
        main.cex = 2,
        cat.pos = 0,
        main = paste("Overlapping ",x, " Upregulated DEGs Across DEA Methods", sep = ''),
        filename = paste('Harmony/Specific/combined/',x,'/dea_method_shared_up_DEG_venn.png', sep = ''),
        height = 5000, 
        width = 5000) 
      
      venn.diagram(
        x = list(wilcox_ex_down, mast_ex_down),
        category.names = c("Wilcoxon" , "MAST"),
        print.mode=c("raw","percent"),
        lwd = 1,
        cex = 2,
        cat.cex = 2,
        main.cex = 2,
        cat.pos = 0,
        main = paste("Overlapping ",x, " Downregulated DEGs Across DEA Methods", sep = ''),
        filename = paste('Harmony/Specific/combined/',x,'/dea_method_shared_down_DEG_venn.png', sep = ''),
        height = 5000, 
        width = 5000)  
      
      print("done venn 2")
      
      meth_com_degs <- meth_com_degs %>% slice_max(n= 10000000, order_by = abs(wilcox_logFC))
      
      mast_p =  meth_com_degs$mast_p_val
      ind_1_p =  meth_com_degs$p_val_1
      ind_2_p =   meth_com_degs$p_val_2
      
      if (length(meth_com_degs$gene_name) > 0){
        pdf(paste("Harmony/Specific/combined/",x,"/ind1_v_ind2_p_value_scatterplot.pdf", sep = "")) 
        
        # creating scatterplot 
        data <- data.table()
        data$individual_1_p_value <- -log10(ind_1_p)
        data$individual_2_p_value <- -log10(ind_2_p)
        data <- as.data.frame(data)
        
        p <- ggplot(data) + 
          geom_point(aes(x= individual_1_p_value, y= individual_2_p_value))+
          xlab("individual 1 -log10(p value)")+
          ylab("individual 2 -log10(p value)")+
          xlim(0, 50)+
          ylim(0, 50)
        
        print(p)
        
        dev.off()
        
        print("done scatter 1")
        
        #Barplot for top 50 DEGs in celltype: wilcox avg log fold change, MAST log fold change
        
        pdf(paste("Harmony/Specific/combined/",x,"/ind1_v_mast_p_value_scatterplot.pdf", sep = "")) 
      
        data <- data.table()
        data$individual_1_p_value <- -log10(ind_1_p)
        data$mast_p_value <- -log10(mast_p)
        data <- as.data.frame(data)
        
        p <- ggplot(data) + 
          geom_point(aes(x= individual_1_p_value, y= mast_p_value))+
          xlab("individual 1 -log10(p value)")+
          ylab("mast -log10(p value)")+
          xlim(0, 50)+
          ylim(0, 50)
        
        print(p)
        
        dev.off()
        print("done scatter 2")
        
        pdf(paste("Harmony/Specific/combined/",x,"/ind2_v_mast_p_value_scatterplot.pdf", sep = ""))
        
        data <- data.table()
        data$individual_2_p_value <- -log10(ind_2_p)
        data$mast_p_value <- -log10(mast_p)
        data <- as.data.frame(data)
        
        p <- ggplot(data) + 
          geom_point(aes(x= individual_2_p_value, y= mast_p_value))+
          xlab("individual 2 -log10(p value)")+
          ylab("mast -log10(p value)")+
          xlim(0, 50)+
          ylim(0, 50)
        
        plot(p)
        
        dev.off()
        print("done scatter 3")
        
          for (gene_test in meth_com_degs$gene_name){
              a<- subset(meth_com_degs, gene_name == gene_test)$mast_p_val
              b<- subset(meth_com_degs, gene_name == gene_test)$p_val_1
              c<- subset(meth_com_degs, gene_name == gene_test)$p_val_2
            
              new_row <- data.table("cell_type" = x, "gene" = gene_test, "mast_p_value" = -log10(a), "individual_1_p_value" = -log10(b), "individual_2_p_value" = -log10(c))
              p_val_DT <- rbindlist(list(p_val_DT, new_row))
            } #for 
          } #if
      
        } #if deg list length
      
      } #p adj method: BH or BC
  } #condition for dea
} #cell type

#p value correlation table

p_val_DT <- as.data.frame(p_val_DT)[-1,]

pdf(paste("Harmony/Specific/combined/giant_p_value_scatterplot.pdf", sep = ""))
  
  p1 <- ggplot(p_val_DT)+
        geom_point(aes(x = individual_1_p_value, y = individual_2_p_value))+ 
        facet_wrap(vars(cell_type))
  
  p2 <- ggplot(p_val_DT, aes(x = individual_1_p_value, y = mast_p_value))+
    geom_point()
  
  p2 <- p2 + facet_wrap(vars(cell_type))
  
  p3 <- ggplot(p_val_DT, aes(x = individual_2_p_value, y = mast_p_value))+
    geom_point()
  
  p3 <- p3 + facet_wrap(vars(cell_type))
   
  giant_p <- p1 + p2 + p3 + plot_layout(widths = c(50, 50, 50)) #+ plot_layout(heights = c(50, 50, 50))
  
  #dev.new(width=50, height=50)
  giant_p

dev.off()
