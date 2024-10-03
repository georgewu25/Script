library(WGCNA)

res_enr_final <- read.csv("/output/WGCNA/deepsplit2/functional_enrichment_padj0.1_overlap5.csv")

module_yellow <- res_enr_final[res_enr_final$query. == "yellow", ]

module_yellow_sorted <- module_yellow[order(module_yellow$padj, decreasing = F), ]

module_leading_edges <- unlist(lapply(module_yellow_sorted$genes.overlap, function(x) {
  unlist(strsplit(x, "\\|"))
}))



expression_df <- read.csv("/output/WGCNA/normalized_expression.csv")

chooseTopHubInEachModule(expression_df[,-1], module_color)

connectivity_allClusters <- intramodularConnectivity.fromExpr(expression_df[,-1], module_color, 
                                                              corFnc = "bicor", corOptions = "use = 'p'",
                                                              weights = NULL,
                                                              distFnc = "dist", distOptions = "method = 'euclidean'",
                                                              networkType = "signed", power = 9,
                                                              scaleByMax = FALSE,
                                                              ignoreColors = if (is.numeric(colors)) 0 else "grey",
                                                              getWholeNetworkConnectivity = TRUE)

rownames(connectivity_allClusters) <- colnames(expression_df[,-1])
connectivity_allClusters$module_color <- module_color
connect <- connectivity_allClusters[order(connectivity_allClusters$module_color,-connectivity_allClusters$kWithin),]
