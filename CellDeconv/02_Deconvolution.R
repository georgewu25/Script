library("Seurat")
library(BayesPrism)

runBayesPRISM <- function(){
  #Loads in sc and bulk data, initializes cell type and state vectors, filters for uninformative genes,
  #Creates bulk data matrix necessary for creating PRISM object (sample x gene matrix)
  bk.data = t(pseudobulkMix)
  
  #Creates sc data matrix necessary for creating prism object (sc x gene matrix)
  sc.data = t(as.matrix(inputSeurat@assays$RNA@counts))
  
  #Initializes cell.type labels (lower resolution of annotation)
  cell.type.labels = as.character(inputSeurat@meta.data[[cellTypeSeurat]])
  
  #Initializes cell.state as higher resolution annotation
  cell.state.labels = as.character(inputSeurat@meta.data[[cellTypeSeurat]])
  
  #Filter outlier genes that are highly expressed but noninformative
  sc.data.filtered = cleanup.genes(input = sc.data,
                                   input.type = "count.matrix",
                                   species = "hs",
                                   gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "act", "hb", "MALAT1"),
                                   exp.cells=5)
  #Filters for protein coding genes
  sc.data.filtered.pc = select.gene.type(sc.data.filtered, gene.type = "protein_coding")
  
  
  # Initailizes PRISM object ##
  myPrism = new.prism(
    reference = sc.data.filtered.pc,
    mixture = bk.data,
    input.type = "count.matrix",
    cell.type.labels = cell.type.labels,
    cell.state.labels = cell.state.labels,
    key = NULL,
    outlier.cut = 0.01,
    outlier.fraction = 0.1
  )
  
  bp.res <- run.prism(prism = myPrism, n.cores=16)
  resultsPRISM = get.fraction(bp=bp.res, which.theta="first", state.or.type="type")
  resultsPRISM = t(resultsPRISM); resultsPRISM = data.frame(resultsPRISM)
  finalResults = resultstodf(resultsPRISM, trueProportionList, deconvolution)
  write.csv(finalResults, outputFile, row.names = FALSE)
}


######CIBERSORT
CoreAlg <- function(X, y){
  
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #Main function
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #Normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

#Permutations
doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}


#main function
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){
  
  X <- sig_matrix
  Y <- mixture_file
  
  P <- perm #number of permutations
  
  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}
  
  #quantile normalization
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #Intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #Standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #Empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}
  
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures
  while(itor <= mixtures){
    
    y <- Y[,itor]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    result <- CoreAlg(X, y)
    
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    #Calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}
    
    itor <- itor + 1
    
  }
  
  write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}

runCIBERSORT <- function(){
  sharedGenes = intersect(rownames(referenceProfile), rownames(pseudobulkMix))
  referenceProfile = referenceProfile[sharedGenes,]
  pseudobulkMix = pseudobulkMix[sharedGenes,]
  CIBERSORTresults <<- CIBERSORT(referenceProfile, pseudobulkMix)
  
  CIBERSORTresults = CIBERSORTresults[, 1:(ncol(CIBERSORTresults)-3)]
  CIBERSORTresults = data.frame(CIBERSORTresults)
  CIBERSORTresults = t(CIBERSORTresults)
  resultsdf = resultstodf(data.frame(CIBERSORTresults),trueProportionList,"CIBERSORT")
  resultsdf$CellType = vapply(resultsdf$CellType, function(x) substr(x, 1, nchar(x)), FUN.VALUE = character(length = 1L))
  write.csv(resultsdf, outputFile, row.names = FALSE)
}


runDeTREM <- function(){
  sc = inputSeurat
  sc_expr = as.matrix(sc@assays$RNA@counts)
  
  #Creates Annotated Dataframe of input's cellType and sampID
  celltype_column=cellTypeSeurat
  sample_column=sample_metadata
  sc_meta=sc@meta.data[,c(celltype_column,sample_column)]
  sc_meta=AnnotatedDataFrame(sc_meta)
  colnames(sc_meta)=c("Celltype","sampID")
  
  #Creates Expression set
  sc.eset = ExpressionSet(assayData = sc_expr, phenoData = sc_meta)
  bulk.eset = ExpressionSet(assayData = pseudobulkMix)
  estimates=DeTREM(bulk.eset = bulk.eset, sc.eset = sc.eset, clusters = 'Celltype',samples = 'sampID', verbose = F)
  
  #Estimated Proportions
  estimatedProportions=t(estimates$Est.prop.weighted)
  write.csv(resultstodf(data.frame(estimatedProportions), trueProportionList, "DeTREM"), outputFile, row.names = FALSE)
}


runDSA <- function() {
  bulkExpr = ExpressionSet(pseudobulkMix)
  
  #Converts gene names
  geneList <- rownames(exprs(bulkExpr))
  
  if (substring(geneList[1], 1, 4) != "ENSG"){
    ensembl_id <- TransSymboltoEnsemblVers(geneList)
    
    for (i in 1:length(geneList)){
      if (geneList[i] %in% ensembl_id$hgnc_symbol){
        geneList[i] = ensembl_id$ensembl_gene_id_version[ensembl_id$hgnc_symbol == geneList[i]][1]
      }
    }
    fData(bulkExpr) <- data.frame(Geneid = geneList, Symbol = rownames(exprs(bulkExpr)))
    rownames(bulkExpr) <- geneList
  } else {
    symbolVec <- c()
    symbol_id <- TransEnsembltoSymbol(rownames(bulkExpr)); symbol_id = data.table(symbol_id)
    symbol_id = removeDuplicateSymbols(symbol_id); symbol_id = removeDuplicateEnsembl(symbol_id)
    bulkExpr = bulkExpr[rownames(bulkExpr) %in% symbol_id$ensembl_gene_id]
    
    for (x in 1:nrow(bulkExpr)){
      symbolVec <- c(symbolVec, symbol_id[symbol_id$ensembl_gene_id[x]== rownames(bulkExpr)]$hgnc_symbol)
    }
    
    fData(bulkExpr) <- data.frame(Geneid = rownames(exprs(bulkExpr)), Symbol = symbolVec)
  }
  
  
  exprs(bulkExpr) = exprs(bulkExpr) + 0.00001
  genes = fData(bulkExpr)
  
  #Marker list
  markers = BRETIGEA::markers_df_human_brain #Character vector mapping marker gene name to cell type
  markers = markers[markers$cell %in% unique(celltype),]
  markers = markers[markers$markers %in% genes$Symbol, ] #Keeps marker genes that are present in sample
  
  
  markers$Geneid = genes$Geneid[match(markers$markers, genes$Symbol)] #Convert marker gene symbosl to Ensembl ids.
  
  markers = markers[markers$Geneid %in% rownames(bulkExpr), ] #Keeps the marker rows where the gene id and symbol matches with the bulkExpr
  
  markerSubset = lapply(split(markers$Geneid, markers$cell)[unique(celltype)], function(x, n) x[1:min(n, length(x))], n=5) #Obtains subset of gene markers 
  a = unlist(markerSubset); a = a[duplicated(a)]
  if(length(a) > 0) markerSubset = lapply(markerSubset, function(x, a) x[!x %in% a], a = a)#Removes duplicates
  m1 = MarkerList(markerSubset)
  
  #Runs DSA
  res2 <- ged(exprs(bulkExpr), m1, method = 'DSA')
  
  write.csv(resultstodf(data.frame(coef(res2)),trueProportionList,"Minghui"), outputFile, row.names = FALSE)
}


runDWLS <- function(){
  #Prepares reference profiles for DWLS by creating a matrix of single-cell expression and cell-type vector
  scdata <- as.matrix(inputSeurat@assays[["RNA"]]@counts)
  
  referenceDWLSMatrix = referenceProfile
  
  #Generates list of trimmed data (bulk and sc expression values) for each sample 
  trimmatrix <- list()
  for (i in 1:ncol(pseudobulkMix)){
    trimmatrix[[i]] = trimData(referenceDWLSMatrix,pseudobulkMix[,i])
  }
  names(trimmatrix) <- colnames(pseudobulkMix)
  
  
  results <<- list()
  
  for (i in unique(names(trimmatrix))){
    #Creates list entry for results
    prop = solveDampenedWLS(trimmatrix[[i]]$sig,trimmatrix[[i]]$bulk)
    results[[i]] <<- prop[sort(names(prop))]
  }
  
  results <<- data.frame(results)
  write.csv(resultstodf(results,trueProportionList,"DWLS"), outputFile, row.names = FALSE)
}


runMUSIC <- function() {
  #Subsamples sc reference if too large
  if (ncol(inputSeurat) > 60000){
    inputSeurat = inputSeurat[,sample(x = 1:ncol(inputSeurat),size = 60000)]
  }

  #Develops expression object for SCDC
  counts <- as.matrix(inputSeurat@assays[["RNA"]]@counts); #querycounts <- as.matrix(querySeurat@assays[["RNA"]]@counts)
  inputeset <- ExpressionSet(counts, fdata = rownames(counts)); #queryeset <- ExpressionSet(querycounts, fdata = rownames(querycounts))
  
  inputeset@phenoData@data[["cluster"]] = as.character(inputSeurat@meta.data[[cellTypeSeurat]])
  inputeset@phenoData@data[["sample"]] = as.character(inputSeurat@meta.data[[sample_metadata]])
  inputeset@phenoData@varMetadata = rbind(data.frame("labelDescription" = NA),data.frame("labelDescription" = NA))

  #Runs MUSIC 
  inputeset.sce = SingleCellExperiment(list(counts=as.matrix(inputSeurat@assays$RNA@counts)), colData = inputSeurat@meta.data)
  
  resultsMUSIC <<- music_prop(bulk.mtx = pseudobulkMix, sc.sce = inputeset.sce, clusters = cellTypeSeurat, samples = sample_metadata, select.ct = levels(inputSeurat@meta.data[[cellTypeSeurat]]))
  
  write.csv(resultstodf(data.frame(t(resultsMUSIC[["Est.prop.weighted"]])), kprop = trueProportionList,"MUSIC"),outputFile,row.names = FALSE)
}


runRLR <- function(){
  #Runs RLR
  proportionsRLR <<- do.call(cbind.data.frame,lapply(apply(pseudobulkMix,2,function(x) MASS::rlm(x ~ referenceProfile, maxit=1000)), function(y) y$coefficients[-1]))
  proportionsRLR <<- apply(proportionsRLR, 2, function(x) ifelse(x < 0, 0, x))
  proportionsRLR <<- apply(proportionsRLR, 2, function(x) x/sum(x))
  rownames(proportionsRLR) <<- unname(vapply(rownames(proportionsRLR), function(x) substring(x, 17, nchar(x)), FUN.VALUE = character(length = 1L)))
  
  write.csv(resultstodf(data.frame(proportionsRLR),trueProportionList,"RLR"), outputFile, row.names = FALSE)
}


runSCDC <- function() {
  #Subsamples sc reference if too large
  if (ncol(inputSeurat) > 60000){
    inputSeurat = inputSeurat[,sample(x = 1:ncol(inputSeurat),size = 60000)]
  }
  #Develops expression object for SCDC
  counts <- as.matrix(inputSeurat@assays[["RNA"]]@counts); #querycounts <- as.matrix(querySeurat@assays[["RNA"]]@counts)
  inputeset <- ExpressionSet(counts, fdata = rownames(counts)); #queryeset <- ExpressionSet(querycounts, fdata = rownames(querycounts))
  
  inputeset@phenoData@data[["cluster"]] = as.character(inputSeurat@meta.data[[cellTypeSeurat]])
  inputeset@phenoData@data[["sample"]] = as.character(inputSeurat@meta.data[[sample_metadata]])
  inputeset@phenoData@varMetadata = rbind(data.frame("labelDescription" = NA),data.frame("labelDescription" = NA))
  
  
  #Quality controlled expression set
  inputset.qc <- SCDC_qc_ONE(inputeset,ct.varname = "cluster", sample = "sample", scsetname = 'reference', ct.sub = levels(inputSeurat@meta.data[[cellTypeSeurat]]))

  #Runs SCDC deconvolution method
  resultsSCDC <<- SCDC_prop_ONE(bulk.eset = pseudobulkMix, sc.eset = inputset.qc$sc.eset.qc, 
                            ct.varname = "cluster", sample = "sample",
                            ct.sub = levels(inputSeurat@meta.data[[cellTypeSeurat]]))
  
  write.csv(resultstodf(data.frame(t(resultsSCDC[["prop.est.mvw"]])), kprop = trueProportionList,"SCDC"),outputFile,row.names = FALSE)
}

