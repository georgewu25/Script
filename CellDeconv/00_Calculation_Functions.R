library(biomaRt)
library(data.table)

GetMartGenes<-function()biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

TransSymboltoEnsemblVers<-function(hgnc_symbols){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'hgnc_symbol'),
                                               filters = 'hgnc_symbol', 
                                               values = hgnc_symbols, 
                                               mart = GetMartGenes())))
}

TransEnsembltoSymbol<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                                               filters = 'ensembl_gene_id', 
                                               values = ensembl_ids, 
                                               mart = GetMartGenes())))
}

convertFormat <- function(readDF){
  readDF = data.frame(readDF)
  rownames(readDF) = readDF[[1]]
  readDF = readDF[,2:ncol(readDF)]
  readDF
}

updateTrueProportions <- function(resultsPath,truePropPath){
  trueProportionList = convertFormat(read.csv(truePropPath))
  resultsDF = data.frame(read.csv(resultsPath))
  
  for (r in 1:nrow(resultsDF)){
    sample = resultsDF[r,]$Sample
    ct = resultsDF[r,]$CellType
    trueProp = trueProportionList[ct,sample]
    resultsDF[r,]$KnownProportions = trueProp
  }
  resultsDF
}

redoTruePropFolder <- function(folderDir, truePropPath){
  outputFiles = list.files(folderDir)
  outputFiles_fulldir = vapply(outputFiles, function(x) paste(folderDir, x, sep = ""), FUN.VALUE = character(length=1L))
  
  for (dir in outputFiles_fulldir){
    adjustedDF = updateTrueProportions(dir, truePropPath)
    write.csv(adjustedDF, dir, row.names = FALSE)
  }
}

convertTruePropFormat <- function(df){
  df = data.table(df)

  samples = unique(df$individualID)
  cellTypes = unique(df$cell_type)
  truePropDF = matrix(nrow = length(cellTypes), ncol = length(samples))
  
  colnames(truePropDF) = samples
  rownames(truePropDF) = cellTypes
  
  for (sample in colnames(truePropDF)){
    for (ct in rownames(truePropDF)){
      trueProp = df[individualID == sample][cell_type == ct]
      if (nrow(trueProp)>0){
        truePropDF[ct,sample] = trueProp$n.cells
      } else {
        truePropDF[ct,sample] = 0
      }
    }
  }
  truePropDF = data.frame(truePropDF)
  apply(truePropDF,2,function(x) x/sum(x))
}

convertSeuratGeneNamestoID <- function(seuratObj, featureTranslationDF) {
  library(Seurat)
  library(data.table)
  #converts gene names to ENSG
  featureTranslationDF = featureTranslationDF[gene_name %nin% featureTranslationDF$gene_name[duplicated(featureTranslationDF$gene_name)]]
  countsMatrix = as.matrix(seuratObj@assays$RNA@counts); countsMatrix = countsMatrix[rownames(countsMatrix) %in% featureTranslationDF$gene_name,]
  rownames(countsMatrix) = unname(vapply(rownames(countsMatrix), function(x) featureTranslationDF[gene_name == x]$gene_id, character(length=1L)))
  
  #Creates new seurat object
  newSeuratObj = CreateSeuratObject(counts = countsMatrix)
  newSeuratObj = NormalizeData(newSeuratObj)
  newSeuratObj@meta.data = seuratObj@meta.data
  
  newSeuratObj
}
removeDuplicateSymbols <- function(ensemblsymboldf){
  symbolsVec <- c()
  duplicatedSymbols = unique(ensemblsymboldf$hgnc_symbol[duplicated(ensemblsymboldf$hgnc_symbol)])
  for (x in ensemblsymboldf$hgnc_symbol){
    if (!(x %in% duplicatedSymbols)){
      symbolsVec <- c(symbolsVec, x)
    }
  }
  ensemblsymboldf = data.table(ensemblsymboldf)
  ensemblsymboldf[hgnc_symbol %in% symbolsVec]
}

removeDuplicateEnsembl <- function(ensemblsymboldf){
  symbolsVec <- c()
  duplicatedSymbols = unique(ensemblsymboldf$ensembl_gene_id[duplicated(ensemblsymboldf$ensembl_gene_id)])
  for (x in ensemblsymboldf$ensembl_gene_id){
    if (!(x %in% duplicatedSymbols)){
      symbolsVec <- c(symbolsVec, x)
    }
  }
  ensemblsymboldf = data.table(ensemblsymboldf)
  ensemblsymboldf[ensembl_gene_id %in% symbolsVec]
}


##########


# result = Data frame with id as column and cell as row
# kprop (known proportion) = Data.frame with id as column and cell as row
# method = String of method
resultstodf <- function(result, kprop, method){
  
  cellTypevec <- c()
  samplevec <- c()
  knownproportions <- c()
  estimatedproportions <- c()
  methodvec <- c()
  
  #Iterates across sample vectors
  for (i in 1:length(result)){
    
    #Iterates within sample vector
    for (k in 1:length(result[[i]])){
     
      celltype = rownames(result)[k]
      sample = names(result)[i]
      
      cellTypevec <- append(cellTypevec, celltype)
      samplevec <- append(samplevec, sample)
      estimatedproportions <- append(estimatedproportions, result[celltype,sample])
      knownproportions <- append(knownproportions, kprop[celltype,sample])
      methodvec <- append(methodvec, method)
    }
  }
  
  data.frame(CellType = cellTypevec, Sample = samplevec, KnownProportions = knownproportions, Estimatedproportions = estimatedproportions, Method = methodvec)
  
}

resultstoMatrix <- function(filePath){
  resultsDF = data.table(read_csv(filePath))
  samples = unique(resultsDF$Sample); ct = unique(resultsDF$CellType)
  
  resultsMatrix = matrix(nrow = length(ct), ncol = length(samples))
  rownames(resultsMatrix) = ct; colnames(resultsMatrix) = samples
  
  for (s in samples){
    for (c in ct){
      resultsMatrix[c,s] = resultsDF[CellType == c][Sample == s]$Estimatedproportions
    }
  }
  return(resultsMatrix)
}

addMetaDataSinai <- function(resultsDF, metadata){
  resultsDF = data.table(resultsDF)
  metadata = data.table(metadata)
  metadata = metadata[ID %in% unique(resultsDF$Sample),]
  resultsDF[["Region"]] = metadata$Region
  resultsDF[["CERAD"]] = metadata$CERAD
  
  return(resultsDF)
}



runSeuratPCApipeline <- function(seuratObj, pcaGenes = NA, umapDims = 1:20){
  library(Seurat)
  
  #Normalizes and selects variable genes
  seuratObj <- NormalizeData(seuratObj)
  seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
  seuratObj <- ScaleData(seuratObj, features = rownames(seuratObj))
  
  #Runs PCA and UMAP
  if (any(is.na(pcaGenes))){
    seuratObj <- RunPCA(seuratObj, features = VariableFeatures(seuratObj))
  } else {
    seuratObj <- RunPCA(seuratObj, features = pcaGenes)
  }
  
  seuratObj <- RunUMAP(seuratObj, dims = umapDims)
  seuratObj
}

removeDuplicates <- function(ensemblsymboldf){
  symbolsVec <- c()
  duplicatedSymbols = unique(ensemblsymboldf$hgnc_symbol[duplicated(ensemblsymboldf$hgnc_symbol)])
  for (x in ensemblsymboldf$hgnc_symbol){
    if (!(x %in% duplicatedSymbols)){
      symbolsVec <- c(symbolsVec, x)
    }
  }
  ensemblsymboldf = data.table(ensemblsymboldf)
  ensemblsymboldf[hgnc_symbol %in% symbolsVec]
}

calculate_overlap <- function(DEGs_subglia, seuratmarkers){
  num_overlap <- list()
  
  for (i in 1:length(unique(seuratmarkers$cluster))){
    cluster_num = i - 1
    num_overlap[[i]] = vapply(unique(DEGs_subglia$microglia_type), 
                              function(x) sum(seuratmarkers[cluster == cluster_num]$gene %in% DEGs_subglia[microglia_type == x]$DEGs),
                              FUN.VALUE = integer(length = 1L))
  }
  df_overlap = do.call(cbind, num_overlap)
  colnames(df_overlap) = unique(seuratmarkers$cluster)
  df_overlap
}

OR<-function(set1,set2,size_universe){
  if(any(duplicated(set1))){
    set1<-unique(set1)
  }
  if(any(duplicated(set2))){
    set2<-unique(set2)
  }
  phyper(q=sum(set1%in%set2)-1, 
         #number of white balls drawn without replacement from an urn which contains both black and white balls. 
         #<=> number of tcf4 target in DEGs. "-1" because normally give P[X > x] but we want P[X >= x])
         m=length(set2), #the number of white balls in the urn. <=> number of DEGs 
         n=size_universe-length(set2), #the number of black balls in the urn. <=> number of genes tested - number of DEGs
         k=length(set1), #the number of balls drawn from the urn <=> numbers of tcf4 targets 
         lower.tail=FALSE)  # if TRUE (default), probabilities are P[X ≤ x] (under-representation), otherwise, P[X > x] (over-representation).
  
}

runOR_acrossastrocyte  <- function(DEGs_state, seuratmarkers, queryseuratobj){
  num_overlap <- list()
  
  #Iterate across cluster
  for (i in 1:length(unique(seuratmarkers$cluster))){
    #Keep track of cluster number
    cluster_num = i - 1
    
    num_overlap[[i]] = vapply(unique(DEGs_state$cell_state), 
                              function(x) OR(DEGs_state[cell_state == x]$DEG, seuratmarkers[cluster == paste("Cluster", cluster_num)]$gene, nrow(queryseuratobj)),
                              FUN.VALUE = double(length = 1L))

  }
  df_overlap = do.call(cbind, num_overlap)
  colnames(df_overlap) = unique(seuratmarkers$cluster)
  df_overlap
}

CalcProp <- function(seuratDEG){
  seuratDEG$cluster = as.character(seuratDEG$cluster)
  countVec <- c()
  
  for (x in unique(seuratDEG$cluster)){
    countVec <- c(countVec, nrow(seuratDEG[cluster == x]))
  }
  return(countVec)
}


computeDirectedKLD <- function(knownPropVector, estimatedPropVector){

  library(LaplacesDemon)
  KLD(knownPropVector, estimatedPropVector, base = 2)$sum.KLD.px.py
}


createResultsDF_ROSMAP <- function(results,trueProportion){

  resultsDF = data.table(results)
  donorList = colnames(trueProportion); cellTypeList = rownames(trueProportion)
  
  donorVec = c(); cellTypeVec = c(); propVec = c()

  for (s in donorList){
    for (c in cellTypeList){
      donorVec = c(donorVec, s)
      cellTypeVec = c(cellTypeVec, c)
      propVec = c(propVec, trueProportion[c,s])
    }
  }
  print(paste(length(donorVec),length(cellTypeVec),length(propVec)))
  resultsSC = data.frame(CellType = cellTypeVec, Sample = donorVec, KnownProportions = propVec, Estimatedproportions = propVec,
                         Method = rep("Zhou Single Cell", length(donorVec)),pearson = rep(0,length(donorVec)), RMSE = rep(0,length(donorVec)),
                         averagepearson = rep(0,length(donorVec)), Dataset = rep("ROSMAP",length(donorVec)))
  
  rbind(resultsDF, data.table(resultsSC))
}


computeDivergenceSC <- function(results, trueProportion){

  resultsDF = createResultsDF_ROSMAP(results, trueProportion)
  resultsSC = resultsDF[Method == "Zhou Single Cell"]
  
  KLDmatrix = matrix(nrow = length(unique(results$Method)), ncol = length(unique(results$CellType)))
  rownames(KLDmatrix) = unique(results$Method)
  colnames(KLDmatrix) = unique(results$CellType)

  #Iterates across methods and cell types and computes KLD
  for (c in colnames(KLDmatrix)){
    for (m in rownames(KLDmatrix)){
      #Estimates PDF from SC object 
      scPDF = density(resultsSC[CellType == c]$Estimatedproportions,from = 0, to = 1, n = 1000)
      
      #Estimates PDF
      methodPDF = density(resultsDF[Method == m][CellType == c]$Estimatedproportions,from = 0, to = 1, n = 1000)
      
      #Computes KLD
      KLDmatrix[m, c] = computeDirectedKLD(scPDF$y, methodPDF$y)
    
    }
  }

  return(KLDmatrix)
}

updateAstCS <- function(wholeSEAAD,astSEAAD){

  astNames = colnames(subset(wholeSEAAD, subset = simple_FCT == "astrocyteofthecerebralcortex"))
  missingAstrocyteNames = astNames[astNames %nin% colnames(astSEAAD)]
  
  wholeSEAAD = wholeSEAAD[,colnames(wholeSEAAD) %nin% missingAstrocyteNames]
  
  newAnnotation = wholeSEAAD@meta.data$simple_FCT
  
  subsetSEAAD = astSEAAD[,colnames(subset(wholeSEAAD, subset = simple_FCT == "astrocyteofthecerebralcortex"))]
  cellStates = subsetSEAAD@meta.data$CellStates
  names(cellStates) = colnames(subsetSEAAD)
  
  
  for (i in 1:length(newAnnotation)){
    if (newAnnotation[i] == "astrocyteofthecerebralcortex"){
      cellName = colnames(wholeSEAAD)[i]
      newAnnotation[i] = cellStates[cellName]
    }
  }
  wholeSEAAD@meta.data$simpleFCTandAstCS = newAnnotation
  return(wholeSEAAD)
}

