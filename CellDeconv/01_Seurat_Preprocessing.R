library("Seurat")
library("Biobase")
library("rjson")
library(pbapply)


seuratNormalization <- function(inputMatrix){
  #Performs seurat normalization scheme for bulk matricies (Column is identity value and row is gene)
  
  #Begins with column normalization (dividing by total count of a sample)
  scaleFactor = 10000 #Scale factor used by default in seurat
  inputMatrix = pbapply(inputMatrix, 2, function(x) log1p((x/sum(x)) * scaleFactor))
  
  return(inputMatrix)
  
}

convertFormat <- function(readDF){
  readDF = data.frame(readDF)
  rownames(readDF) = readDF[[1]]
  readDF = readDF[,2:ncol(readDF)]
  readDF
}



preprocessData <- function(){

  inputSeurat <<- readRDS(refFile)
  
  DefaultAssay(inputSeurat) <- "RNA"
  
  #Normalize
  inputSeurat <<- NormalizeData(inputSeurat, normalization.method = "LogNormalize",assay = "RNA")
  
  #Edit cell type names 
  celltype <<- as.character(inputSeurat@meta.data[[cellTypeSeurat]])
  inputSeurat@meta.data[[cellTypeSeurat]] <<- as.factor(celltype)

  #Creates referenceMatrix
  if (referenceNormalization == "Seurat"){
    referenceProfile <<- AverageExpression(inputSeurat,group.by = cellTypeSeurat)[["RNA"]]
  } else if (referenceNormalization == "None"){
    referenceProfile <<- AverageExpression(inputSeurat, group.by = cellTypeSeurat, slot = "counts")[["RNA"]]
    #referenceProfile <<- cbind(referenceProfile[,1:(((1:ncol(referenceProfile))[colnames(referenceProfile) == "Excitatory_Neuron"]) -1), drop = FALSE],
                              #referenceProfile[,(((1:ncol(referenceProfile))[colnames(referenceProfile) == "Excitatory_Neuron"]) + 1):ncol(referenceProfile), drop = FALSE])
  }

  
  #Loads in bulk matrix if specified. Otherwise, creates matrix from "querySeurat" object 
  if (pseudobulkmixpath != ""){
    
    #Loads in matrix
    pseudobulkMix <<- fread(pseudobulkmixpath)
    pseudobulkMix <<- as.matrix(convertFormat(pseudobulkMix))
    #pseudobulkMix <<- as.matrix(apply(pseudobulkMix, 2,  function(x) log2(((x/sum(x)) + 1) * 10^6)))

    #Performs normalization
    #if (bulkNormalization == "Seurat"){
      #pseudobulkMix <<- seuratNormalization(pseudobulkMix)
    #}
    
    #Loads in true proportion List
    if (trueProportionPath != ""){
     
      trueProportionList = convertFormat(read.csv(trueProportionPath))
    }else{
      trueProportionList=list()
    }
    
  } else {
 
    
    querySeurat <<- readRDS(queryFile)
    querySeurat <<- NormalizeData(querySeurat, normalization.method = "LogNormalize",assay = "RNA")
    
    querySeurat@meta.data[[sample_metadata]] <<- as.factor(as.character(querySeurat@meta.data[[sample_metadata]]))
    querySeurat@meta.data[[cellTypeSeurat]] <<- as.factor(as.character(querySeurat@meta.data[[cellTypeSeurat]]))
    
    #Splits seurat object to determine proportions
    splitSeurat = SplitObject(querySeurat,split.by = sample_metadata)
    
    
    trueProportion = as.data.frame.matrix(table(querySeurat@meta.data[[cellTypeSeurat]], querySeurat@meta.data[[sample_metadata]]))
    trueProportion = data.frame(apply(trueProportion,2,function(x) x/sum(x)))
    trueProportionList = trueProportion
    
    
    #Creates  pseudobulkMix
    if (bulkNormalization == "Seurat"){
      pseudobulkMix <<- AggregateExpression(querySeurat,group.by = sample_metadata)[["RNA"]]
    } else if (bulkNormalization == "None") {
      pseudobulkMix <<- AggregateExpression(querySeurat, group.by = sample_metadata, normalization.method = NULL, slot = "counts")[["RNA"]]
    }
    
  }
  sharedGenes = rownames(pseudobulkMix)[rownames(pseudobulkMix) %in% rownames(referenceProfile)]
  pseudobulkMix <<- data.frame(pseudobulkMix)
  referenceProfile <<- data.frame(referenceProfile)
  pseudobulkMix <<- pseudobulkMix[sharedGenes,]
  referenceProfile <<- referenceProfile[sharedGenes,]
  
  pseudobulkMix <<- as.matrix(pseudobulkMix)
  referenceProfile <<- as.matrix(referenceProfile)
  
  return(trueProportionList)

}
