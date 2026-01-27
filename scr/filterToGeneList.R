filterToGeneList <- function(counts.mat, samp.list, g.list){
  
  #' helper function to filter an expression matrix based on a list of samples and list of genes
  #'
  #' @param counts.mat expression matrix with patients as rows and genes as column
  #' @param samplist list of samples
  #' @param g.list list of genes
  
  counts.mat.gene <- counts.mat[which(rownames(counts.mat) %in% samp.list), which(colnames(counts.mat) %in% g.list)]
  
  return(counts.mat.gene)
}