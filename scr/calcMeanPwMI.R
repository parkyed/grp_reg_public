calcMeanPwMI <- function(g.list, counts.mat){
 
  #' Filter a gene expression array on a set of genes, and calculate the mean pair-wise mutual information between the genes
  #'
  #' @param counts.mat expression matrix with patients as rows and genes as column
  #' @param g.list list of genes to include in calculation
  
  x.mat.fil.gene <- filterToGeneList(counts.mat = counts.mat,
                                     samp.list = rownames(counts.mat),
                                     g.list = g.list)
  
  # use mpmi to calculate the matrix of pairwise mutual information
  mi.res <- mpmi::cmi(x.mat.fil.gene)
  
  # extract the matrix
  mi.mat <- mi.res$mi
  
  # look at the upper triangle only, including the diagonals to remove duplicate values, and take arithmetic mean
  M <- mi.mat[upper.tri(mi.mat, diag = F)] %>% mean
  
  return(M)
}