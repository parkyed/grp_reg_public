glmnetExtractFeatures <- function(cvfit, x.train){
  
  #' Extract the selected features from a fit glmnet model, using the training set as a reference
  #' 
  #' Function takes a fit cross validated glmnet model, and extract a table of error values and number of features selected for each
  #' lambda value, as well as the indices and names of the features selected by the model at each lambda value
  #'
  #' @param cvfit fit cross validated glmnet object
  #' @param x.train the expression matrix used to train the model. row = samples, columns = features. Just X, no labels in the matrix
  
  # create data frame of performance metrics for each value of lambda
  (lambda.features.df <- data.frame('lambda' = cvfit$lambda,            # all lambda values in solution path
                                    'class.error' = cvfit$cvm,          # cross validation error calculated by glmnet
                                    'num.features' = cvfit$nzero))      # number of selected features (non zero coefficients)
  
  
  # create coefficients object from the model fitted on the full training dataset, at all lambda values
  
  (sparse.coefs.matrix   <- coef(cvfit$glmnet.fit))
  
  # check that the coefs matrix is correct size
  if(!identical((sparse.coefs.matrix@Dim[1]-1), as.numeric(ncol(x.train)))){stop()}
  
  # get the indices of the selected features for each value of lambda, removing the intercept which is indexed at zero
  # use summary to get the ijx format of the matrix to get j, the column number of the value. Just use j, and not i, since i is row number
  # but not the index of the feature, as the sparse matrix indices are 0-based and not 1-based
  # then just remove the intercept after splitting, and the indices will be correct when applied to the feature matrix
  
  ij.matrix <- data.frame(i = sparse.coefs.matrix@i, j = summary(sparse.coefs.matrix)$j)
  
  (fs.idx.all.lambda <- lapply(split(ij.matrix, ij.matrix$j), function(x){x[-1,]}))
  
  # get the names of the features selected (valuable if function called on a training set inside outer cross validation loop
  
  (fs.fids.all.lambda <- sapply(fs.idx.all.lambda, function(x){colnames(x.train)[x$i]}))
  
  return(list('results.df' = lambda.features.df, 'lambda.val'= cvfit$lambda.min, 'fs.idx'= fs.idx.all.lambda, 'sel.feature' = fs.fids.all.lambda))
}