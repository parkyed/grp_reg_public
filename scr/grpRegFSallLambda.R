grpRegFSallLambda <- function(train, test, inner.k, pen, grp.lst, alpha=0.4, lambda.min=0.05, lambda.path = NULL, tau=0.3, fmax=NA){
  
  #' Train, predict and feature extraction function
  #' 
  #' Trains group regression model on training set,
  #' Extracts selected features at each lambda value in the trained model
  #' Makes predictions on the test set, at all values of lambda.
  #' Returns a dataframe of error and features selected at leach lambda, best lamda (min error or feature range)
  #' Set up to be mapped to k-fold split dataset
  #'
  #' @param train the input matrix of counts - training set - either log2, vst or rlog transformed. row = samples, columns = features. last column is labels
  #' @param test the matrix of counts to predict on - test set - either log2, vst or rlog transformed. row = samples, columns = features. last column is labels
  #' @param inner.k the number of folds in the inner cross validation loop for the group regression function, can be set to loocv
  #' @param pen the type of group penalty to fit the model with (see grpreg options)
  #' @param grp.lst named list of the gene groups, named by pathway or group name
  #' @param alpha determine the balance between the group penalty and a ridge penalty
  #' @param lambda.min the minimum value for lambda in the solution path
  #' @param tau tuning parameter that applies to the gel penalty only
  #' @param fmax the maximum number of features required - used in calculating lambda.best for given number of features
  
  # split the covariate matrix from the class labels
  
  y.train  <- train$y.true
  y.test   <- test$y.true
  x.train  <- train[ , !names(train) %in% c("y.true")]
  x.test   <- test[ , !names(test) %in% c("y.true")]
  
  # determine number of folds for the inner cross validation loop to calculate the classification error at each lambda value
  
  if(inner.k=="loocv"){inner.k <- length(y.train)}
  
  # fit cross validated group regression model model (updated _EP version for overlapping groups)
  
  cv.grp.reg.fit <- cv.grpregOverlap_EP(X = as.matrix(x.train),
                                        y = y.train,
                                        group = grp.lst,
                                        penalty= pen,
                                        family= "binomial",
                                        nlambda = 100,
                                        lambda = lambda.path,
                                        lambda.min = lambda.min,        
                                        k = inner.k,
                                        alpha = alpha,
                                        tau = tau)
  
  # extract the dataframe of features selected at each value of lambda
  # extract the best value of lambda - at either min error, or target max features where provided
  
  fs.res <- GrpRegExtractFeaures(cv.grp.reg.fit, fmax=fmax)
  
  # predict the log odds values (linear predictors) on the test set using the model, at all lambda values, and at best lambda
  
  (y.pred.lbest.lodds <- predict.cv.grpregOverlap_EP(cv.grp.reg.fit, X = as.matrix(x.test), type="link", lambda = fs.res$lambda.best))
  
  (y.pred.all.lambda.lodds <- predict.cv.grpregOverlap_EP(cv.grp.reg.fit, X = as.matrix(x.test), type="link", lambda = cv.grp.reg.fit$lambda))
  
  # return the features dataframe, predictions as log odds, and y.test for scoring
  
  return(list(#'cv.fit.train' = cv.grp.reg.fit, 
    'fs.results' = fs.res,
    'y.pred.l.best' = y.pred.lbest.lodds,
    'y.pred.all.lambda' = y.pred.all.lambda.lodds,
    'y.true' = y.test))
}