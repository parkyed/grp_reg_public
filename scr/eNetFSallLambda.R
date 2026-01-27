eNetFSallLambda <- function(train, test, inner.k = 10, alpha=0.95, lambda.path = NULL, penalties = NULL){
  
  #' Nested cross-validation feature selection using elastic net, returning predictions on the test set, at all values of lambda
  #'
  #' @param train expression matrix for training, including labels - log2 or vst transformed. row = samples, columns = features. last column is labels
  #' @param test expression matrix test set, including labels - log2 or vst transformed. row = samples, columns = features. last column is labels
  #' @param inner.k the number of folds in the inner cross validation loop for the regression function, can be set to loocv
  #' @param alpha determine the balance between the lasso penalty and a ridge penalty
  #' @param lambda.path a pre-defined path of lambda values to be tested (this ensures that same lamda values are used for each fold)
  
  # split the covariate matrix from the class labels
  
  y.train  <- train$y.true
  y.test   <- test$y.true
  x.train  <- train[ , !names(train) %in% c("y.true")]
  x.test   <- test[ , !names(test) %in% c("y.true")]
  
  # determine number of folds for the inner cross validation
  
  if(inner.k=="loocv"){inner.k <- length(y.train)}
  
  # set the default penalty factor if none given
  
  if(is.null(penalties)){penalties = rep(1, ncol(x.train))}
  
  # fit cross validated elastic net model
  
  cvfit <- cv.glmnet(x= as.matrix(x.train),
                     y = y.train,
                     family='binomial',
                     alpha=alpha,
                     lambda = lambda.path,
                     penalty.factor = penalties,
                     type.measure = 'class',
                     nfolds = inner.k, grouped=FALSE, keep = TRUE, parallel = FALSE)
  
  # extract and identify the selected features for each value of lambda
  
  fs.res <- glmnetExtractFeatures(cvfit, x.train)
  
  
  # predict the log odds values (linear predictors) on the test set using the model, at all lambda values, and at best lambda
  
  (y.pred.lbest.lodds <- predict(cvfit, newx = as.matrix(x.test), s = cvfit$lambda.min))
  
  (y.pred.all.lambda.lodds <- predict(cvfit, newx = as.matrix(x.test), s = cvfit$lambda))
  
  # return the predictions as log odds, along with the best lambda value and features selected on the train set
  
  return(list('fs.results' = fs.res,
              'y.pred.l.best' = y.pred.lbest.lodds,
              'y.pred.all.lambda' = y.pred.all.lambda.lodds,
              'y.true' = y.test,
              'best.lambda' = fs.res$lambda.val,
              'results.df' = fs.res$results.df,
              'fs.idx' = fs.res$fs.idx,
              'fs.fids' = fs.res$sel.feature))
  
}