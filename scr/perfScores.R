# ======== Helper Functions ===========

# manual calculation of the predictions based on the linear predictors (log odds from the model)

inv.logit <- function(x) exp(x)/(1+exp(x))

# calculating class predictions from log loss  / linear predictors for a given cut off threshold

LogLoss2ClassPred <- function(y.pred.log.loss, t=0.5){
  class.pred <- (inv.logit(y.pred.log.loss) > t) * 1
  class.pred <- factor(class.pred, levels = c(0,1))
  return(class.pred)
}

# calculation of classification error, sensitivity, specificity, etc. using caret

SensSpecScores <- function(cv.pred.df, t,  pclass = "1"){
  y.pred.class <- LogLoss2ClassPred(cv.pred.df$y.pred, t=t)
  caret.conf.mat <- caret::confusionMatrix(y.pred.class, factor(cv.pred.df$y.true, levels=c(0,1)), positive = pclass)
  (accuracy <- caret.conf.mat$overall['Accuracy'])
  (sensitivity <- caret.conf.mat$byClass['Sensitivity'])
  (specificity <- caret.conf.mat$byClass['Specificity'])
  (precision <- caret.conf.mat$byClass['Precision'])
  (f1 <- caret.conf.mat$byClass['F1'])
  return(list('accuracy' = accuracy, 'sensitivity' = sensitivity, 'specificity' = specificity, 'precision' = precision, 'f1' = f1))
}

# ======== Scoring Functions ===========

perfScoresAllLambda <- function(y.pred, y.true){
  
  #' Calculate performance metrics given y.pred and y.true
  #' 
  #' Returns single column data frame of scores (auc, accuracy, sensitivity, specificity, precision, f1) for given y.pred and true
  #'
  #' @param y.pred vector of predictions given as log odds
  #' @param y.true vector of true class labels
  
  (cv.pred.df <- data.frame(y.pred, y.true))
  colnames(cv.pred.df) <- c("y.pred", "y.true")
  (perf.scores.test <- assess.glmnet(cv.pred.df$y.pred, newy = cv.pred.df$y.true, family = "binomial"))
  (sens.spec.scores <- SensSpecScores(cv.pred.df, 0.5,  pclass = "1"))
  (auprc <- MLmetrics::PRAUC(cv.pred.df$y.pred, cv.pred.df$y.true)) # added auprc
  (perf.score.df <- t(round(data.frame(c('auc'= perf.scores.test$auc, 'auprc' = auprc, lapply(sens.spec.scores, unname))),3))) # added auprc
  return(perf.score.df)
}


perfScoresDFAllFolds <- function(cv.preds.all.folds, calc.means = F){
  
  #' Calculate average performance over all folds for at each lambda value
  #' 
  #' Combines the predictions made by the repeated k-fold cross validation across all folds
  #' Then scores the predictions together at each value of lambda. Note: assumption that lambda values in the solution path
  #' are the same for each folds, which is only true is a fixed set of lambda values provided as input to model training
  #' This would be solved by training the model using a fixed set of lambda values, reused over every fold
  #'
  #' @param cv.preds.all.folds results of nested cross validation, including y.pred for each lambda, for every fold.
  
  # extract the predictions from each fold
  
  y.preds.list <- sapply(cv.preds.all.folds, function(x){x$y.pred.all.lambda})
  
  # convert resulting matrices to dataframes
  
  y.preds.list <- lapply(y.preds.list, as.data.frame)
  
  # row bind all the data frames - use dplyr::bind_rows so that columns that don't match are set as NA
  
  y.preds.all.lambda <- Reduce(dplyr::bind_rows, y.preds.list)
  
  # get indices of columns not containing NA values
  
  not.na.idx <- which(colSums(is.na(y.preds.all.lambda))==0) # added 12.11.24
  
  # remove all columns that contain NA values for the prediction

  y.preds.all.lambda <- y.preds.all.lambda %>% dplyr::select(all_of(not.na.idx)) # added 12.11.24
  
  # get the true class labels for the examples in every fold, and combine
  (y.true <- sapply(cv.preds.all.folds, function(x){x$y.true}) %>% unlist())
  
  # calculate the mean number of features selected at each value of lambda # EXCLUDE FOR NOW, AS RELIES ON SOLUTION FOR ALL LAMBDA
  if(calc.means){
    mean.features <- sapply(cv.preds.all.folds, function(x){lapply(x$fs.results$sel.feature, length)}) %>%
       as.numeric %>%
       matrix(ncol = length(cv.preds.all.folds)) %>%
       rowMeans
  } else {
    mean.features = rep(NA, length(y.preds.all.lambda))
  }
  
  # score the log odds predictions at each value of lambda
  (scores.all.lambda <- y.preds.all.lambda %>% as.data.frame() %>% map(perfScoresAllLambda, y.true))
  (auc.all.lambda <- lapply(scores.all.lambda, function(x){x["auc",]}) %>% unlist %>% unname)
  (auprc.all.lambda <- lapply(scores.all.lambda, function(x){x["auprc",]}) %>% unlist %>% unname) # added auprc
  (accuracy.all.lambda <- lapply(scores.all.lambda, function(x){x["accuracy",]}) %>% unlist %>% unname)
  (sensitivity.all.lambda <- lapply(scores.all.lambda, function(x){x["sensitivity",]}) %>% unlist %>% unname)
  (specificity.all.lambda <- lapply(scores.all.lambda, function(x){x["specificity",]}) %>% unlist %>% unname)
  (f1.all.lambda <- lapply(scores.all.lambda, function(x){x["f1",]}) %>% unlist %>% unname)
  
  # combine scores in output data frame
  (scores.df.all.lambda <- data.frame(lambda = cv.preds.all.folds[[1]]$fs.results$results.df$lambda[not.na.idx], # added 12.11.24 to allow for incomplete solution path
                                      mean.features = mean.features,
                                      auc = auc.all.lambda,
                                      auprc = auprc.all.lambda, # added auprc
                                      accuracy = accuracy.all.lambda,
                                      sensitivity = sensitivity.all.lambda,
                                      specificity = specificity.all.lambda,
                                      f1 = f1.all.lambda)
  ) 
  
  return(scores.df.all.lambda)
}

