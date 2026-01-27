createTrainTestSplits <- function(counts.mat, y.true, k, r){
  
  # AT SOME POINT BUILD SET SEED INTO THIS SO IT ALWAYS SELECTS THE SAME SEED
  
  #' Generate repeated k fold traning and test sets
  #' 
  #' Combines X matrix and Y labels vector into a single dataframe
  #' Stratified repeated k-fold cross valiation. Stratification based on y.true to ensure same proportion of each class in each fold
  #' Repeated k-fold -to enable averaging of output performance over larger number of train / val splits
  #'
  #' @param counts.mat input matrix, with rows as samples and columns as variables (genes)
  #' @param y.true labels vector
  #' @param k the number of folds int he k-fold partitioning (e.g. 10)
  #' @param r the number of times to repeat the k-fold partitioning
  
  ## set the seed # NEED TO INCLUDE seed AS AN ARGUMENT TO THE FUNCTION
  # if (!missing(seed)) 
  #   set.seed(seed) 
  
  # split the covariate matrix from the class labels
  
  # combine expression matrix with class labels for stratified k-fold sampling
  
  dat.mat <- cbind(counts.mat, y.true) %>% as.data.frame()
  
  # create an k-fold rsplit object that defines train and test sets 
  
  set.seed(456)
  
  data.splits <- rsample::vfold_cv(dat.mat, v = k, repeats = r, strata = 'y.true')
  
  # extract the train and test sets as a dataframe
  
  train.test.data <-  data.splits %>% 
    dplyr::mutate(train = map(splits, ~rsample::analysis(.x)),
                  test = map(splits, ~rsample::assessment(.x))) %>% 
    dplyr::select(train, test)
  
  return(train.test.data)
}