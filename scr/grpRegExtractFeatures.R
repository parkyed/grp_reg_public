# =================== function to extract selected features and groups ======================

# main grp reg extract features function - building in the group indices selection


GrpRegExtractFeaures <- function(cv.grp.reg.fit, fmax=NA){
  
  #' Feature extraction function for group regression objects
  #' 
  #' Extracts the selected features and groups at each value of lambda in the solution path for a fitted grpreg object.
  #' Returns a dataframe of error, number of features and number of groups selected by lambda, the selected features and selected groups for all lambda, and the best value of lambda, based either on minimum error, or on the input min-max feature range if provided
  #'
  #' @param cv.grp.reg.fit a fit grpreg model
  #' @param fmax the maximum number of features in the target feature range
  
  
  # extract the genes with non-zero coefficients from the beta matrix for each value of lambda for the model fit on the training data
  # intercept removed first
  
  (gene.by.lambda <- apply(cv.grp.reg.fit$fit$beta[-1,], 2, function(x){subset(x, subset= x !=0)}))
  
  # get the gene names from the named vectors of coefficients
  
  (sel.genes <- lapply(gene.by.lambda, names))
  
  # create vector of cross validation errors from model that matches the length of all lambda values
  # the $fit$lambda is the input lambda values - all values in the solution path, the $lambda are those acutally calculated - the NAs fill the gap
  
  (cve.vec <- c(cv.grp.reg.fit$cve, rep(NA, (length(cv.grp.reg.fit$fit$lambda)-length(cv.grp.reg.fit$lambda)))))
  
  # get the beta latents table from the fitted model, exclude the intercept
  
  (beta.lat.coefs <- cv.grp.reg.fit$fit$beta.latent[-1,])
  
  # convert to named list - each items is a lambda, with all coefs at that lambda
  
  beta.lat.coefs.lst <- lapply(seq_len(ncol(beta.lat.coefs)), function(i) beta.lat.coefs[,i])
  
  names(beta.lat.coefs.lst) <- colnames(beta.lat.coefs)
  
  # map to group index identification helper function
  
  (sel.grps <- lapply(beta.lat.coefs.lst, latentBetaSelGroups))
  
  # create dataframe of performance metrics for each value of lambda - to identify best lambda, based on either error or # features
  
  (lambda.features.df <- data.frame('lambda' = cv.grp.reg.fit$fit$lambda,         # all lambda values in solution path
                                    'cross.val.error' = cve.vec,                  # cross validation error calculated by grpreg
                                    'num.features' = lengths(gene.by.lambda),     # calculated number of features at each lambda based on the non-zero betas
                                    'n.sel.grps' = lengths(sel.grps)))            # number of selected groups based on the predict function
  
  if(!is.na(fmax)){
    
    # filter the data frame to only the lambda values that result in a number of features in the fmin fmax range
    
    (lambda.features.df <- dplyr::filter(lambda.features.df, num.features <= fmax))
    
    # best lambda is the minimum of values that select number of features in the desired range => model with number of features closest to the fmax
    
    (lambda.best  <- min(lambda.features.df$lambda))
    
  } else {
    
    # the best lambda is the one that minimises the error
    
    lambda.best <- cv.grp.reg.fit$lambda.min
  }
  
  # return the error and the features selected by the model
  
  return(list('results.df' = lambda.features.df,
              'lambda.best'= lambda.best, 
              'sel.feature' = sel.genes, 
              'sel.grps' = sel.grps)
  )
}

# helper function to extract indices of the groups from the named list of latent beta coefficients

latentBetaSelGroups <- function(coef.lst){
  
  (lst.filtered <- coef.lst[coef.lst !=0])
  
  (grp.names <- names(lst.filtered))
  
  (grp.names.vec <- lapply(grp.names, stringr::str_split_i, pattern = '_', i = 1) %>% unlist)
  
  (grp.idx <- as.numeric(gsub("grp", "", grp.names.vec)) %>% unique)
  
  return(grp.idx)
}



# function to extract selected genes and pathways for given number of features range
# calls the GrpRegExtractFeaures above for a specific number of features - min / max

maxFeatExtract <- function(grpreg.fit.model, max.feat){
  
  #' identify selected features, by group for given number of features
  #'
  #' @param grpreg.fit.model a fit grpreg model
  #' @param max.feat the maximum number of features in the target feature range
  
  feat.extract <- GrpRegExtractFeaures(grpreg.fit.model, fmax=max.feat)
  
  (feat.names <- feat.extract$sel.feature[[as.character(round(feat.extract$lambda.best,4))]])
  
  (feat.groups <- grpreg.fit.model$fit$group[feat.extract$sel.grps[[as.character(round(feat.extract$lambda.best,4))]]])
  
  (feat.names.groups <- lapply(feat.groups, FUN = function(x){x[x %in% feat.names]}))
  
  return(list('n.features' = length(feat.names), 'feat.names' = feat.names, 'feat.names.groups' = feat.names.groups))
}
