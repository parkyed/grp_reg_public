## function: cross-validation
# ------------------------------------------------------------------------------
cv.grpregOverlap_EP <- function(X, y, group, ...) {
  fit <- grpregOverlap_EP(X=X, y=y, group=group, returnX = TRUE, ...)
  cvfit <- cv.grpreg(X = fit$X.latent, y = y, group = fit$grp.vec, ...)
  cvfit$fit <- fit
  val <- structure(cvfit, class = c('cv.grpregOverlap', 'cv.grpreg'))
  val
}
# ------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
## function: overlapping group selection based on R Package 'grpreg' 

grpregOverlap_EP <- function(X, y, group, 
                             family=c("gaussian","binomial", "poisson", 'cox'), 
                             returnX.latent = FALSE,
                             returnOverlap = FALSE,
                             ...) {
  
  # Error checking
  if (is.matrix(X)) {
    tmp <- try(X <- as.matrix(X), silent=TRUE)
    if (class(tmp)[1] == "try-error")  {
      stop("X must be a matrix or able to be coerced to a matrix")
    }   
  }
  if (storage.mode(X)=="integer") X <- 1.0*X
  
  incid.mat <- incidenceMatrix(X, group) # group membership incidence matrix
  over.mat <- over.temp <- Matrix(incid.mat %*% t(incid.mat)) # overlap matrix
  grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat)) # group index vector
  X.latent <- expandX_EP(X, group)
  
  diag(over.temp) <- 0
  if (all(over.temp == 0)) {
    cat("Note: There are NO overlaps between groups at all!", "\n") 
    cat("      Now conducting non-overlapping group selection ...")
  }
  
  family <- match.arg(family)
  if (family != 'cox') {
    fit <- grpreg(X = X.latent, y = y, group = grp.vec, family = family, ...) # fit the grpreg model here
  } else {
    ## survival analysis
    fit <- grpsurv(X = X.latent, y = y, group = grp.vec, ...)
  }
  
  fit$beta.latent <- fit$beta # fit$beta from grpreg is latent beta
  fit$beta <- gamma2beta(gamma = fit$beta, incid.mat, grp.vec, family = family)
  fit$incidence.mat <- incid.mat
  fit$group <- group
  fit$grp.vec <- grp.vec # this is 'group' argument in Package 'grpreg'
  fit$family <- family
  if (returnX.latent) {
    fit$X.latent <- X.latent
  } 
  if (returnOverlap) {
    fit$overlap.mat <- over.mat
  }
  
  if (family != 'cox') {
    # get results, store in new class 'grpregOverlap', and inherited from 'grpreg'
    val <- structure(fit,
                     class = c('grpregOverlap', 'grpreg'))
  } else {
    val <- structure(fit, 
                     class = c("grpsurvOverlap", "grpregOverlap"))
  }
  val
}
# -------------------------------------------------------------------------------

## function: convert latent beta coefficients (gamma's) to non-latent beta's

# -------------------------------------------------------------------------------
gamma2beta<- function(gamma, incidence.mat, grp.vec, family) {
  # gamma: matrix, ncol = length(lambda), nrow = # of latent vars.
  p <- ncol(incidence.mat) # number of genes in dataset
  J <- nrow(incidence.mat) # number of gene groups (including orphan groups)
  beta <- matrix(0, ncol = ncol(gamma), nrow = p)
  
  if (family != 'cox') {
    intercept <- gamma[1, , drop = FALSE] # extract the intercept as that isn't partitioned into groups
    gamma <- gamma[-1, , drop = FALSE]
  } else {
    # Cox model doesn't have an intercept
    gamma <- gamma
  }
  
  for (i in 1:J) {
    ind <- which(incidence.mat[i, ] == 1)
    beta[ind, ] <- beta[ind, ] + gamma[which(grp.vec == i), , drop = FALSE]  # understand this better, basically adds up the gammas to be betas
  }
  if (family != 'cox') {
    beta <- rbind(intercept, beta)                 # add back the intercept once beta has been reconstructed from gamma
    rownames(beta) <- c("(Intercept)", colnames(incidence.mat))
  } else {
    rownames(beta) <- colnames(incidence.mat)
  }
  beta
}
# -------------------------------------------------------------------------------


## function: expand design matrix X to overlapping design matrix (X.latent)

# -------------------------------------------------------------------------------
expandX_EP <- function(X, group) {
  incidence.mat <- incidenceMatrix(X, group) # group membership incidence matrix
  over.mat <- Matrix(incidence.mat %*% t(incidence.mat), sparse = TRUE, 
                     dimnames=list(rownames(incidence.mat),
                                   rownames(incidence.mat))) # overlap matrix
  # dimnames = dimnames(incidence.mat)) # THIS WAS THE LINE CAUSING AN ERROR - REPLACED WITH THE ABOVE
  grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat)) # group index vector
  
  # expand X to X.latent
  X.latent <- NULL
  names <- NULL
  
  ## the following code will automatically remove variables not included in 'group'
  ## solution is to ensure that the group input includes orphan genes as their own group - that way, included in X.latent
  
  for(i in 1:nrow(incidence.mat)) {
    idx <- incidence.mat[i,]==1  # index
    X.latent <- cbind(X.latent, X[, idx, drop=FALSE])
    names <- c(names, colnames(incidence.mat)[idx])
    #     colnames(X.latent) <- c(colnames(X.latent), colnames(X)[incidence.mat[i,]==1])
  }
  colnames(X.latent) <- paste('grp', grp.vec, '_', names, sep = "")
  X.latent
}
# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------

## function: incidence matrix: I[i, j] = 1 if group i contains variable j.

incidenceMatrix <- function(X, group) {
  n <- nrow(X)
  p <- ncol(X)
  if (! is.list(group)) {
    stop("Argument 'group' must be a list of integer indices or character names of variables!")
  }
  J <- length(group)
  
  # create grp.mat of zeros
  grp.mat <- Matrix(0, nrow = J, ncol = p, sparse = TRUE, 
                    dimnames=list(as.character(rep(NA, J)),
                                  as.character(rep(NA, p))))  
  
  # create temporary names for columms and rows of grp.mat
  if(is.null(colnames(X))) {
    colnames(X) <- paste("V", 1:ncol(X), sep="")    
  }
  if (is.null(names(group))) {
    names(group) <- paste("grp", 1:J, sep="")
  }
  
  if (is.numeric(group[[1]])) {
    for (i in 1:J) {
      ind <- group[[i]]
      grp.mat[i, ind] <- 1
      colnames(grp.mat)[ind] <- colnames(X)[ind]
    }
  } else { ## character, names of variables
    for (i in 1:J) {
      grp.i <- as.character(group[[i]]) # set the names of each gene in the group as a character to the grp.i variable
      ind <- colnames(X) %in% grp.i     # create logical of the with membership of group i over all genes
      grp.mat[i, ] <- 1*ind             # create a row for this group in the grp.mat matrix, 1 => gene in group
      colnames(grp.mat)[ind] <- colnames(X)[ind] # add in the column names to the matrix (not sure why this is done inside the loop, order doesn't change)
    }
  }
  
  # add group names to be the rownames of the matrix
  rownames(grp.mat) <- as.character(names(group))
  
  # check grp.mat has at least some genes present across the groups - this just checks if nothing has happened
  if (all(grp.mat == 0)) {
    stop("The names of variables in X don't match with names in group!")
  }
  
  return(grp.mat)
}
# -------------------------------------------------------------------------------


## function: predict cv.grpregOverlap_EP

## cross validation function that calls predict.grpregOverlap_EP instead of the original function
# -------------------------------------------------------------------------------
predict.cv.grpregOverlap_EP <- function(object, X, type=c("link", "response", "class", "coefficients", "vars", "groups", "nvars", "ngroups", "norm"),
                                        latent = FALSE, lambda = object$lambda.min,
                                        which=object$min, ...) {
  type <- match.arg(type)
  predict.grpregOverlap_EP(object$fit, X=X, type=type, latent=latent, lambda=lambda,
                           which=which)
}


coef.cv.grpregOverlap <- function(object, latent = FALSE, lambda = object$lambda.min,
                                  which = object$min, ...) {
  coef(object$fit, lambda=lambda, latent=latent, which=which, ...)
}
# -------------------------------------------------------------------------------


## function: predict.grpregOverlap_EP

# function adjusted to use ExpandX_EP resolving dimensions issue
# -------------------------------------------------------------------------------
predict.grpregOverlap_EP <- function(object, X, 
                                     type=c("link", "response", "class", 
                                            "coefficients", "vars", "groups", 
                                            "nvars", "ngroups", "norm"), 
                                     latent = FALSE, lambda, 
                                     which=1:length(object$lambda), ...) {
  family <- object$family
  if (!missing(X) && is.character(X)) {
    type <- X
    X <- NULL
  }
  type <- match.arg(type)
  beta <- coef(object=object, lambda=lambda, latent = latent, 
               which=which, drop = FALSE, ...)
  if (type == 'coefficients') return(beta)
  if (length(dim(object$beta)) == 2) {
    if (type=="vars") {
      if (family != 'cox') {
        return(drop(apply(beta[-1, , drop=FALSE]!=0, 2, FUN=which)))
      } else {
        return(drop(apply(beta != 0, 2, FUN=which)))
      }
    }
    if (type=="nvars") {
      if (family != 'cox') {
        v <- drop(apply(beta[-1, , drop=FALSE]!=0, 2, FUN=which))
      } else {
        v <- drop(apply(beta != 0, 2, FUN=which))
      }
      if (is.list(v)) {
        res <- sapply(v, length)
      } else {
        res <- length(v)
      }
      return(res)
    }
    if (type == 'norm') {
      beta <- coef(object=object, lambda=lambda, latent = TRUE, 
                   which=which, drop = FALSE, ...)
      if (!latent) {
        cat("The returned is the L2 norm of the latent coefficients! Set latent = 'TRUE' to avoid this warning message.\n\n")
      } 
      if (family != 'cox') {
        return(drop(apply(beta[-1, , drop=FALSE], 2, 
                          function(x) tapply(x, object$grp.vec, function(x){sqrt(sum(x^2))}))))
      } else {
        return(drop(apply(beta, 2, function(x) tapply(x, object$grp.vec, function(x){sqrt(sum(x^2))}))))
      }
    }
  } else {
    if (type=="vars") 
      stop("Predicting type 'vars' not implemented with multivariate outcomes\n")
    if (type=="nvars") {
      return(drop(apply(beta[,-1, , drop=FALSE]!=0, 3, FUN=sum)))
    }
    if (type == 'norm') {
      beta <- coef(object=object, lambda=lambda, latent = TRUE, 
                   which=which, drop = FALSE, ...)
      if (!latent) {
        cat("The returned is the L2 norm of the latent coefficients. Set latent = 'TRUE' to avoid this warning message.\n\n")
      } 
      if (family != 'cox') {
        return(drop(apply(beta[, -1, , drop=FALSE], 3, function(x) apply(x, 2, function(x){sqrt(sum(x^2))})))) 
      } else {
        return(drop(apply(beta, 3, function(x) apply(x, 2, function(x){sqrt(sum(x^2))})))) 
      }
    }
  }
  if (!missing(X) && !is.null(X)) {
    X <- expandX_EP(X, object$group) # THIS IS THE ONLY LINE THAT HAS CHANGED FROM THE ORIGINAL
  }
  if (latent) {
    cat("Only latent 'coefficients', 'vars', 'nvars', 'norm' can be returned! Set latent = 'FALSE' to suppress this message.\n\n")
  }
  obj.new <- object
  obj.new$group <- object$grp.vec
  obj.new$beta <- object$beta.latent
  if (family != 'cox') {
    class(obj.new) <- 'grpreg'
  } else {
    class(obj.new) <- 'grpsurv'
  }
  return(predict(obj.new, X=X, type=type, lambda=lambda, which=which, ...))
}
# -------------------------------------------------------------------------------

## function: coef.grpregOverlap, coef for grpregOverlap
# -------------------------------------------------------------------------------

coef.grpregOverlap <- function(object, lambda, latent = FALSE, 
                               which=1:length(object$lambda), drop=TRUE, ...) {
  family <- object$family
  obj.new <- object
  obj.new$beta <- object$beta.latent
  class(obj.new) <- 'grpreg'
  
  ## latent beta
  beta <- coef(object = obj.new, lambda = lambda, 
               which = which, drop=FALSE, ...)  
  if (!latent) {
    beta <- gamma2beta(beta, incidence.mat = object$incidence.mat, 
                       grp.vec = object$grp.vec, family = family)
  }
  if (drop) return(drop(beta)) else return(beta)
}
# -------------------------------------------------------------------------------



