# functions that takes the sample information, gene list and counts matrix and runs the DE analysis
# inputs are limma inputs - exp.mat is a matrix like object containing log expression values for a series of arrays
# with rows corresponding to genes and columns corresponding to samples
# outputs table of fc and p.vals to be used in a volcano plot


calcDECondition <- function(exp.mat, samp.info, gene.lst, res.number = 100){

  exp.mat.f <- exp.mat[gene.lst, samp.info$sampleID]

  # create simple design matrix, that allows for the paired samples by patients study id and the time factor

  condition <- factor(samp.info$condition)

  (design <- stats::model.matrix(~condition))

  # fit expression matrix to a linear model

  fit <- limma::lmFit(exp.mat.f, design)

  fit <- limma::eBayes(fit)

  # Generate a table with the output metrics for all the genes. Use BH correction for multiple testing

  de.res <- limma::topTable(fit, number = res.number, adjust = "BH")

  return(de.res)
}

calcDEDay <- function(exp.mat, samp.info, gene.lst){

  # exp.mat.f <- exp.mat[samp.info$sampleID, gene.lst]
  exp.mat.f <- exp.mat[gene.lst, samp.info$sampleID]

  patient.id <- factor(samp.info$patient.study.id)

  day <- factor(samp.info$day)

  (design <- stats::model.matrix(~patient.id + day))

  # fit <- limma::lmFit(t(exp.mat.f), design)
  fit <- limma::lmFit(exp.mat.f, design)

  fit <- limma::eBayes(fit)

  de.res <- limma::topTable(fit, coef = "day3", number = Inf, adjust = "BH")

  return(de.res)
}
