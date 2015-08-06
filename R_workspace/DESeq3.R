library(DESeq2)

# Modified functions from package DESeq2 for improved performance

rlogData <- function(object, intercept, betaPriorVar) {
  if (is.null(mcols(object)$dispFit)) {
    stop("first estimate dispersion with a design of formula(~ 1)")
  }
  samplesVector <- as.character(seq_len(ncol(object)))
  if (!missing(intercept)) {
    if (length(intercept) != nrow(object)) {
      stop("intercept should be as long as the number of rows of object")
    }
  }
  if (!"allZero" %in% names(mcols(object))) {
    mcols(object)$allZero <- rowSums(counts(object)) == 0
  }
  if (!"baseMean" %in% names(mcols(object))) {
    mcols(object)$baseMean <- rowMeans(counts(object,normalized=TRUE))
  }
  
  # make a design matrix with a term for every sample
  # this would typically produce unidentifiable solution
  # for the GLM, but we add priors for all terms except
  # the intercept
  samplesVector <- factor(samplesVector,levels=unique(samplesVector))
  if (missing(intercept)) {
    samples <- factor(c("null_level",as.character(samplesVector)),
                      levels=c("null_level",levels(samplesVector)))
    modelMatrix <- model.matrix(~samples)[-1,]
    modelMatrixNames <- colnames(modelMatrix)
    modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
  } else {
    # or we want to set the intercept using the
    # provided intercept instead
    samples <- factor(samplesVector)
    if (length(samples) > 1) {
      modelMatrix <- model.matrix(~ 0 + samples)
    } else {
      modelMatrix <- matrix(1,ncol=1)
      modelMatrixNames <- "samples1"
    }
    modelMatrixNames <- colnames(modelMatrix)
    if (!is.null(normalizationFactors(object))) { 
      nf <- normalizationFactors(object)
    } else {
      sf <- sizeFactors(object)
      nf <- matrix(rep(sf,each=nrow(object)),ncol=ncol(object))
    }
    # if the intercept is not finite, these rows
    # were all zero. here we put a small value instead
    intercept <- as.numeric(intercept)
    infiniteIntercept <- !is.finite(intercept)
    intercept[infiniteIntercept] <- -10
    normalizationFactors(object) <- nf * 2^intercept
    # we set the intercept, so replace the all zero
    # column with the rows which were all zero
    # in the previous dataset
    mcols(object)$allZero <- infiniteIntercept
  }
  
  # only continue on the rows with non-zero row sums
  objectNZ <- object[!mcols(object)$allZero,]
  stopifnot(all(!is.na(mcols(objectNZ)$dispFit)))
  
  # if a prior sigma squared not provided, estimate this
  # by the matching upper quantiles of the
  # log2 counts plus a pseudocount
  if (missing(betaPriorVar)) {
    logCounts <- log2(counts(objectNZ,normalized=TRUE) + 0.5)
    baseMean <- rowMeans(counts(objectNZ,normalized=TRUE))
    logFoldChangeMatrix <- logCounts - log2(baseMean + 0.5)
    logFoldChangeVector <- as.numeric(logFoldChangeMatrix)
    betaPriorVar <- matchUpperQuantileForVariance(logFoldChangeVector)
  }
  stopifnot(length(betaPriorVar)==1)
  
  lambda <- 1/rep(betaPriorVar,ncol(modelMatrix))
  # except for intercept which we set to wide prior
  if ("Intercept" %in% modelMatrixNames) {
    lambda[which(modelMatrixNames == "Intercept")] <- 1e-6
  }
  
  fit <- fitNbinomGLMs(object=objectNZ, modelMatrix=modelMatrix,
                       lambda=lambda, renameCols=FALSE,
                       alpha_hat=mcols(objectNZ)$dispFit,
                       betaTol=1e-4, useOptim=FALSE,
                       useQR=TRUE)
  normalizedDataNZ <- t(modelMatrix %*% t(fit$betaMatrix))
  
  normalizedData <- buildMatrixWithZeroRows(normalizedDataNZ, mcols(object)$allZero)
  
  # add back in the intercept, if finite
  if (!missing(intercept)) {
    normalizedData <- normalizedData + ifelse(infiniteIntercept, 0, intercept)
  }
  colnames(normalizedData) <- colnames(object)
  attr(normalizedData,"betaPriorVar") <- betaPriorVar
  if ("Intercept" %in% modelMatrixNames) {
    fittedInterceptNZ <- fit$betaMatrix[,which(modelMatrixNames == "Intercept"),drop=FALSE]
    fittedIntercept <- buildMatrixWithNARows(fittedInterceptNZ, mcols(object)$allZero)
    fittedIntercept[is.na(fittedIntercept)] <- -Inf
    attr(normalizedData,"intercept") <- fittedIntercept
  }
  normalizedData
}


rlog <- function(object, blind=TRUE, fast=FALSE,intercept, betaPriorVar, B) {
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  if (blind) {
    design(object) <- ~ 1
    object <- estimateDispersionsGeneEst(object, quiet=TRUE)
    object <- estimateDispersionsFit(object, quiet=TRUE)
  }
  if (is.null(mcols(object)$dispFit)) {
    object <- estimateDispersionsGeneEst(object, quiet=TRUE)
    object <- estimateDispersionsFit(object, quiet=TRUE)
  }
  if (!missing(intercept)) {
    if (length(intercept) != nrow(object)) {
      stop("intercept should be as long as the number of rows of object")
    }
  }
  if (!missing(B)) {
    if (length(B) != nrow(object)) {
      stop("B should be as long as the number of rows of object")
    }
    if (!all(B >= 0 & B <= 1)) {
      stop("B should be defined between 0 and 1")
    }
  }
  if (fast) {
    rld <- rlogDataFast(object, intercept, betaPriorVar, B)
  } else {
    rld <- rlogData(object, intercept, betaPriorVar)
  }
  se <- SummarizedExperiment(
    assays = rld,
    colData = colData(object),
    rowData = rowData(object),
    exptData = exptData(object))
  attr(se,"betaPriorVar") <- attr(rld, "betaPriorVar")
  if (!is.null(attr(rld,"intercept"))) {
    mcols(se)$rlogIntercept <- attr(rld,"intercept")
  }
  if (!is.null(attr(rld,"B"))) {
    attr(se,"B") <- attr(rld,"B")
  }
  se
}


# convenience function for testing the log likelihood
# for a count matrix, mu matrix and vector disp
nbinomLogLike <- function(counts, mu, disp) {
  rowSums(matrix(dnbinom(counts, mu=mu,size=1/disp, log=TRUE),ncol=ncol(counts)))
}


# Unexported, low-level function for fitting negative binomial GLMs
#
# Users typically call \code{\link{nbinomWaldTest}} or \code{\link{nbinomLRT}}
# which calls this function to perform fitting.  These functions return
# a \code{\link{DESeqDataSet}} object with the appropriate columns
# added.  This function returns results as a list.
#
# object a DESeqDataSet
# modelMatrix the design matrix
# modelFormula a formula specifying how to construct the design matrix
# alpha_hat the dispersion parameter estimates
# lambda the 'ridge' term added for the penalized GLM on the log2 scale
# renameCols whether to give columns variable_B_vs_A style names
# betaTol control parameter: stop when the following is satisfied:
#   abs(dev - dev_old)/(abs(dev) + 0.1) < betaTol
# maxit control parameter: maximum number of iteration to allow for
#   convergence
# useOptim whether to use optim on rows which have not converged:
#   Fisher scoring is not ideal with multiple groups and sparse
#   count distributions
# useQR whether to use the QR decomposition on the design matrix X
# forceOptim whether to use optim on all rows
# warnNonposVar whether to warn about non positive variances,
#   for advanced users only running LRT without beta prior,
#   this might be desirable to be ignored.
#
# return a list of results, with coefficients and standard
# errors on the log2 scale
fitNbinomGLMs <- function(file, useOptim=TRUE, forceOptim=FALSE, warnNonposVar=TRUE){
  # continue processing the results obtained from cluster
  load(file)
  
  modelMatrixNames <- colnames(modelMatrix)
  
  mu <- normalizationFactors * t(exp(modelMatrix %*% t(betaRes$beta_mat)))
  dispersionVector <- rep(dispersions(object), times=ncol(object))
  logLike <- nbinomLogLike(counts(object), mu, dispersions(object))
  
  # test for stability
  rowStable <- apply(betaRes$beta_mat,1,function(row) sum(is.na(row))) == 0
  
  # test for positive variances
  rowVarPositive <- apply(betaRes$beta_var_mat,1,function(row) sum(row <= 0)) == 0
  
  # test for convergence, stability and positive variances
  betaConv <- betaRes$iter < maxit
  
  # here we transform the betaMatrix and betaSE to a log2 scale
  betaMatrix <- log2(exp(1))*betaRes$beta_mat
  colnames(betaMatrix) <- modelMatrixNames
  colnames(modelMatrix) <- modelMatrixNames
  # warn below regarding these rows with negative variance
  betaSE <- log2(exp(1))*sqrt(pmax(betaRes$beta_var_mat,0))
  colnames(betaSE) <- paste0("SE_",modelMatrixNames)
  
  # switch based on whether we should also use optim
  # on rows which did not converge
  rowsForOptim <- if (useOptim) {
    which(!betaConv | !rowStable | !rowVarPositive)
  } else {
    which(!rowStable | !rowVarPositive)
  }
  
  if (forceOptim) {
    rowsForOptim <- seq_along(betaConv)
  }
  
  # not tested section
  if (length(rowsForOptim) > 0) {
    scaleCols <- apply(modelMatrix,2,function(z) max(abs(z)))
    x <- sweep(modelMatrix,2,scaleCols,"/")
    lambdaColScale <- lambda / scaleCols^2
    lambdaColScale <- ifelse(lambdaColScale == 0, 1e-6, lambdaColScale)
    lambdaLogScaleColScale <- lambdaLogScale / scaleCols^2
    large <- 30
    for (row in rowsForOptim) {
      betaRow <- if (rowStable[row]) {
        betaMatrix[row,] * scaleCols
      } else {
        beta_mat[row,] * scaleCols
      }
      betaRow <- pmin(pmax(betaRow, -large), large)
      nf <- normalizationFactors[row,]
      k <- counts(object)[row,]
      alpha <- alpha_hat[row]
      objectiveFn <- function(p) {
        mu_row <- as.numeric(nf * 2^(x %*% p))
        logLike <- sum(dnbinom(k,mu=mu_row,size=1/alpha,log=TRUE))
        logPrior <- sum(dnorm(p,0,sqrt(1/lambdaColScale),log=TRUE))
        -1 * (logLike + logPrior)
      }
      o <- optim(betaRow, objectiveFn, method="L-BFGS-B",lower=-30, upper=30)
      ridge <- if (length(lambdaLogScale) > 1) {
        diag(lambdaLogScaleColScale)
      } else {
        as.matrix(lambdaLogScaleColScale,ncol=1)
      }
      # if we converged, change betaConv to TRUE
      if (o$convergence == 0) {
        betaConv[row] <- TRUE
      }
      # with or without convergence, store the estimate from optim
      betaMatrix[row,] <- o$par / scaleCols
      # calculate the standard errors
      mu_row <- as.numeric(nf * 2^(x %*% o$par))
      w <- diag((mu_row^-1 + alpha)^-1)
      xtwx <- t(x) %*% w %*% x
      xtwxRidgeInv <- solve(xtwx + ridge)
      sigma <- xtwxRidgeInv %*% xtwx %*% xtwxRidgeInv
      # warn below regarding these rows with negative variance
      betaSE[row,] <- log2(exp(1)) * sqrt(pmax(diag(sigma),0)) / scaleCols
      # store the new mu vector
      mu[row,] <- mu_row
      logLike[row] <- sum(dnbinom(k, mu=mu_row, size=1/alpha, log=TRUE))
    }
  }
  
  nNonposVar <- sum(rowSums(betaSE == 0) > 0)
  if (warnNonposVar & nNonposVar > 0) warning(nNonposVar,"rows had non-positive estimates of variance for coefficients, likely due to rank deficient model matrices without betaPrior")
  
  list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
       betaSE = betaSE, mu = mu, betaIter = betaRes$iter,
       deviance = betaRes$deviance, modelMatrix=modelMatrix, 
       nterms=ncol(modelMatrix), hat_diagonals=betaRes$hat_diagonals)
}


# convenience function for building larger matrices
# by filling in 0 rows
buildMatrixWithZeroRows <- function(m, zeroRows) {
  mFull <- matrix(0, ncol=ncol(m), nrow=length(zeroRows))
  mFull[!zeroRows,] <- m
  mFull
}


# convenience function for building larger matrices
# by filling in NA rows
buildMatrixWithNARows <- function(m, NARows) {
  mFull <- matrix(NA, ncol=ncol(m), nrow=length(NARows))
  mFull[!NARows,] <- m
  mFull
}


# looking at the values of x which are large
# in absolute value, find the zero-centered Normal distribution
# with the matching quantile, and return the variance
# of that Normal distribution
matchUpperQuantileForVariance <- function(x, upperQuantile=.05) {
  sdEst <- quantile(abs(x), 1 - upperQuantile) / qnorm(1 - upperQuantile/2)
  unname(sdEst)^2
}


rlogData <- function(file, object, intercept, betaPriorVar) {
  # only continue on the rows with non-zero row sums
  objectNZ <- object[!mcols(object)$allZero,]
  stopifnot(all(!is.na(mcols(objectNZ)$dispFit)))
  
  # if a prior sigma squared not provided, estimate this
  # by the matching upper quantiles of the
  # log2 counts plus a pseudocount
  if (missing(betaPriorVar)) {
    logCounts <- log2(counts(objectNZ,normalized=TRUE) + 0.5)
    baseMean <- rowMeans(counts(objectNZ,normalized=TRUE))
    logFoldChangeMatrix <- logCounts - log2(baseMean + 0.5)
    logFoldChangeVector <- as.numeric(logFoldChangeMatrix)
    betaPriorVar <- matchUpperQuantileForVariance(logFoldChangeVector)
  }
  stopifnot(length(betaPriorVar)==1)
  
  fit <- fitNbinomGLMs(file, useOptim=FALSE)
  normalizedDataNZ <- t(modelMatrix %*% t(fit$betaMatrix))
  
  normalizedData <- buildMatrixWithZeroRows(normalizedDataNZ, mcols(object)$allZero)
  
  # add back in the intercept, if finite
  if (!missing(intercept)) {
    normalizedData <- normalizedData + ifelse(infiniteIntercept, 0, intercept)
  }
  colnames(normalizedData) <- colnames(object)
  attr(normalizedData,"betaPriorVar") <- betaPriorVar
  if ("Intercept" %in% modelMatrixNames) {
    fittedInterceptNZ <- fit$betaMatrix[,which(modelMatrixNames == "Intercept"),drop=FALSE]
    fittedIntercept <- buildMatrixWithNARows(fittedInterceptNZ, mcols(object)$allZero)
    fittedIntercept[is.na(fittedIntercept)] <- -Inf
    attr(normalizedData,"intercept") <- fittedIntercept
  }
  normalizedData
}


# Apply a 'regularized log' transformation
rlog <- function(file, object, blind=TRUE) {
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  if (blind) {
    design(object) <- ~ 1
    object <- estimateDispersionsGeneEst(object, quiet=TRUE)
    object <- estimateDispersionsFit(object, quiet=TRUE)
  }
  if (is.null(mcols(object)$dispFit)) {
    object <- estimateDispersionsGeneEst(object, quiet=TRUE)
    object <- estimateDispersionsFit(object, quiet=TRUE)
  }
  if (!missing(intercept)) {
    if (length(intercept) != nrow(object)) {
      stop("intercept should be as long as the number of rows of object")
    }
  }
  if (!missing(B)) {
    if (length(B) != nrow(object)) {
      stop("B should be as long as the number of rows of object")
    }
    if (!all(B >= 0 & B <= 1)) {
      stop("B should be defined between 0 and 1")
    }
  }

  rld <- rlogData(file, object, intercept, betaPriorVar)
  se <- SummarizedExperiment(
    assays = rld,
    colData = colData(object),
    rowData = rowData(object),
    exptData = exptData(object))
  attr(se,"betaPriorVar") <- attr(rld, "betaPriorVar")
  if (!is.null(attr(rld,"intercept"))) {
    mcols(se)$rlogIntercept <- attr(rld,"intercept")
  }
  if (!is.null(attr(rld,"B"))) {
    attr(se,"B") <- attr(rld,"B")
  }
  se
}