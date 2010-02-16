####################################################
### Authors: Simone Padoan.
### Email: simone.padoan@eofl.ch.
### Institute: EPFL.
### File name: Statistics.r
### Description:
### This file contains a set of procedures
### for the computation of statistics useful for
### statistical hypothesis testing.
### Last change: 10/02/2010.
####################################################


### Procedures are in alphabetical order.



### Statistical hypothesis testing 

HypoTest <- function(object1, object2, ..., statistic)
  {
    ### Check input parameters:
    if(any(missing(object1), missing(object2)))
      stop('Models one and two must be specified')

    if(missing(statistic))
      stop('Insert the type of statistic use in the hypothesis test')

    objects <- as.list(substitute(list(...)))[-1]
    objects <- sapply(objects,function(x) deparse(x))
    if(!length(objects)) objects <- NULL

    models <- c(deparse(substitute(object1)),
                deparse(substitute(object2)), objects)

    nummod <- length(models)
    numpar <- NULL
    numst <- length(statistic)
    lmodels <- vector('list', nummod)
    W <- pvalue <- df <- nu <- double(nummod - 1)
       
    ### Perform the test:
    for(i in 1:nummod)
      {
        model <- get(models[i], envir = parent.frame())
        if(!inherits(model, "extremal"))
          stop("Use only with 'extremal' objects")

        numpar <- c(numpar, length(model$param))
        lmodels[[i]] <- model
        if(i > 1)
          {
            j <- i - 1
            if((!all(names(lmodels[[j]]) %in% names(lmodels[[i]]))) &&
               (!identical(lmodels[[j]]$model, lmodels[[i]]$model)) &&
               (!TestCorr(lmodels[[j]]$corrmodel, lmodels[[i]]$corrmodel)))
              stop('models are not nested')

            df[j] <- length(lmodels[[j]]$param) - length(lmodels[[i]]$param)
            if(df[j] <= 0)
              stop('model are not nested')
            nu[j] <- df[j]
            namespar <- names(lmodels[[j]]$param)[!names(lmodels[[j]]$param) %in%
                                                      names(lmodels[[i]]$param)]

            ### Composite likelihood ratio statistic (or Wilks test):
            if(grepl('Wilks', statistic))
              {
                W[j] <- 2 * (lmodels[[j]]$logLik - lmodels[[i]]$logLik)
                if(W[j] < 0)
                  stop('negative value of the statistic')
                varcov <- lmodels[[j]]$varcov
                iH <- solve(lmodels[[j]]$sensmat)
                lambda <- eigen(solve(iH[namespar, namespar]) %*%
                                varcov[namespar, namespar])$values
                if(grepl('-RJ', statistic))
                  W[j] <- W[j] / mean(lambda)
                if(grepl('-S', statistic))
                  {
                    nu[j] <- sum(lambda) / sum(lambda^2)
                    W[j] <- nu[j] * W[j] / (df[j] * mean(lambda))
                  }
              }
            ### Wald-type statistic:
            if(grepl('Wald', statistic))
              W[j] <- t(lmodels[[j]]$param[namespar] - lmodels[[i]]$fixed[namespar]) %*%
                solve(lmodels[[j]]$varcov)[namespar, namespar] %*%
                  (lmodels[[j]]$param[namespar] - lmodels[[i]]$fixed[namespar])
            
            ### Score-type statistic (or Rao test):
            if(grepl('Rao', statistic))
              {
                corrmodel <- CheckCorrModel(lmodels[[j]]$corrmodel)
                param <- c(lmodels[[i]]$param, lmodels[[i]]$fixed)
                param <- param[sort(names(param))]
                numflag <- length(lmodels[[i]]$param) + length(lmodels[[i]]$fixed)
                eps <- (.Machine$double.eps)^(1/3)
                flag <- rep(0, numflag)
                flag[names(param) %in% namespar] <- 1
                gradient <- double(length(namespar))
                lags <- dist(lmodels[[j]]$coord)
                model <- CheckModel(lmodels[[j]]$model)
                numdata <- nrow(lmodels[[j]]$data)
                numcoord <- nrow(lmodels[[j]]$coord)
                
                .C('Gradient', as.integer(corrmodel), as.double(lmodels[[j]]$unifdata), as.double(eps),
                   as.integer(flag), as.double(lags), as.integer(model), as.integer(numdata),
                   as.integer(length(flag)), as.integer(length(namespar)), as.integer(numcoord),
                   as.double(param), gradient, PACKAGE='ExtremalProc', DUP = FALSE, NAOK=TRUE)

                flag <- rep(0, numflag)
                flag[names(param) %in% names(lmodels[[j]]$param)] <- 1
                numpar <- length(lmodels[[j]]$param)
                dimvar <- numpar^2
                vsmat <- double(2 * dimvar)
                
                .C('SquaredScore', as.integer(corrmodel), as.double(lmodels[[j]]$unifdata), as.double(eps),
                   as.integer(flag), as.double(lags), as.integer(model), as.integer(numdata),
                   as.integer(length(flag)), as.integer(numpar), as.integer(numcoord),
                   as.double(param), vsmat, PACKAGE='ExtremalProc', DUP = FALSE, NAOK=TRUE)

                # Set sensitivity matrix:
                sensmat <- matrix(vsmat[1 : dimvar], ncol=numpar)
                # Set variability matrix:
                varmat <- matrix(vsmat[(dimvar + 1) : (2 * dimvar)], ncol=numpar)
                dimnames(sensmat) <- dimnames(varmat) <- list(names(lmodels[[j]]$param), names(lmodels[[j]]$param))
                # Godambe matrix:
                G <- sensmat %*% varmat %*% sensmat
                
                W[j] <- t(gradient) %*% solve(sensmat)[namespar, namespar]  %*% solve(solve(G)[namespar, namespar]) %*%
                  solve(sensmat)[namespar, namespar] %*% gradient
              }

            pvalue[j] <- pchisq(W[j], df = nu[j], lower.tail = FALSE)
          }
      }

    table <- data.frame(numpar, c(NA, df), c(NA, nu),c(NA, W), c(NA, pvalue))
    dimnames(table) <- list(models, c("Num.Par", "Diff.Par", "Df","Chisq", "Pr(>chisq)"))
    structure(table, heading = c("Statistical Hypothesis Test Table\n"),
              class = c("data.frame"))

    
  }

TestCorr <- function(model1, model2)
  {
    result <- FALSE

    if(identical(model1, model2))
      result <- TRUE

    if(identical(model1, 'gencauchy'))
      if(identical(model2, 'cauchy'))
        result <- TRUE

    if(identical(model1, 'stable'))
      if(identical(model2, 'gauss') ||
         identical(model2, 'exponential'))
        result <- TRUE

    return(result)
  }
