####################################################
### Authors: Simone Padoan.
### Email: simone.padoan@eofl.ch.
### Institute: EPFL.
### File name: CompositeLikelihood.r
### Description:
### This file contains a set of procedures
### for fitting extremal models.
### Last change: 20/01/2010.
####################################################


### Procedures are in alphabetical order.



### Composite log-likelihood procedure for Extremal processes

CompLikelihood <- function(corrmodel, data, fixed, lags, model, numcoord, numdata, param)
  {
    result <- double(1)
    result <- -1.0e8

    if(!CheckParamRange(param))
      return(result)

    param <- c(param, fixed)
    param <- param[sort(names(param))]

    result <- .C('CompLikelihood', as.integer(corrmodel), as.double(data), as.double(lags),
                 as.integer(model), as.integer(numdata), as.integer(numcoord), as.double(param),
                 res=double(1), PACKAGE='ExtremalProc', DUP = FALSE, NAOK=TRUE)$res
    
    return(result)
  }


### Fit multivariate extremal models 

fExtremal <- function(coord=NULL, corrmodel=NULL, data=NULL, fixed=list(nugget=0),
                      method='L-BFGS-B', model=NULL, parscale=TRUE, start=NULL)
  {
    ### Initialization global variables:
    fExtremal <- NULL
    namemodel <- model
    namecorrmodel <- corrmodel
    model <- CheckModel(model)
    corrmodel <- CheckCorrModel(corrmodel)

    ### Check the parameters given in input:
     if(is.null(model))
      {
        cat('The model inserted is not a multivariate spatial model available\n')
        return(fExtremal)
      }
    
    if(is.null(corrmodel))
      {
        cat('the model paramater/s need to be completely inserted\n')
        return(fExtremal)
      }

    if(is.null(coord) || !is.matrix(coord))
      {
        cat('insert a suitable set of coordinates\n')
        return(fExtremal)
      }
    
    if(is.null(data) || !is.matrix(data))
      {
        cat('insert a suitable matrix of data\n')
        return(fExtremal)
      }

    if(!is.null(start) & !is.list(start))
      {
        cat('insert starting values as a list of parameters\n')
        return(fExtremal)
      }

    if(!is.null(start) & !is.list(fixed))
      {
        cat('insert fixed values as a list of parameters\n')
        return(fExtremal)
      }

    ### Initialization parameters:
    
    initparam <- InitParam(fixed, namecorrmodel, namemodel, parscale,
                           method=='L-BFGS-B', start)

    if(!is.null(initparam$error))
      {
        cat(initparam$error)
        return(fExtremal)
      }

    ### Set further global variables:
    
    numcoord <- ncol(data)       # number of coordinates
    numdata <- nrow(data)        # number of observations
    lags <- c(dist(coord))       # distances between pair of sites


    ### Transform data to uniform marginals:
    unifdata <- Gev2Unif(data)
    
    ### Model fitting with uniform marginals:
    if(method=='L-BFGS-B')
      fitted <- optim(initparam$param, CompLikelihood, corrmodel=corrmodel, data=unifdata,
                      fixed=initparam$fixed, lags=lags, method=method, model=model,
                      numcoord=numcoord, numdata=numdata, control=list(fnscale=-1, factr=1,
                      pgtol=1e-14, maxit = 1e8, parscale=initparam$parscale),
                      lower=initparam$lower, upper=initparam$upper, hessian=TRUE)
    else
      fitted <- optim(initparam$param, CompLikelihood, corrmodel=corrmodel, data=unifdata,
                      fixed=initparam$fixed, lags=lags, method=method, model=model,
                      numcoord=numcoord, numdata=numdata, control=list(fnscale=-1,
                      reltol=1e-14, maxit=1e8, parscale=initparam$parscale), hessian=TRUE)

### Start computation of the variance-covariance matrix:

    param <- c(fitted$par, initparam$fixed)
    param <- param[sort(names(param))]
    eps <- (.Machine$double.eps)^(1/3)
    numflag <- initparam$numparam + initparam$numfixed
    
    squaredscore <- .C('SquaredScore', as.integer(corrmodel), as.double(unifdata),
                       as.double(eps), as.integer(initparam$flag), as.double(lags),
                       as.integer(model), as.integer(numdata), as.integer(numflag),
                       as.integer(numcoord), as.integer(initparam$numparam),
                       as.double(param), sqscore=double(length=initparam$numparam^2),
                       PACKAGE='ExtremalProc', DUP = FALSE, NAOK=TRUE)$sqscore

    squaredscore <- matrix(squaredscore, ncol=initparam$numparam,
                           nrow=initparam$numparam)
    ### End computation of the variance covariance matrix
    
    ihessian <- try(solve(fitted$hessian), silent = TRUE)
    
    if(fitted$convergence == 0)
      fitted$convergence <- 'Successful'
    else
      if(fitted$convergence == 1)
        fitted$convergence <- 'Iteration limit reached'
      else
        fitted$convergence <- "Optimization may have failed"

    if(!is.matrix(ihessian) || !is.matrix(squaredscore))
      {
        warning("observed information matrix is singular")
        varcov <- 'none'
        std.err <- 'none'
      }
    else
      {
        penalty <- squaredscore %*% ihessian
        CLIC <- -2 * (fitted$value - sum(diag(penalty)))
        varcov <- ihessian %*% penalty
        dimnames(varcov) <- list(initparam$namesparam, initparam$namesparam)
        std.err <- diag(varcov)
        if(any(std.err < 0))
          std.err <- 'none'
        else
          {
            std.err <- sqrt(std.err)
            names(std.err) <- initparam$namesparam
          }
      }
    
    fExtremal <- list(CLIC = CLIC,
                      coord = coord,
                      convergence = fitted$convergence,
                      corrmodel = namecorrmodel,
                      data = data,
                      fixed = initparam$fixed,
                      iterations = fitted$counts,
                      logLik = fitted$value,
                      message = fitted$message,
                      model = namemodel,
                      param = fitted$par,
                      std.err = std.err,
                      varcov = varcov)

    class(fExtremal) <- "Extremal"

    return(fExtremal)
  }

### Summary of fitted models:

summary.Extremal <- function(fitted)
  {
    numdata <- nrow(fitted$data)
    numcoord <- ncol(fitted$data)
    numparam <- length(fitted$param)

    cat('\nResults of fitting Extremal processes.\n')
    cat('\nModel: ', fitted$model, '\n')
    cat('Covariance model: ', fitted$corrmodel, '\n')
    cat('Number of coordinates: ', numcoord, '\n')
    cat('Number of observation per location: ', numdata, '\n\n')
    cat('Log composite likelihood:      ', fitted$logLik, '\n')
    cat('Composite likelihood criterion: ', fitted$CLIC,'\n')
    cat('\nEstimated parameters:\n')
    print(fitted$param)
    cat('\nStandard errors:\n')
    print(fitted$std.err)
    cat('\nVariance-covariance matrix of the estimates:\n')
    print(fitted$varcov)

  }

