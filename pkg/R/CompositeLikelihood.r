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

fExtremal <- function(coord, corrmodel, data, fixed=list(nugget=0), hessian=TRUE,
                      method='L-BFGS-B', model, parscale=TRUE, start=NULL)
{
    call <- match.call()
    
    ### Check the parameters given in input:

    if(missing(model) || is.null(model))
       stop('The model inserted is not an Extremal process available\n')
    
    if(missing(corrmodel) || is.null(model))
      stop('the model paramater/s need to be completely inserted\n')
 
    if(missing(coord) || !is.matrix(coord))
      stop('insert a suitable set of coordinates\n')
     
    if(missing(data) || !is.matrix(data))
      stop('insert a suitable matrix of data\n')
  
    if(!is.null(start) & !is.list(start))
      stop('insert starting values as a list of parameters\n')
   
    if(!is.null(fixed) & !is.list(fixed))
      stop('insert fixed values as a list of parameters\n')
 
    ### Initialization global variables:
     
    fExtremal <- NULL
    namemodel <- model
    namecorrmodel <- corrmodel
    model <- CheckModel(model)
    corrmodel <- CheckCorrModel(corrmodel)

  
    ### Initialization parameters:
    
    initparam <- InitParam(namecorrmodel, fixed, namemodel, parscale,
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

    ### Some checks of the output from the optimization procedure:

    if(fitted$convergence == 0)
      fitted$convergence <- 'Successful'
    else
      if(fitted$convergence == 1)
        fitted$convergence <- 'Iteration limit reached'
      else
        fitted$convergence <- "Optimization may have failed"


    ### Start computation of the variance-covariance matrix:

    param <- c(fitted$par, initparam$fixed)
    param <- param[sort(names(param))]
    eps <- (.Machine$double.eps)^(1/3)
    numflag <- initparam$numparam + initparam$numfixed
    dimvar <- initparam$numparam^2
    vsmat <- double(2 * dimvar)
    
    .C('SquaredScore', as.integer(corrmodel), as.double(unifdata),
       as.double(eps), as.integer(initparam$flag), as.double(lags),
       as.integer(model), as.integer(numdata), as.integer(numflag),
       as.integer(initparam$numparam), as.integer(numcoord),
       as.double(param), vsmat, PACKAGE='ExtremalProc',
       DUP = FALSE, NAOK=TRUE)

    # Set sensitivity matrix:
    sensmat <- matrix(vsmat[1 : dimvar], ncol=initparam$numparam)
    # Set variability matrix:
    varmat <- matrix(vsmat[(dimvar + 1) : (2 * dimvar)],
                        ncol=initparam$numparam)

    if(hessian)
      sensmat <- -fitted$hessian
    
    isensmat <- try(solve(sensmat), silent = TRUE)
    
    if(!is.matrix(isensmat) || !is.matrix(varmat))
      {
        warning("observed information matrix is singular")
        varcov <- 'none'
        std.err <- 'none'
      }
    else
      {
        penalty <- varmat %*% isensmat
        CLIC <- -2 * (fitted$value - sum(diag(penalty)))
        varcov <- isensmat %*% penalty
        dimnames(varcov) <- dimnames(sensmat) <- dimnames(varmat) <- list(initparam$namesparam, initparam$namesparam)
        std.err <- diag(varcov)
        if(any(std.err < 0))
          std.err <- 'none'
        else
          {
            std.err <- sqrt(std.err)
            names(std.err) <- initparam$namesparam
          }
      }

    ### End computation of the variance covariance matrix

     
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
                      sensmat = sensmat,
                      std.err = std.err,
                      unifdata = unifdata,
                      varcov = varcov,
                      varmat = varmat)

    structure(c(fExtremal, call = call), class = c("extremal"))
  }

### Summary of fitted models:

print.extremal <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    numdata <- nrow(x$data)
    numcoord <- ncol(x$data)
    numparam <- length(x$param)

    cat('\nResults of fitting Extremal processes.\n')
    cat('\nModel: ', x$model, '\n')
    cat('Covariance model: ', x$corrmodel, '\n')
    cat('Number of coordinates: ', numcoord, '\n')
    cat('Number of observation per location: ', numdata, '\n\n')
    cat('Log composite likelihood:      ', x$logLik, '\n')
    cat('Composite likelihood criterion: ', x$CLIC,'\n')
    cat('\nEstimated parameters:\n')
    print.default(format(x$param, digits = digits), print.gap = 2, 
        quote = FALSE)
    cat('\nStandard errors:\n')
    print.default(format(x$std.err, digits = digits), print.gap = 2, 
        quote = FALSE)
    cat('\nVariance-covariance matrix of the estimates:\n')
    print.default(format(x$varcov, digits = digits), print.gap = 3,
                  quote = FALSE)
    invisible(x)
  }

