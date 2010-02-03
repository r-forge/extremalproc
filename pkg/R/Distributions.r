####################################################
### Authors: Simone Padoan.
### Email: simone.padoan@eofl.ch.
### Institute: EPFL.
### File name: Distribution.r
### Description:
### This file contains a set of procedures
### that define several distributions and
### related functions.
### Last change: 20/01/2010.
####################################################


### Procedures are in alphabetical order.

fGev <- function(data, method='Nelder-Mead', std.err=FALSE)
  {
    ### Starting values:
    start <- c(1, 1, .01)

    ### fit the log-likelihood of the GEV model:
    fitted <- optim(start, GevLoglik, data=data, method=method,
                    control=list(fnscale=-1, reltol=1e-14, maxit=1e8),
                    hessian=std.err)
    return(fitted)
  }


Gev2Unif <- function(data)
  {
    dimdata <- dim(data)
    
    if(length(dimdata)==2)
      {
        numdata <- dimdata[1]
        numcoord <- dimdata[2]
        est <- apply(data, 2, MLE)
      }

    if(length(dimdata)==3)
      {
        numdata <- dimdata[3]
        numcoord <- prod(dimdata[1:2])
      }

    result <- .C('FromDistToDist', as.double(data), as.double(rep(0, numcoord)),
                 as.double(rep(1, numcoord)), as.double(rep(1, numcoord)),
                 as.double(est[1,]), as.integer(numdata), as.integer(numcoord),
                 as.double(est[2,]), as.double(est[3,]), as.integer(1),
                 res=double(numdata * numcoord), PACKAGE='ExtremalProc',
                 DUP = FALSE, NAOK=TRUE)$res
    
    dim(result) <- dimdata
    
    return(result)
  }

GevLoglik <- function(data, param)
  {
    ### Initialization variables:
    loc <- param[1]
    numdata <- length(data)
    scale <- param[2]
    shape <- param[3]

    ### Compute the log-likelihood:
    result <- .C('GevLogLik', as.double(data), as.double(loc),
                 as.integer(numdata), as.double(scale),
                 as.double(shape), res=double(1),
                 PACKAGE='ExtremalProc', DUP = FALSE,
                 NAOK=TRUE)$res

    return(result)
  }

MLE <- function(x)
  return(as.double(fGev(x)$par))
