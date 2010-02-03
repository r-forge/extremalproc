####################################################
### Authors: Simone Padoan.
### Email: simone.padoan@eofl.ch.
### Institute: EPFL.
### File name: Utility.r
### Description:
### This file contains a set of procedures
### for supporting all the other functions.
### Last change: 20/01/2010.
####################################################



### Procedures are in alphabetical order.

CheckCorrModel <- function(corrmodel)
  {
    CheckCorrModel <- NULL
    # Correlation function are in alphabetical order
    CheckCorrModel <- switch(corrmodel,
                      cauchy=1,
                      exponential=2,
                      gauss=3,
                      gencauchy=4,
                      stable=5,
                      whittlematern=6)

    
    return(CheckCorrModel)
  }

CheckModel <- function(model)
  {
    CheckModel <- NULL
    
    CheckModel <- switch(model,
                         Extremal_g=1,
                         Extremal_t=2)

    return(CheckModel)   
  }

CheckParam <- function(corrmodel, namesparam, numparam)
  {
    for(i in 1 : numparam)
      {
        if(corrmodel=='exponential' || corrmodel=='gauss')
          if(is.null(switch(namesparam[i],
                            df=1,
                            nugget=2,
                            scale=3)))
            return(FALSE)

        if(corrmodel=='stable')
          if(is.null(switch(namesparam[i],
                            df=1,
                            nugget=2,
                            power=3,
                            scale=4)))
            return(FALSE)

        if(corrmodel=='cauchy')
          if(is.null(switch(namesparam[i],
                            df=1,
                            nugget=2,
                            power2=3,
                            scale=4)))
            return(FALSE)
        
        if(corrmodel=='gencauchy')
          if(is.null(switch(namesparam[i],
                            df=1,
                            nugget=2,
                            power1=3,
                            power2=4,
                            scale=5)))
            return(FALSE)

        if(corrmodel=='whittlematern')
          if(is.null(switch(namesparam[i],
                            df=1,
                            nugget=2,
                            scale=3,
                            smooth=4)))
            return(FALSE)
      }

    return(TRUE)
  }

CheckParamRange <- function(param)
  {

    if(!is.na(param['df'])) if(param['df'] <= 0) return(FALSE)
    if(!is.na(param['nugget'])) if(param['nugget'] < 0 || param['nugget'] > 1) return(FALSE)
    if(!is.na(param['power'])) if(param['power'] <=0 || param['power'] > 2) return(FALSE)
    if(!is.na(param['power1'])) if(param['power1'] <=0 || param['power1'] > 2) return(FALSE)
    if(!is.na(param['power2'])) if(param['power2'] <= 0) return(FALSE)
    if(!is.na(param['scale'])) if(param['scale'] <= 0) return(FALSE)
    if(!is.na(param['smooth'])) if(param['smooth'] <= 0) return(FALSE)

    return(TRUE)
  }


InitParam <- function(fixed, corrmodel, model, parscale, paramrange, start)
  {    
    ### Parameters initalization:
    error <- NULL
    df <- smooth <- 1
    nugget <- .1
    scale <- 10
    numfixed <- numstart <- 0

    if(corrmodel=='cauchy')
      {
        param <- c(nugget, smooth, scale)
        namesparam <- c('nugget', 'power2', 'scale')
      }
    
    if(corrmodel=='exponential' || corrmodel=='gauss')
      {
        param <- c(nugget, scale)
        namesparam <- c('nugget', 'scale')
      }

    if(corrmodel=='gencauchy')
      {
        param <- c(nugget, smooth, smooth, scale)
        namesparam <- c('nugget', 'power1', 'power2','scale')

        if(model=='Extremal_g')
          if(is.list(fixed)) fixed$power2 <- 1 else fixed <- list(power2=1)
      }

    if(corrmodel=='stable')
      {
        param <- c(nugget, smooth, scale)
        namesparam <- c('nugget', 'power', 'scale')
      }

    if(corrmodel=='whittlematern')
      {
        param <- c(nugget, scale, smooth)
        namesparam <- c('nugget', 'scale', 'smooth')
      }

    param <- c(df, param)
    namesparam <- c('df', namesparam)
    names(param) <- namesparam

    if(model=='Extremal_g')
      if(is.list(fixed)) fixed$df <- df else fixed <- list(df=df)

    numparam <- length(param)
    flag <- rep(1, numparam)
    namesflag <- namesparam
    
    ### first check if the starting parameters are
    ### corrects and update then the  parameter vector
    
    if(!is.null(start))
      {
        namesstart <- names(start)
        numstart <- length(namesstart)
        start <- as.numeric(start)
        names(start) <- namesstart
        
        if(!CheckParam(corrmodel, namesstart, numstart))
          {
            error <- 'some names of the starting parameters is/are not correct\n'
            return(list(error=error))
          }
        
        if(!CheckParamRange(start))
          {
            error <- 'some starting values are out of the range\n'
            return(list(error=error))
          }

        for(i in 1 : numstart)
          param[namesstart[i]] <- start[namesstart[i]]
      }

    ### first check if the fixed parameters are
    ### corrects and update then the  parameter vector

    if(!is.null(fixed))
      {
        namesfixed <- names(fixed)
        numfixed <- length(namesfixed)
        fixed <- as.numeric(fixed)
        names(fixed) <- namesfixed
     
        if(!CheckParam(corrmodel, namesfixed, numfixed))
          {
            error <- 'some names of the fixed parameters is/are not correct\n'
            return(list(error=error))
          }
        
        if(!CheckParamRange(fixed))
          {
            error <- 'some fixed values are out of the range\n'
            return(list(error=error))
          }

        if(numfixed==numparam)
          {
            error <- 'the are not parameters left to estimate\n'
            return(list(error=error))
          }
        
        for(i in 1 : numfixed)
          {
            flag[namesflag==namesfixed[i]] <- 0
            param <- param[!namesparam==namesfixed[i]]
            namesparam <- names(param)
          }
        numparam <- length(param)
      }

    ### check the consistence between fixed and starting values
    
    if(numstart > 0 && numfixed > 0)
      for(i in 1 : numstart)
        for(j in 1 : numfixed)
          if(namesstart[i]==namesfixed[j])
            {
              error <- ('some fixed parameter name/s is/are matching with starting parameter name/s\n')
              return(list(error=error))
            }

    ### set the scale of the parameters:

    
    if(is.logical(parscale))
      parscale <- SetParamScale(numparam, param, parscale, namesparam)
    else
      {
        if(!is.list(parscale))
          {
            error <- ('the scale of the parameters need to be a list\n')
            return(list(error=error))
          }
        else
          {
            namesparscale <- names(parscale)
            numparscale <- length(namesparscale)
            
            if(!CheckParam(corrmodel, namesparscale, numparscale))
              {
                error <- 'some names of the parameters scales are not correct\n'
                return(list(error=error))
              }
            parscale <- SetParamScale(numparam, param, parscale, namesparam)
          }
      }

    ### set the range of the parameters if its the case

    if(paramrange)
      paramrange <- SetRangeParam(namesparam, numparam)
    else
      paramrange <- list(lower=NULL, upper=NULL)
    

    return(list(error=error, flag=flag, fixed=fixed, lower=paramrange$lower,
                namesparam=namesparam, numparam=numparam, numfixed=numfixed,
                param=param, parscale=parscale, upper=paramrange$upper))
  }

InitSimu <- function(corrmodel, model, numblock, numsim, param, coord)
  {
    error <- NULL
    namecorr <- corrmodel

    model <- CheckModel(model)

    if(is.null(model))
      {
        error <- ('The name of the process is not correct\n')
        return(list(error=error))
      }

    if(is.null(numsim))
      {
        error <- ('the sample size need to be a positive integer\n')
        return(list(error=error))
      }
    
    if(!is.null(numsim))
      if(any(numsim < 0, !is.numeric(numsim)))
        {
          error <- ('the sample size need to be a positive integer\n')
          return(list(error=error))
        }

    if(is.null(coord) || !is.matrix(coord))
      {
        error <- ('insert a suitable set of coordinates\n')
        return(list(error=error))
      }
 
    if(is.null(numblock))
      {
        error <- ('the block size need to be a positive integer\n')
        return(list(error=error))
      }

    if(!is.null(numblock))
      if(any(numblock < 0, !is.numeric(numblock)))
        {
          error <- ('the block size need to be a positive integer\n')
          return(list(error=error))
        }

    if(!is.list(param))
      {
        error <- ('model parameters need to be passed as list')
        return(list(error=error))
      }

     if(is.null(param$nugget))
      {
        error <- ('the nugget parameter is required for the model definition, please insert it\n')
        return(list(error=error))
      }
        
    if(is.null(param$scale))
      {
        error <- ('the scale parameter is required from all the correlation models, please insert it\n')
        return(list(error=error))
      }

    if(model == 2 && is.null(param$df))
      {
        error <- ('insert the degree of freedom of the Student t distribution\n')
        return(list(error=error))
      }

    corrmodel <- CheckCorrModel(corrmodel)
    
    if(is.null(corrmodel))
      {
        error <- ('Insert a valid correlation model\n')
        return(list(error=error))
      }
    
    if(is.null(param$mean))
      param$mean <- 0

    if(model == 1)
      param$df <- 1

    param$variance <- 1 - param$nugget

    namesparam <- names(param)
    param <- as.numeric(param)
    names(param) <- namesparam
    namesparam <- grep('mean', namesparam, value=TRUE, invert=TRUE)
    namesparam <- grep('variance', namesparam, value=TRUE, invert=TRUE)
    numparam <- length(namesparam)

    if(!CheckParam(namecorr, namesparam, numparam))
      {
        error <- ('paramater/s of the correlation model need to be completely inserted\n')
        return(list(error))
      }

    namesparam <- grep('df', namesparam, value=TRUE, invert=TRUE)
    namesparam <- grep('nugget', namesparam, value=TRUE, invert=TRUE)
    namesparam <- grep('scale', namesparam, value=TRUE, invert=TRUE)
    namesparam <- sort(namesparam)
    
    corrparam <- c(param['mean'], param['variance'], 0, param['scale'], param[namesparam])


    return(list(corrmodel=corrmodel, corrparam=corrparam, error=error,
                model=model, param=param))
  }


RescalingCoord <- function(corrmodel, numblock, param, coord)
  {
    coordscaled <- NULL
    
    if(corrmodel=='exponential')
      coordscaled <- coord / log(numblock)

    if(corrmodel=='gauss')
      coordscaled <- coord / sqrt(log(numblock))
    
    if(corrmodel=='stable')
      coordscaled <- coord / (log(numblock))^(1 / param[5])
 
    if(corrmodel=='gencauchy')
      coordscaled <- coord / (log(numblock))^(1 / param[5])
    
    return(coordscaled)
  }

    
SetParamScale <- function(numparam, param, parscale, namesparam)
  {    
    # default values:  

    scales <- c(100, 1, 2, 2, 10, 300, 10)
    names(scales) <- c('df', 'nugget', 'power', 'power1', 'power2', 'scale', 'smooth')
    
    if(is.logical(parscale))
      {
        if(!parscale)
          scales <- rep(1, 7)
      }
    else
      {
        # set given values:

        if(!is.null(parscale$df))
          scales['df'] <- parscale$nugget

        if(!is.null(parscale$nugget))
          scales['nugget'] <- parscale$nugget

        if(!is.null(parscale$power2))
          scales['power2'] <- parscale$power2
    
        if(!is.null(parscale$scale))
          scales['scale'] <- parscale$scale

        if(!is.null(parscale$smooth))
          scales['smooth'] <- parscale$smooth
      }

    return(scales[namesparam])
  }

SetRangeParam <- function(namesparam, numparam)
  {
    low <- 1e-12
    lower <- NULL
    upper <- NULL
    
    for(i in 1 : numparam)
      {
        if(namesparam[i]=='df')
          {
            lower <- c(lower, low)
            upper <- c(upper, Inf)
          }

        if(namesparam[i]=='nugget')
          {
            lower <- c(lower, 0)
            upper <- c(upper, 1)
          }

        if(namesparam[i]=='power')
          {
            lower <- c(lower, low)
            upper <- c(upper, 2)
          }

        if(namesparam[i]=='power1')
          {
            lower <- c(lower, low)
            upper <- c(upper, 2)
          }

        if(namesparam[i]=='power2')
          {
            lower <- c(lower, low)
            upper <- c(upper, Inf)
          }

        if(namesparam[i]=='scale')
          {
            lower <- c(lower, low)
            upper <- c(upper, Inf)
          }
        if(namesparam[i]=='smooth')
          {
            lower <- c(lower, low)
            upper <- c(upper, Inf)
          }

      }
    return(list(lower=lower, upper=upper))
  }
