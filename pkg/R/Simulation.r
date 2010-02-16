####################################################
### Authors: Simone Padoan.
### Email: simone.padoan@eofl.ch.
### Institute: EPFL.
### File name: Simulation.r
### Description:
### This file contains a set of procedures
### for the simulation of extremal processes and
### related functions.
### Last change: 21/01/2010.
####################################################


### Procedures are in alphabetical order.

### Simulation of extremal processes:

rExtremal <- function(coord, corrmodel, fitted=NULL, grid=FALSE, model,
                      numblock=300, numsim=1, param)
  {
    result <- NULL

    ### In case the input is passed as fitted model: 
    if(!is.null(fitted))
      {
        if(!inherits(fitted, "extremal"))
          stop("Use only with 'extremal' objects")
        
        corrmodel <- fitted$corrmodel        
        model <- fitted$model
        param <- fitted$param
        
        if(missing(coord))
          coord <- fitted$coord

        if(!is.null(fitted$fixed))
          {
            fixed <- fitted$fixed
            param <- c(fixed, param)
          }
        
        param <- as.list(param)
      }

    ### Check the parameters given in input:
    initsimu <- InitSimu(corrmodel, model, numblock, numsim, param, coord)

    if(!is.null(initsimu$error))
      stop(initsimu$error)
  
    numcoord <- numgrid <- nrow(coord)
    if(grid) numcoord <- numgrid^2


    ### Simulation stage:
    DeleteAllRegisters()
    RFparameters(Storing=TRUE, PrintLevel=1)

    ### Simulation of Extremal-Gaussian random fields: 
    if(model=='Extremal_g')
      {
        coordresc <- RescalingCoord(corrmodel, numblock, initsimu$corrparam, coord)
        
        for(i in 1 : numsim)
          {
            replicates <- double(numcoord)
            ### Simulation Gaussian random field:
            gaussrf <- GaussRF(x = coordresc[, 1], y = coordresc[, 2], model = corrmodel,
                               param = initsimu$corrparam, grid = grid, n = numblock, pch='')
            
            ### Compute compotentwise maxima:
              .C('ComputeMaxima', as.double(initsimu$param['df']), replicates,
                 as.integer(initsimu$model), as.integer(numblock),
                 as.integer(numcoord), gaussrf, PACKAGE='ExtremalProc',
                 DUP = FALSE, NAOK=TRUE)
            
            result <- c(result, replicates)
          }       
      }

    ### Simulation Extremal-t random fields:
    if(model=='Extremal_t')
      {
 
        for(i in 1 : numsim)
          {
            replicates <- double(numcoord)
            # Simulation Gaussian random field (for the Student t process):
            gaussrf <- GaussRF(x = coord[, 1], y = coord[, 2], model = corrmodel,
                            param = initsimu$corrparam, grid = grid, n = numblock, pch='')
            
            # Compute compotentwise maxima:
              .C('ComputeMaxima', as.double(initsimu$param['df']), replicates,
                 as.integer(initsimu$model), as.integer(numblock),
                 as.integer(numcoord), gaussrf, PACKAGE='ExtremalProc',
                 DUP = FALSE, NAOK=TRUE)
            
            result <- c(result, replicates)        
          }
      }

    ### Set the format of the simulated random field:
    if(grid)
      {
        if(numsim > 1)
          dim(result) <- c(numgrid, numgrid, numsim)
        else
          dim(result) <- c(numgrid, numgrid) 
      }
    else
      result <- matrix(result, ncol=numcoord, nrow=numsim, byrow=TRUE)
        
    return(result)
  }
