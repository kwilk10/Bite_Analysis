## Common Functions that are used in FDA analysis
## Includes:
## data_wide_func --> turns data long data frame into wide data frame
##                    Each column become a single cycle
## fda_func --> turns raw cycle data into functional data objects
## plot_beta_n --> taken from fda package plot function in github 
##                and adjusted for our plotting needs
## plt_func_a --> runs an fAnova to compare across animals
## plt_func_f --> runs an fAnova to compare across food
## plt_func_int --> runs an fAnova to compare interactions (food and animal)
##
##
## Updated 03/28/19 -- KSW
#############################################################################

library(readr); library(pastecs); library(stringr); library(data.table)
library(dplyr); library(fda)

setwd('~/Desktop/FDA')

normalized = function(x) ( (x- max(x))/(max(x)-min(x) ) )

data_wide_func = function(data,food_type, axis, animal_num = NULL){
  
  if(axis == "X.TMJdata.rz"){
    ## First normalize the axis data so that all cycles peak at 1
    #data_s = data[,lapply(.SD, normalized), .SDcols = axis, by = uniq_id ]
    #bd = data[,.(ref_frame, Cycle, animal, animal_id, food, date, ID, lab)]
  
    #data = cbind(data_s, bd)
  }
  
  if(!is.null(animal_num)){
    data = data[animal_id == animal_num ]
  }
  data = data[food == food_type]
  
  
  data = data %>% select('lab','uniq_id','Cycle', axis)
  
  #data = data[, Cycle_food := .GRP, .(uniq_id)][,.(cyc_id, Cycle_food, X.TMJdata.rz)]
  data = data[, Cycle_food := .GRP, .(uniq_id)]
  data = data %>% select('lab', 'Cycle_food', axis)
  
  #data_wide = dcast(data, cyc_id ~ Cycle_food, value.var = "X.TMJdata.rz")
  #data_wide = data_wide %>% dplyr::select(-cyc_id)
  
  data_wide = dcast(data, lab ~ Cycle_food, value.var = axis)
  data_wide = data_wide %>% dplyr::select(-lab)
  
  if(!is.null(animal_num)){
    cols = sprintf('%s_%s_%d', animal_num, food_type, 1:dim(data_wide)[2] )
    colnames(data_wide) = cols
  }
  else{
    cols = sprintf('%s_%d', food_type, 1:dim(data_wide)[2] )
    colnames(data_wide) = cols
  }
  return(data_wide)
}

## This outputs a regular functional data object
## The arguments are your data, in wide format (so 1 column for each cycle)
# splines --> number of splines you want to use to fit your data (I think...)
# So if your splines = cycle length, then you are overfitting 
# and have a spline for each data point. 
# The smaller our splines, the less time the next function takes to run
# we can play around with this. 

fda_func = function(data, splines){
  ## Our parameters
  #smaller lambda --> overfitting 
  lambda = 1e-12
  # number of splits (see knots above, this is the same as the number of knots we had)
  norder = 6
  #split 0 to 1 into number of splines we choose
  samples = seq(0,1,length = splines)
  #samples = seq(0,1,length = cycle_length)
  
  # Define the number of basis functions we want
  nbasis = length(samples) + norder -2
  # Create basic basis functions with our parameters
  mybasis = create.bspline.basis(c(0,1), nbasis, norder, samples)
  myfdPar = fdPar(mybasis, 4, lambda)
  # Dataset has to be in matrix form
  dataset = as.matrix(data)
  ## Lets make it a fd object based off of our basis functions! woohoo! 
  myfd = smooth.basis(seq(0,1,length = 163), dataset, myfdPar)$fd
  
}

## This is a function to plot the regression outputs for FDA
# ylim = NULL --> y limits of plot are preset to fit the data
# In order to set y limits to desired values, ylim = c(min, max)
# plot_type = 1 --> outputs overall mean
# plot_type = 2 --> output regression comparison 

plotbeta_n = function(betaestlist, betastderrlist, argvals=NULL,
                    xlab="", plot_type = 1, ylim = NULL, ...)
{
  #  PLOTBETA plots a functional parameter along with confidence
  #  limits
  #  Arguments
  #  BETAESTLIST    ... A list object containing one or more functional
  #                     parameter objects or functional data objects.
  #  BETASTDERRLIST ... A list object containing functional data objects
  #                     for the standard error of the objects in
  #                     BETAESTLIST.
  
  #  Last modified 12 December 2008
  
  #  check BETAESTLIST
  
  if (inherits(betaestlist, "fdPar") || inherits(betaestlist, "fd")) {
    betaestlist = list(betaestlist)
  }
  
  if (!inherits(betaestlist, "list")) {
    stop("BETAESTLIST is not a list, fd, or fdpar object.")
  }
  
  #  check BETASTDERRLIST
  
  if (inherits(betastderrlist, "fd")) {
    betastderrlist = list(betastderrlist)
  }
  
  if (!inherits(betastderrlist, "list")) {
    stop("BETASTDERRLIST is not a list, or fd object.")
  }
  
  #  get range
  
  if (is.fdPar(betaestlist[[1]])) {
    rangeval = betaestlist[[1]]$fd$basis$rangeval
  } else {
    if (is.fd(betaestlist[[1]])) {
      rangeval = betaestlist[[1]]$basis$rangeval
    } else {
      stop(paste("A list does not contain either a functional parameter ",
                 "or a functional data object."))
    }
  }
  
  if (is.null(argvals)) {
    argvals = seq(rangeval[1],rangeval[2],len=51)
  }
  n = length(argvals)
  p = length(betaestlist)
  
  par(ask=F)
  #for (j in 1:p) {
  
  if(plot_type == 1){
    t = 1}
  if(plot_type == 2){
    t = p
  }
  
  for(j in t){
    if (is.fdPar(betaestlist[[j]])) {
      betavec = eval.fd(argvals, betaestlist[[j]]$fd)
    } else {
      if (is.fd(betaestlist[[j]])) {
        betavec = eval.fd(argvals, betaestlist[[j]])
      } else {
        stop(
          "BETAESTLIST does not contain a functional parameter or data object.")
      }
    }
    betastderr = eval.fd(argvals, betastderrlist[[j]])
    betavecp   = betavec + 2*betastderr
    betavecm   = betavec - 2*betastderr
    zeroval  = c(0,0)
    if(is.null(ylim) == FALSE){
      plot(argvals, betavec, type="l", xlab=xlab, ylab="",
           xlim=rangeval, 
           ylim = ylim, ...)
      
    }
    
    else{
      plot(argvals, betavec, type="l", xlab=xlab, ylab="",
           xlim=rangeval, 
           ylim=c(min(betavecm),max(betavecp)), ...)
      
    }
    lines(rangeval, zeroval,lty=3, col=2)
    lines(argvals, betavec,  col=1, lwd=2)
    lines(argvals, betavecp, col=1, lwd=1)
    lines(argvals, betavecm, col=1, lwd=1)
    lines(rangeval, zeroval, col=1, lty=3)
    for (i in 1:n) {
      lines(c(argvals[i],argvals[i]),c(betavecm[i],betavecp[i]))
    }
    #title(paste("Regression function ",j))
  }
  
}

####################################################################################

# This function takes fda data and raw data and runs an fAnova to compare across animals

####################################################################################


plt_func_a = function(fda_data, full_data, animal_num, type, 
                      stderrList = NULL,
                      ylim = NULL, plot_type = 1){
  
  s_basis <- fda_data$basis
  
  basismat <- eval.basis(seq(0,1,length=163),s_basis)
  ## this maps the observations to the coefficients (??)
  y2cMap <- solve(crossprod(basismat)) %*% t(basismat)
  
  
  reps = dim(full_data)[2]
  
  indx = which(colnames(full_data) %like% sprintf('%i_', animal_num)) 
  
  print(indx)
  
  gmat <- matrix(0, reps, 2)
  # Fist column = 1s
  gmat[ ,1] <- 1
  
  gmat[indx,2] <- 1
  ## Set up xfdlist for regression
  #this is for our covariates
  #First covariate is just the intercept so a bunch of 1's
  #Next covariate is carrot, 3rd if almond, 4th is apple
  p <- 2
  xfdlist <- vector("list",p)
  for (j in 1:p) xfdlist[[j]] <- gmat[,j]
  
  
  bwtlist = list(fdPar(s_basis,2,0), fdPar(s_basis,2,0))
  fRegressList <- fRegress(fda_data, xfdlist, bwtlist)
  print(names(fRegressList))
  
  # get predicted functions
  yhatfdobj <- fRegressList$yhatfdobj
  
  # compute residual matrix and get covariance of residuals
  cycletime <- seq(0,1, length=163)
  #predicted from regression
  yhatmat <- eval.fd(cycletime, yhatfdobj$fd)
  #from the actual data
  ymat <- eval.fd(cycletime, fda_data)
  #errors 
  speciesmat <- ymat[,1:reps] - yhatmat[,1:reps]
  
  SigmaE <- var(t(speciesmat))
  
  if(is.null(stderrList)){
    print('running stderrList')
    stderrList <- fRegress.stderr(fRegressList, y2cMap, SigmaE)
    print('finished running stderrList!')
    file = sprintf('by_animal/stderrList_%s_%s', type, animal_num)
    saveRDS(stderrList, file)
  }
  
  betastderrlist <- stderrList$betastderrlist
  
  #plotbeta(fRegressList$betaestlist, betastderrlist)
  xlabel = sprintf('%s_%s',type, animal_num )
  #plot.new()
  #plot.window(xlim=c(0,1), ylim=c(-6,6))
  #plotbeta_n(fRegressList$betaestlist, betastderrlist, plot_type = plot_type, ylim = ylim)
  
  return(list(fRegressList$betaestlis, betastderrlist))
  
}

####################################################################################

# This function takes fda data and raw data and runs an fAnova to compare across food

####################################################################################


plt_func_f = function(fda_data, full_data, food, type, 
                      stderrList = NULL,
                      ylim = NULL, plot_type = 1){
  
  s_basis <- fda_data$basis
  
  basismat <- eval.basis(seq(0,1,length=163),s_basis)
  ## this maps the observations to the coefficients (??)
  y2cMap <- solve(crossprod(basismat)) %*% t(basismat)
  
  
  reps = dim(full_data)[2]
  
  indx = which(colnames(full_data) %like% sprintf('%s_', food)) 

  
  print(indx)
  
  gmat <- matrix(0, reps, 2)
  # Fist column = 1s
  gmat[ ,1] <- 1
  
  gmat[indx,2] <- 1
  ## Set up xfdlist for regression
  #this is for our covariates
  #First covariate is just the intercept so a bunch of 1's
  #Next covariate is carrot, 3rd if almond, 4th is apple
  p <- 2
  xfdlist <- vector("list",p)
  for (j in 1:p) xfdlist[[j]] <- gmat[,j]
  
  
  bwtlist = list(fdPar(s_basis,2,0), fdPar(s_basis,2,0))
  fRegressList <- fRegress(fda_data, xfdlist, bwtlist)
  print(names(fRegressList))
  
  # get predicted functions
  yhatfdobj <- fRegressList$yhatfdobj
  
  # compute residual matrix and get covariance of residuals
  cycletime <- seq(0,1, length=163)
  #predicted from regression
  yhatmat <- eval.fd(cycletime, yhatfdobj$fd)
  #from the actual data
  ymat <- eval.fd(cycletime, fda_data)
  #errors 
  speciesmat <- ymat[,1:reps] - yhatmat[,1:reps]
  
  SigmaE <- var(t(speciesmat))
  
  if(is.null(stderrList)){
    print('running stderrList')
    stderrList <- fRegress.stderr(fRegressList, y2cMap, SigmaE)
    print('finished running stderrList!')
    file = sprintf('by_animal/stderrList_%s_%s', type, food)
    saveRDS(stderrList, file)
  }
  
  betastderrlist <- stderrList$betastderrlist
  
  #plotbeta(fRegressList$betaestlist, betastderrlist)
  xlabel = sprintf('%s_%s',type, food )
  #plot.new()
  #plot.window(xlim=c(0,1), ylim=c(-6,6))
  #plotbeta_n(fRegressList$betaestlist, betastderrlist, plot_type = plot_type, ylim = ylim)
  
  return(list(fRegressList$betaestlis, betastderrlist))
  
}
 
####################################################################################

# This function takes fda data and raw data and runs an fAnova to compare interactions (food and animal)

####################################################################################

plt_func_int = function(fda_data, full_data, animal_num = NULL, food = NULL, type, 
                   stderrList = NULL){
  
  s_basis <- fda_data$basis
  
  basismat <- eval.basis(seq(0,1,length=163),s_basis)
  ## this maps the observations to the coefficients (??)
  y2cMap <- solve(crossprod(basismat)) %*% t(basismat)
  
  
  reps = dim(full_data)[2]
  
  indx = which(colnames(full_data) %like% sprintf('%s_%i', food, animal_num)) 
  
  gmat <- matrix(0, reps, 2)
  # Fist column = 1s
  gmat[ ,1] <- 1
  
  gmat[indx,2] <- 1
  
  ## Set up xfdlist for regression
  #this is for our covariates
  #First covariate is just the intercept so a bunch of 1's
  #Next covariate is carrot, 3rd if almond, 4th is apple
  p <- 2
  xfdlist <- vector("list",p)
  for (j in 1:p) xfdlist[[j]] <- gmat[,j]
  
  
  bwtlist = list(fdPar(s_basis,2,0), fdPar(s_basis,2,0))
  fRegressList <- fRegress(fda_data, xfdlist, bwtlist)
  print(names(fRegressList))
  
  # get predicted functions
  yhatfdobj <- fRegressList$yhatfdobj
  
  # compute residual matrix and get covariance of residuals
  cycletime <- seq(0,1, length=163)
  #predicted from regression
  yhatmat <- eval.fd(cycletime, yhatfdobj$fd)
  #from the actual data
  ymat <- eval.fd(cycletime, fda_data)
  #errors 
  speciesmat <- ymat[,1:reps] - yhatmat[,1:reps]
  
  SigmaE <- var(t(speciesmat))
  
  if(is.null(stderrList)){
    print('running stderrList')
    stderrList <- fRegress.stderr(fRegressList, y2cMap, SigmaE)
    print('finished running stderrList!')
    file = sprintf('by_food/stderrList_%s_%s%s', type, food, animal_num)
    saveRDS(stderrList, file)
  }
  
  
  betastderrlist <- stderrList$betastderrlist
  
  
  #plotbeta(fRegressList$betaestlist, betastderrlist)
  #xlabel = sprintf('%s_%s',type, animal_num )
  #plotbeta_n(fRegressList$betaestlist, betastderrlist, plot_type = plot_type, ylim = ylim)
  
  
  return(list(fRegressList$betaestlis, betastderrlist))
  
  
  
  
}

