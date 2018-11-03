##############################
# Run some FDA and regression by Animal number         #
#
# This script takes the interprolated data from the pigs_inter.R script
# And creates funcational data objects on the data and plots the regression
# curves with a confidence bound showing how each of the animals
# differ from the overall mean. 
# 
# Data needed: the pigs_cycles_int.csv
#   If previously run, it is recommended that the user use
#   the saved stderrList objects
#
#####################################

library(stringr); library(data.table)
library(dplyr); library(fda)
#setwd('~/Desktop/FDA')
setwd('~/Desktop/FDA')


data_wide_func = function(data,food_type = NULL, animal_num){
  
  #if(!is.null(animal_num)){
  #  data = data[animal_id == animal_num ]
  #}
  data = data[animal_id == animal_num]
  print(head(data))
  data = data[,.(cyc_id, uniq_id, Cycle, X.TMJdata.rz)]
  data = data[, Cycle_food := .GRP, .(uniq_id)][,.(cyc_id, Cycle_food, X.TMJdata.rz)]
  data_wide = dcast(data, cyc_id ~ Cycle_food, value.var = "X.TMJdata.rz")
  data_wide = data_wide %>% dplyr::select(-cyc_id)
  cols = sprintf('%s_%d', animal_num, 1:dim(data_wide)[2] )
  colnames(data_wide) = cols
  return(data_wide)
}


fda_func = function(data, splines){
  ## Our parameters
  #smaller lambda --> overfitting 
  lambda = 1e-12
  # number of splits (see knots above, this is the same as the number of knots we had)
  norder = 6
  #split 0 to 1 into 84 values (length of our cycles )
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

plt_func = function(fda_data, full_data, animal_num, type, food = NULL, stderrList = NULL){
  
  carrot_wide = full_data %>% dplyr:: select(starts_with('carrot'))
  almond_wide = full_data %>% dplyr:: select(starts_with('almond'))
  apple_wide = full_data %>% dplyr:: select(starts_with('apple'))
  
  pig5_wide = full_data %>% dplyr:: select(starts_with('5'))
  pig6_wide = full_data %>% dplyr:: select(starts_with('6'))
  
  pig9_wide = full_data %>% dplyr:: select(starts_with('9'))
  pig10_wide = full_data %>% dplyr:: select(starts_with('10'))
  print(head(pig10_wide))
  
  
  
  
  
  
  s_basis <- fda_data$basis
  
  basismat <- eval.basis(seq(0,1,length=163),s_basis)
  ## this maps the observations to the coefficients (??)
  y2cMap <- solve(crossprod(basismat)) %*% t(basismat)
  
  pig5_reps = dim(pig5_wide)[2]
  pig6_reps = dim(pig6_wide)[2]
  pig9_reps = dim(pig9_wide)[2]
  pig10_reps = dim(pig10_wide)[2]
  
  
  reps = dim(full_data)[2]
  
  ## Groups names --> Names of animals 
  groupnames = c('Pig5','Pig6','Pig9','Pig10')
  # indices for groups
  
  fiveindex <- seq(1,pig5_reps)
  sixindex <- seq(pig5_reps+1,pig5_reps + pig6_reps)
  nineindex <- seq(sixindex[pig6_reps]+1, sixindex[pig6_reps] + pig9_reps)
  tenindex <- seq(nineindex[pig9_reps]+1, nineindex[pig9_reps] + pig10_reps)
  
  gmat <- matrix(0, reps, 2)
  # Fist column = 1s
  gmat[ ,1] <- 1
  # Second column = groups
  if(animal_num == 5){
    gmat[fiveindex, 2] <- 1
  }
  
  if(animal_num == 6){ 
    gmat[sixindex,2] <- 1
  }
  
  if(animal_num == 9){
    gmat[nineindex, 2] <-1
  }
  
  if(animal_num == 10){
    gmat[tenindex, 2] <-1
  }
  
  
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
  #add an xlab = xlabel to the following to add label to the plots
  # can adjust to anything you would like
  plotbeta(fRegressList$betaestlist, betastderrlist)
  
  
  
  
  
}


## Read our interprolated data in
pigs = fread('pig_cycles_int.csv')

## Add a column indicating place in cycle
## Cycle length = 163
## # of cycles = 1122
id = rep(1:163, 1122)
pigs[,cyc_id := id]


data_wide_func(pigs, animal_num = 5, )
pig6 = data_wide_func(pigs, animal_num = 6)
pig9 = data_wide_func(pigs, animal_num = 9)
pig10 = data_wide_func(pigs, animal_num = 10)

full_data = cbind(pig5, pig6, pig9, pig10)

animal_fda = fda_func(full_data, 20)

## These are the previously saved standard deviation list objects
  # These are needed to create the confidence bounds for regression
five_fd = readRDS('by_animal/stderrList_fda_5')
six_fd = readRDS('by_animal/stderrList_fda_6')
nine_fd = readRDS('by_animal/stderrList_fda_9')
ten_fd = readRDS('by_animal/stderrList_fda_10')

# This outputs the plots
plt_func(animal_fda, full_data, animal_num = 5, 'fda', food = NULL, stderrList = five_fd)
plt_func(animal_fda, full_data, animal_num = 6, 'fda', food = NULL, stderrList = six_fd)
plt_func(animal_fda, full_data, animal_num = 9, 'fda', food = NULL, stderrList = nine_fd)
plt_func(animal_fda, full_data, animal_num = 10, 'fda', food = NULL, stderrList = ten_fd)
