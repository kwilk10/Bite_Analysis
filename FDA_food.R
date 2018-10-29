##############################
# Run some FDA by Food Type          #
#
# This script takes the interprolated data from the pigs_inter.R script
# And creates funcational data objects and continuous registration on the data
# It is important to note that continuous registration TAKES A LONG TIME
# It is highly recommended that the user does NOT use this function locally!! 
# 
# Outputs are two R data objects
# --> 1 functional data object, 1 continuously registered functional data object
#
#####################################


library(pastecs); library(stringr); library(data.table)
library(dplyr); library(fda)
#setwd('~/Desktop/FDA')
setwd('~/FDA')


data_wide_func = function(data,food_type, animal_num = NULL){

  if(!is.null(animal_num)){
    data = data[animal_id == animal_num ]
  }
  data = data[food == food_type]
  print(head(data))
  data = data[,.(cyc_id, uniq_id, Cycle, X.TMJdata.rz)]
  data = data[, Cycle_food := .GRP, .(uniq_id)][,.(cyc_id, Cycle_food, X.TMJdata.rz)]
  data_wide = dcast(data, cyc_id ~ Cycle_food, value.var = "X.TMJdata.rz")
  data_wide = data_wide %>% dplyr::select(-cyc_id)
  cols = sprintf('%s_%d', food_type, 1:dim(data_wide)[2] )
  colnames(data_wide) = cols
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



## Continous registration function
### TAKE A LONG TIME: DO NOT RUN LOCALLY ON LARGE DATASETS
cont_reg = function(myfd, data, cycle_length){
  ## Now lets do continuous registration
  lambda <- 1
  norder = 6
  ## This should be the same as the nbasis above but are now pulling it from our fd object
  nbasis <- myfd$basis$nbasis
  ntrials <- dim(data)[2]
  #Take the mean of the functional data objects to be our target fd object
  y0fd <- mean.fd(myfd)
  # functions to be registered to y0fd (the target object)
  yfd = myfd
  # vectors for our warping function (need to look at this more...)
  samples = seq(0,1,length = cycle_length)
  y0vec = eval.fd(samples, y0fd)
  yvec <- eval.fd(samples, yfd)
  # Coefficeints for a warping function
  mybasis = create.bspline.basis(c(0,1), nbasis, norder, samples)
  coef0 <- matrix(0, nrow = nbasis, ncol = ntrials)
  Wfd0 <- fd(coef0, mybasis)
  WfdPar <- fdPar(Wfd0, 2, lambda)
  #The actual registration we are doing. 
  #interlim --> number of iterations
  #dbglev --> controls amount of information printed (could be 0, 1, or 2)
  reglist <- register.fd(y0fd, yfd, WfdPar, iterlim = 10, dbglev = 1)

}


#carrot_cr = cont_reg(carrot_fda, carrot_wide, 163)
#apple_cr = cont_reg(apple_fda, apple_wide, 163)
#almond_cr = cont_reg(almond_fda, almond_wide, 163)

## Save our continuously registered fda
#saveRDS(carrot_cr, "/by_food/carrot_cr.rds")
#saveRDS(apple_cr, "/by_food/apple_cr.rds")
#saveRDS(almond_cr, "/by_food/almond_cr.rds")


## This function outputs our regression plots with our error confidence bounds
## Arguments include: a functional data object, either the one from above or a continuous registered one
  ## for continuous registered you must input argument as fda_data$regfd
  ## full data is the data that contains all the foods into one large, wide data set
  ## food is the food that you want to plot. I can't get them to work as one, so I'm doing them separately
  ## stderrList is optional. If you don't put it in, then it will calculate it for you. This is what takes forever
    # hence, I have been running this separately and saving it so I can just reuse it to see the plots. 
    # if you don't have it, this function will save it for you.
    # You do have to have a folder called by_food in your working directory
  ## animal number is if you are doing the regression on a single animal. Its optional as well

plt_func = function(fda_data, full_data, type, food, stderrList = NULL, animal_num = NULL){
  
  carrot_wide = full_data %>% dplyr:: select(starts_with('carrot'))
  almond_wide = full_data %>% dplyr:: select(starts_with('almond'))
  apple_wide = full_data %>% dplyr:: select(starts_with('apple'))
 
  
  s_basis <- fda_data$basis
  
  basismat <- eval.basis(seq(0,1,length=163),s_basis)
  ## this maps the observations to the coefficients (??)
  y2cMap <- solve(crossprod(basismat)) %*% t(basismat)
  
  C_reps = dim(carrot_wide)[2]
  Al_reps = dim(almond_wide)[2]
  Ap_reps = dim(apple_wide)[2]
  
  reps = dim(full_data)[2]
  
  ## Groups names --> Names of animals 
  groupnames = c('Carrot','Almond','Apple')
  # indices for groups

  Cindex <- seq(1,C_reps)
  Alindex <- seq(C_reps+1,C_reps + Al_reps)
  Apindex <- seq(Alindex[Al_reps]+1, Alindex[Al_reps] + Ap_reps)
  
  gmat <- matrix(0, reps, 2)
  # Fist column = 1s
  gmat[ ,1] <- 1
  # Second column = groups
  if(food == 'carrot'){
    gmat[Cindex, 2] <- 1
  }
  
  if(food == 'almond'){ 
    gmat[Alindex,2] <- 1
  }
  
  if(food == 'apple'){
    gmat[Apindex, 2] <-1
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
    file = sprintf('by_food/stderrList_%s_%s', type, food)
    saveRDS(stderrList, file)
  }
  
  
  betastderrlist <- stderrList$betastderrlist
  

  #plotbeta(fRegressList$betaestlist, betastderrlist)
  xlabel = sprintf('%s_%s',type, food )
  plotbeta(fRegressList$betaestlist, betastderrlist, xlab = xlabel)
  
  

  
  
}

####################################################################################################
# Okay, run everything up to this point. So far we have just created a bunch of functions. ##
# Now we run everything on our data and see what we get! 
####################################################################################################



## Read our interprolated data in
pigs = fread('pig_cycles_int.csv')

## Add a column indicating place in cycle
## Cycle length = 163
## # of cycles = 1122
id = rep(1:163, 1122)
pigs[,cyc_id := id]


## First off, we are going to create our wide data and bind it all together into one master data.table
carrot_wide = data_wide_func(pigs, 'carrot')
almond_wide = data_wide_func(pigs, 'almond')
apple_wide = data_wide_func(pigs, 'apple')
full_data = cbind(carrot_wide, almond_wide, apple_wide)

## We can also do similar for a single animal
carrot10 = data_wide_func(pigs, 'carrot', 10)
almond10 = data_wide_func(pigs, 'almond', 10)
apple10 = data_wide_func(pigs, 'apple', 10)
full10 = cbind(carrot10, almond10, apple10)

## Okay and now we would create our functional data object on our FULL data set
## We could also do this separately for all the foods and animals and compare means that way.
## We need it on the full data set to run our regression down the line
full_fda = fda_func(full_data, 20)
fda_carrot = fda_func(carrot_wide, 20)
fda_apple = fda_func(apple_wide, 20)
fda_almond = fda_func(almond_wide, 20)
#This is also a step we can plot at
plot(fda_full)
plot(mean(fda_full))

## If you highlight all of these and select run, they should all show up on the same plot
plot(mean(fda_carrot), col = 'green')
par(new = TRUE)
lines(mean(fda_apple), col = 'red')
lines(mean(fda_almond), col = 'purple')


## So, I have not included continuous registration function in here because it takes forever to run
## And I have been running it on the umich servers overnight. 
## But, I think I sent you some and will try to either upload or send you what I am using below

# Our final step at this point is to plot the regression. That is our last function from above
## Here is the continuous registered object have used
## it has 20 splines
full_cr = readRDS('by_food/full_cr_20.rds')

## I have also already run a ton to get the standard deviation matricies because that part takes forever to run
cr_c = readRDS('by_food/stderr_obj/stderrList_continuous_reg_carrot')
ap_c = readRDS('by_food/stderr_obj/stderrList_continuous_reg_apple')
al_c = readRDS('by_food/stderr_obj/stderrList_continuous_reg_almond')

cr_fd = readRDS('by_food/stderr_obj/stderrList_full_fda_carrot')
ap_fd = readRDS('by_food/stderr_obj/stderrList_full_fda_apple')
al_fd = readRDS('by_food/stderr_obj/stderrList_full_fda_almond')

## Okay, now that we have those, super simple! 
## Just run the following and input what functional data object and food that you want to plot! 
plt_func(fda_data = full_cr$regfd, full_data, type = 'continuous_reg', food = 'carrot', stderrList = cr_c)
plt_func(full_cr$regfd, full_data,'continuous_reg','apple', ap_c)
plt_func(full_cr$regfd, full_data,'continuous_reg','almond', al_c)


plt_func(full_fda, full_data,'full_fda','carrot', cr_fd)
plt_func(full_fda, full_data,'full_fda','apple', ap_fd)
plt_func(full_fda, full_data,'full_fda','almond', al_fd)
