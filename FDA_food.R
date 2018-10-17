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


library(readr); library(pastecs); library(stringr); library(data.table)
library(dplyr); library(fda)
#setwd('~/Desktop/FDA')
setwd('~/FDA')

## Read our interprolated data in
pigs = fread('pig_cycles_int.csv')


## Add a column indicating place in cycle
  ## Cycle length = 163
  ## # of cycles = 1122
id = rep(1:163, 1122)
pigs[,cyc_id := id]

## subset by food
carrot = pigs[food == 'carrot']

carrot = carrot[,.(cyc_id, Cycle, X.TMJdata.rz)]
carrot = carrot[, Cycle_carrot := .GRP, .(uniq_id)][,.(cyc_id, Cycle_carrot, X.TMJdata.rz)]
# Put into wide data set so that each cycle is a column
carrot_wide = dcast(carrot, cyc_id ~ Cycle_carrot, value.var = "X.TMJdata.rz")
carrot_wide = carrot_wide %>% dplyr::select(-cyc_id)


apple = pigs[food == 'apple']
apple = apple[,.(cyc_id, uniq_id, Cycle, X.TMJdata.rz)]
apple = apple[, Cycle_apple := .GRP, .(uniq_id)][,.(cyc_id, Cycle_apple, X.TMJdata.rz)]
apple_wide = dcast(apple, cyc_id ~ Cycle_apple, value.var = "X.TMJdata.rz")
apple_wide = apple_wide %>% dplyr::select(-cyc_id)


almond = pigs[food == 'almond']
almond = almond[,.(cyc_id, uniq_id, Cycle, X.TMJdata.rz)]
almond = almond[, Cycle_almond := .GRP, .(uniq_id)][,.(cyc_id, Cycle_almond, X.TMJdata.rz)]
almond_wide = dcast(almond, cyc_id ~ Cycle_almond, value.var = "X.TMJdata.rz")
almond_wide = almond_wide %>% dplyr::select(-cyc_id)


#Create Functional Data
# This creates a very basic functional data object
fda_func = function(x, cycle_length){
  x <- as.matrix(x)
  knots=c(0,.2,.4,.6,.8,1.0)
  mybasis <- create.bspline.basis(rangeval=c(0, 1), nbasis=NULL, 
                                  norder=4, breaks=knots, dropind=NULL, 
                                  quadvals=NULL, values=NULL, names="bspl")
  myfd <- Data2fd(x,argvals=seq(0, 1, len = cycle_length ), basisobj=mybasis)
}


carrot_fda = fda_func(carrot_wide, 163)
apple_fda = fda_func(apple_wide, 163)
almond_fda = fda_func(almond_wide, 163)
# We can plot them for funsies
plot(apple_fda)
plot(almond_fda)



## This creates an oversaturated functional data object 
## With the same number of splines as there are data points
##  In this case, this means one spline for each point in the cycle (163 points)
fda_func2 = function(x, cycle_length){
  ## Our parameters
  #smaller lambda --> overfitting 
  lambda = 1e-12
  # number of splits (see knots above, this is the same as the number of knots we had)
  norder = 6
  #split 0 to 1 into 84 values (length of our cycles )
  samples = seq(0,1,length = cycle_length)
  # Define the number of basis functions we want
  nbasis = length(samples) + norder -2
  # Create basic basis functions with our parameters
  mybasis = create.bspline.basis(c(0,1), nbasis, norder, samples)
  myfdPar = fdPar(mybasis, 4, lambda)
  # Dataset has to be in matrix form
  dataset = as.matrix(x)
  ## Lets make it a fd object based off of our basis functions! woohoo! 
  myfd = smooth.basis(samples, dataset, myfdPar)$fd

}

carrot_fda = fda_func2(carrot_wide, 163)
almond_fda = fda_func2(almond_wide, 163)
apple_fda = fda_func2(apple_wide, 163)

## Saving our functional data objects
saveRDS(carrot_fda, "/by_food/carrot_fda.rds")
saveRDS(apple_fda, "/by_food/apple_fda.rds")
saveRDS(almond_fda, "/by_food/almond_fda.rds")


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


carrot_cr = cont_reg(carrot_fda, carrot_wide, 163)
apple_cr = cont_reg(apple_fda, apple_wide, 163)
almond_cr = cont_reg(almond_fda, almond_wide, 163)

#gives the names of things we can call on reglist
names(carrot_cr)

## Save our continuously registered fda
saveRDS(carrot_cr, "/by_food/carrot_cr.rds")
saveRDS(apple_cr, "/by_food/apple_cr.rds")
saveRDS(almond_cr, "/by_food/almond_cr.rds")


## After they are saved we can read them all in again and not have to run any of the above
carrot_cr = readRDS('by_food/carrot_cr.rds')
carrot_fda = readRDS('by_food/carrot_fda.rds')

apple_cr = readRDS('by_food/apple_cr.rds')
apple_fda = readRDS('by_food/apple_fda.rds')

almond_cr = readRDS('by_food/almond_cr.rds')
almond_fda = readRDS('by_food/almond_fda.rds')

## LETS PLOT EVERYTHING
#plots registered fd objects
plot(mean(carrot_cr$regfd))
#plot original fd obj for comparison
plot(mean(carrot_fda))
#plots the warp functions
plot(almond_cr$warpfd)

FirstDeriv = deriv.fd(apple_cr$regfd, 1)
plot(mean(FirstDeriv))
