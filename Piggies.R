# Read CSV into R
#Locate minima to segment sequence (individual cycle)
#turn into functions(FDA) splines (smoothing, lambda penalty)
#time normalize the second deriv using continuous registration <-- acceleration is second deriv, 
  #ideally in rz axis (vertical)

#Where in the cycle does food impact kinematics?
# (tx, ry, rz)*(D,V,A)

#install.packages('pastecs')

#Read CSV into R
setwd('/Users/maraudersmap/Desktop/FDA/Pigs_Data')
library(readr); library(pastecs); library(stringr); library(data.table)
library(dplyr); library(fda)
Pig10_carrot <- read_csv("pig10_20150826_ctrl_carrot_1_TMJBUTTER25_sm250.csv")
View(Pig10_carrot)
Pig10_carrot = as.data.table(Pig10_carrot)
head(Pig10_carrot)

## flag peaks and pits using turnpoints
pig10_tb = turnpoints(Pig10_carrot$TMJdata1.rz)
pig10_peaks = extract(pig10_tb)

## adding in turnpoints
Pig10_carrot[, cycles := pig10_peaks]
Pig10_carrot[TMJdata1.rz > -10, cycles := 0]

pig10_min = Pig10_carrot[cycles == -1]

plot(Pig10_carrot$TMJdata1.rz,type='l')
abline(v = pig10_min$frame, col = 'red')

cycle = 0
Pig10_carrot$Cycle = NA

## Add column for index to order later
index = 1:nrow(Pig10_carrot)
Pig10_carrot$Index = index

## Create new dataframe to put in rows that indicate beginning of cycle
cycle_start = (Pig10_carrot[cycles == -1])[1]
#
for(i in 1:length(Pig10_carrot$TMJdata1.rz)){
  if(Pig10_carrot$cycles[i] == -1){
    cycle = cycle + 1
    ## Add low point to cycle_start
    cycle_start = rbind(cycle_start, Pig10_carrot[i])
  }
  Pig10_carrot$Cycle[i] = cycle
}

## Only select unique cycles
cycle_start = unique(cycle_start)
## Add Cycles to cycles_start
cycle_start$Cycle = (1:nrow(cycle_start)) - 1

## Bind Pig10_carrot with cycle start, thus adding duplicate row
Pig10_carrot = rbind(cycle_start, Pig10_carrot)
## Order by index 
Pig10_carrot = Pig10_carrot[order(Index)]

plot(Pig10_carrot[Cycle == 1]$TMJdata1.rz, type='l', 
     xlim=c(0, 150))
for(i in 2:17){
  lines(Pig10_carrot[Cycle == i]$TMJdata1.rz)
}

## Calculate length of each cycle
cycle_lengths = Pig10_carrot[,.N, by = Cycle]
cycle_lengths

## Normalize so each frame starts at 0 by Cycle
Pig10_carrot = Pig10_carrot[, frame := frame - min(frame), by = Cycle]
## Noramlize soe each frame ends at 1 by Cycle
Pig10_carrot = Pig10_carrot[, frame := frame/max(frame), by = Cycle]

## Function to interprolate data to make all cycles same length
## Use max cycle of 84
inter_func = function(data, cycle){
  max_cyc = 84
  data = as.data.table(data)
  data = data[Cycle == cycle][,.(frame, TMJdata1.rz)]
  n = max_cyc-dim(data)[1]
  interp = as.data.table(approx(data$frame, data$TMJdata1.rz, n = n+2 ))
  colnames(interp) = c('frame', 'TMJdata1.rz')
  data2 = rbind(interp, data)
  data2 = data2[order(frame)][,Cycle := cycle]
  data2 = data2[2:85]
  return(data2)

}

Pig10_Cy1 = inter_func(Pig10_carrot, 1)

Pig10_Cy1

## These are the only cycles we want to use 
cycles = c(1,2,3,4,5,7,8,9,10,11)

## Loop through and add new cycles to a list
Pig10_inter_list = list()
for(i in cycles){
  Pig10_inter_list[[i]] = inter_func(Pig10_carrot, i)
  
}

## Bind these together 
Pig10_inter = rbindlist(Pig10_inter_list)

## Look at plot of raw data
library(ggplot2)
ggplot(Pig10_inter, aes(x = frame, y = TMJdata1.rz, colour = Cycle)) +
  geom_line()+
  facet_wrap(~Cycle)

#turn each cycle into functional data
## Long to Wide
library(tidyr) ## Used for spread function below

## Adding a column to identify position in cylce
id = rep(1:84, 10)
Pig10_inter[,id := id]
## Only select TMJdata1.rz
Pig10_inter_rz = Pig10_inter[,.(id, Cycle, TMJdata1.rz)]

## Go from long to wide
  ## So now, each cycle is a column 
Pig10_wide = spread(Pig10_inter_rz, Cycle, TMJdata1.rz)
## Rename colnames because they were being a pain in the ass
colnames(Pig10_wide) = c('id',paste0('C', 1:5), paste0('C',7:11))
## Don't want id column to do fda
Pig10_wide = Pig10_wide[,.(C1, C2, C3, C4,C5,C7,C8,C9,C10,C11)]

#Create Functional Data
Pig10_fd = function(x){
  x <- as.matrix(x)
  knots=c(0,.2,.4,.6,.8,1.0)
  mybasis <- create.bspline.basis(rangeval=c(0, 1), nbasis=NULL, 
                                  norder=4, breaks=knots, dropind=NULL, 
                                  quadvals=NULL, values=NULL, names="bspl")
  myfd <- Data2fd(x,argvals=seq(0, 1, len = 84 ), basisobj=mybasis)
}

## Call fda on full wide data set
Pigfd_all = Pig10_fd(Pig10_wide)

## Calling fda on just two separate cycles
Pigfd_1 = Pig10_fd(Pig10_inter[Cycle == 1]$TMJdata1.rz)
Pigfd_2 = Pig10_fd(Pig10_inter[Cycle == 2]$TMJdata1.rz)

plot(Pigfd_all)

#get 2nd derivative of each cycle.
#unregistered acceleration function
accelfd = deriv.fd(Pigfd_all, 2)
plot(accelfd)

#unregistered mean function
accelmeanfd= mean(accelfd)
plot(accelmeanfd)

#problems: we need to use the 2nd derivative and apply a 
#  penalty function to it *PEN2(x) = Z[D2x(t)]2dt* 
#  OR a roughness penalty *ROUGH(f)=Z[D2f(t)]2dt*  in order to smooth it like beth did?
# Right now i'm confused. Why does our acceleration look so bad!?
# Coach thinks we need to take the 4th derivative, and add the roughness penalty to it so the
# displacement data will be smoothed appropriately...?? Does that make sense?

## Continuous registration

## First lets rework our functional data object
  ## We will probably want to wrap this bit into a function
  ## So we can do it down the line with different datasets
  ## but thats for future claire and katherine :)

## Our parameters
#smaller lambda --> overfitting 
lambda = 1e-12
# number of splits (see knots above, this is the same as the number of knots we had)
norder = 6
#split 0 to 1 into 84 values (length of our cycles )
samples = seq(0,1,length = 84)
# Define the number of basis functions we want
nbasis = length(samples) + norder -2
# Create basic basis functions with our parameters
mybasis = create.bspline.basis(c(0,1), nbasis, norder, samples)
myfdPar = fdPar(mybasis, 4, lambda)
# Dataset has to be in matrix form
dataset = as.matrix(Pig10_wide)
## Lets make it a fd object based off of our basis functions! woohoo! 
myfd = smooth.basis(samples, dataset, myfdPar)$fd

## Now lets do continuous registration
lambda <- 1
## This should be the same as the nbasis above but are now pulling it from our fd object
nbasis <- myfd$basis$nbasis
ntrials <- dim(Pig10_wide)[2]
#Take the mean of the functional data objects to be our target fd object
y0fd <- mean.fd(myfd)
# functions to be registered to y0fd (the target object)
yfd = myfd
# vectors for our warping function (need to look at this more...)
y0vec = eval.fd(samples, y0fd)
yvec <- eval.fd(samples, yfd)
# Coefficeints for a warping function
coef0 <- matrix(0, nrow = nbasis, ncol = ntrials)
Wfd0 <- fd(coef0, mybasis)
WfdPar <- fdPar(Wfd0, 2, lambda)
#The actual registration we are doing. 
#interlim --> number of iterations
#dbglev --> controls amount of information printed (could be 0, 1, or 2)
reglist <- register.fd(y0fd, yfd, WfdPar, iterlim = 10, dbglev = 1)

#gives the names of things we can call on reglist
names(reglist)


## LETS PLOT EVERYTHING
#plots registered fd objects
plot(reglist$regfd)
#plot original fd obj for comparison
plot(myfd)
#plots the warp functions
plot(reglist$warpfd)


## Calculate means of original fd and registered fd
origMean <- mean.fd(myfd)
regMean <- mean.fd(reglist$regfd)
#plot them on same plot
plot(origMean)
lines(regMean, col = 2)

##########################################################################

#landmark registration Parameters???
fdobj = Pigfd_all
ximarks = #max close, max open, transitions
x0marks = xmeanmarks
ylambda = 1e-10
returnMatrix=TRUE

#landmark registration function
landmarkreg <- function(fdobj, ximarks, x0marks=xmeanmarks,
                        WfdPar=NULL, monwrd=FALSE, ylambda=1e-10,
                        returnMatrix=FALSE){
  fdobj = Pigfd_all
  ximarks = #points we want to use as landmarks on the x-axis (max close)
  ylambda = 1e-10
  
}

MaxClose = landmarkreg(Pigfd_all)
plot(MaxClose)

  #  Arguments:
  #  FDOBJ   ... functional data object for curves to be registered
  #  XIMARKS ... N by NL array of times of interior landmarks for
  #                 each observed curve
  #  XOMARKS ... vector of length NL of times of interior landmarks for
  #                 target curve
  #  WFDPAR  ... a functional parameter object defining a warping function.  
  #                 If NULL, registration is done using linear interpolation
  #                 of lamdmark times in XIMARKS plotted against corresponding 
  #                 target times in X0MARKS.
  #  MONWRD  ... If TRUE, warping functions are estimated by monotone smoothing,
  #                 otherwise by regular smoothing.  The latter is faster, but
  #                 not guaranteed to produce a strictly monotone warping
  #                 function.  If MONWRD is 0 and an error message results
  #                 indicating nonmonotonicity, rerun with MONWRD = 1.
  #                 Default:  TRUE
  #  YLAMBDA ... smoothing parameter to be used in computing the registered
  #                 functions.  For high dimensional bases, local wiggles may be
  #                 found in the registered functions or its derivatives that are
  #                 not seen in the unregistered functions.  In this event, this
  #                 parameter should be increased.
  #  Returns:
  #  FDREG   ... a functional data object for the registered curves
  #  WARPFD  ... a functional data object for the warping functions
  #  WFD     ... a functional data object for the W functions defining the
  #              warping functions
  #  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
  #               from a call to function BsplineS.  See this function for
  #               enabling this option.
  
