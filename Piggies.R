# Read CSV into R
#Locate minima to segment sequence (individual cycle)
#turn into functions(FDA) splines (smoothing, lambda penalty)
#time normalize the second deriv using continuous registration <-- acceleration is second deriv, 
  #ideally in rz axis (vertical)

#Where in the cycle does food impact kinematics?
# (rx, ry, rz)*(D,V,A)

## This is an example change! Hiiiiiiii


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


##
##
cycle = 0
Pig10_carrot$Cycle = NA

########################################################################
## This is the start of the changes
#####################################################################

## Add column for index to order later
index = 1:nrow(Pig10_carrot)
Pig10_carrot$Index = index

## Create new dataframe to put in rows that indicate beginning of cycle
cycle_start = (Pig10_carrot[cycles == -1])[1]

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



## Calculate lenght of each cycle
cycle_lengths = Pig10_carrot[,.N, by = Cycle]
cycle_lengths

## Normalize so each frame starts at 0 by Cycle
Pig10_carrot = Pig10_carrot[, frame := frame - min(frame), by = Cycle]
## Noramlize soe each frame ends at 1 by Cycle
Pig10_carrot = Pig10_carrot[, frame := frame/max(frame), by = Cycle]

## We really don't need any of this...
Pig10_carrot[Cycle == 1]
length(Pig10_carrot[Cycle == 1]$TMJdata1.rz)
c1_int = as.data.table(approx(Pig10_c1$frame, Pig10_c1$TMJdata1.rz, n = 21 ))
colnames(c1_int) = c('frame', 'TMJdata1.rz')
c1_new = rbind(c1_int, Pig10_c1)
c1_new = c1_new[order(frame)]
c1_new = c1_new[2:84]

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


######
#time normalize each chew cycle to have the same X-axis
# Plot the curves starting at time zero 
times <- Pig10_carrot$frame[Pig10_carrot$Cycle == 1]
times.starting.at.0 <- times - min(times) 
plot(times.starting.at.0, 
     Pig10_carrot[Cycle == 1]$TMJdata1.rz, type='l', 
     xlim=c(0, 140))
for(i in 2:17){
  times <- Pig10_carrot$frame[Pig10_carrot$Cycle == i]
  times.starting.at.0 <- times - min(times) 
  lines(times.starting.at.0, Pig10_carrot[Cycle == i]$TMJdata1.rz)
}

# Plot curves with normalized time 
times <- Pig10_carrot$frame[Pig10_carrot$Cycle == 1]
times.starting.at.0 <- times - min(times) 
normalized.time <- times.starting.at.0 / max(times.starting.at.0)
plot(normalized.time, 
     Pig10_carrot[Cycle == 1]$TMJdata1.rz, type='l')
for(i in 2:17){
  times <- Pig10_carrot$frame[Pig10_carrot$Cycle == i]
  times.starting.at.0 <- times - min(times) 
  normalized.time <- times.starting.at.0 / max(times.starting.at.0)
  lines(normalized.time, Pig10_carrot[Cycle == i]$TMJdata1.rz)
}

plt <- as.data.frame(Pig10_carrot)
for (j in 2:7) {
  times <- plt$frame[plt$Cycle == 1]
  times.starting.at.0 <- times - min(times) 
  normalized.time <- times.starting.at.0 / max(times.starting.at.0)
  pdf.filename <- paste0(
    "FDA/Pig10_carrot_", str_replace(string = names(plt)[j], pattern = "[.]", replacement = "_"), ".pdf")
  pdf(pdf.filename)
  plot(normalized.time,  
       plt[plt$Cycle == 1, j], type='l', ylim=range(plt[, j]), 
       main=names(plt)[j], 
       xlab="Normalized time", ylab="Degrees")
  for(i in 2:17){
    times <- plt$frame[plt$Cycle == i]
    times.starting.at.0 <- times - min(times) 
    normalized.time <- times.starting.at.0 / max(times.starting.at.0)
    lines(normalized.time, plt[plt$Cycle == i, j])
  }
  dev.off()
}

#########
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


Pig10_fd = function(x){
  x <- as.matrix(x)
  print('hello its working')
  
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

#get 1st and 2nd derivatives of each cycle.
#unregistered acceleration function
accelfd = deriv.fd(Pigfd_all, 2)
plot(accelfd)

#unregistered mean function
accelmeanfd= mean(accelfd)
plot(accelmeanfd)

#decide landmark registration or continuous registration
#landmark registration


lambda = 1e-12
norder = 6
samples = seq(0,1,length = 83)
nbasis = length(samples) + norder -2
mybasis = create.bspline.basis(c(0,1), nbasis, norder, samples)
myfdPar = fdPar(mybasis, 4, lambda)
myfd_c1 = smooth.basis(samples, Pig10_wide$C1, myfdPar)$fd
plot(myfd_c1)
myfd_c2 = smooth.basis(samples, Pig10_wide$C2, myfdPar)$fd
lines(myfd_c2)
myfd_c3 = smooth.basis(samples, Pig10_wide$C3, myfdPar)$fd
lines(myfd_c3)
myfd_c4 = smooth.basis(samples, Pig10_wide$C4, myfdPar)$fd
lines(myfd_c4)
myfd_c5 = smooth.basis(samples, Pig10_wide$C5, myfdPar)$fd
lines(myfd_c5)
myfd_c7 = smooth.basis(samples, Pig10_wide$C7, myfdPar)$fd
lines(myfd_c7)
myfd_c8 = smooth.basis(samples, Pig10_wide$C8, myfdPar)$fd
lines(myfd_c8)
myfd_c9 = smooth.basis(samples, Pig10_wide$C9, myfdPar)$fd
lines(myfd_c9)
myfd_c10 = smooth.basis(samples, Pig10_wide$C10, myfdPar)$fd
lines(myfd_c10)
myfd_c11 = smooth.basis(samples, Pig10_wide$C11, myfdPar)$fd
lines(myfd_c11)
## Continuous registration??
lambda <- 1
nbasis <- myfd$basis$nbasis
ntrials <- dim(Pig10_wide)[2]
y0fd <- mean.fd(myfd)
yfd = myfd
y0vec = eval.fd(samples, y0fd)
yvec <- eval.fd(samples, yfd)
coef0 <- matrix(0, nrow = nbasis, ncol = ntrials)
Wfd0 <- fd(coef0, mybasis)
WfdPar <- fdPar(Wfd0, 2, lambda)
reglist <- register.fd(y0fd, yfd, WfdPar, iterlim = 10, dbglev = 1)


