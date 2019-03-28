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

## Use the source for common functions
## This file must either be in current working directory 
## or full file destination must be specified
source('FDA_CoCmmon_funcs.R')

## Read our interprolated data in
pigs = fread('pig_cycles_int_full.csv')

## Add a column indicating place in cycle
## Cycle length = 163
## # of cycles = 1122
id = rep(1:163, 1122)
pigs[,cyc_id := id]

## Use data_wide_func on all pigs
pig10 = data_wide_func(pigs, animal_num = 10, food = 'carrot', axis = 'X.TMJdata.rz')
pig5 = data_wide_func(pigs, animal_num = 5, food = 'carrot', axis = 'X.TMJdata.rz')
pig6 = data_wide_func(pigs, animal_num = 6, food = 'carrot', axis = 'X.TMJdata.rz')
pig9 = data_wide_func(pigs, animal_num = 9, food = 'carrot', axis = 'X.TMJdata.rz')
for(f in c('apple','almond')){
  new_df= data_wide_func(pigs, animal_num = 10, food = f, axis = 'X.TMJdata.rz')
  pig10 = cbind(pig10, new_df)
  
  new_df= data_wide_func(pigs, animal_num = 5, food = f, axis = 'X.TMJdata.rz')
  pig5 = cbind(pig5, new_df)
  
  new_df= data_wide_func(pigs, animal_num = 6, food = f, axis = 'X.TMJdata.rz')
  pig6 = cbind(pig6, new_df)
  
  new_df= data_wide_func(pigs, animal_num = 9, food = f, axis = 'X.TMJdata.rz')
  pig9 = cbind(pig9, new_df)
}

## Full wide-data set (for axis 'X.TMJdata.rz')
full_pigs = cbind(pig5, pig6, pig9, pig10)

animal_fda = fda_func(full_data, 20)


## make the functional data object for our FULL data. This is what we will use
  ## to create our Overall Mean that we are going to compare against
full_fda = fda_func(full_pigs, 20)

plot(mean(full_fda), ylim = c(-20, 0))

## plot_type = 1 will give the overall mean

animal_reg = list()
for(n in c(5, 6, 9, 10)){
  animal_reg[[n]] = plt_func_a(fda_data = full_fda, full_data= full_pigs, animal_num = n, type = 'fda', 
                               stderrList = NULL,
                               ylim = NULL, plot_type = 1)
  
}

## plot overall_mean
plotbeta_n(animal_reg[[5]][[1]], animal_reg[[5]][[2]], plot_type = 1, ylim = c(-20, 0))

## plot how each animal now compare to this overall mean...
  ## 9/10 have waaay bigger CIs, not sure why...
plotbeta_n(animal_reg[[5]][[1]], animal_reg[[5]][[2]], plot_type = 2, ylim = c(-6, 6))

plotbeta_n(animal_reg[[6]][[1]], animal_reg[[6]][[2]], plot_type = 2, ylim = c(-6, 6))

plotbeta_n(animal_reg[[9]][[1]], animal_reg[[9]][[2]], plot_type = 2, ylim = c(-6, 6))

plotbeta_n(animal_reg[[10]][[1]], animal_reg[[10]][[2]], plot_type = 2, ylim = c(-6, 6))


##################### 5 ##########################################
five_fda = fda_func(pig5, 20)

five = plt_func_f(fda_data = five_fda, full_data = pig5, food = 'carrot', type = 'fda', 
                               stderrList = NULL,
                               ylim = NULL, plot_type = 1)

##################### 6 ##########################################
six_fda = fda_func(pig6, 20)

six = plt_func_f(fda_data = six_fda, full_data = pig6, food = 'carrot', type = 'fda', 
                 stderrList = NULL,
                 ylim = NULL, plot_type = 1)

##################### 9 ##########################################
nine_fda = fda_func(pig9, 20)

nine = plt_func_f(fda_data = nine_fda, full_data = pig9, food = 'carrot', type = 'fda', 
                 stderrList = NULL,
                 ylim = NULL, plot_type = 1)

##################### 10 ##########################################
ten_fda = fda_func(pig10, 20)

ten = plt_func_f(fda_data = ten_fda, full_data = pig10, food = 'carrot', type = 'fda', 
                 stderrList = NULL,
                 ylim = NULL, plot_type = 1)

plotbeta_n(five[[1]],five[[2]], plot_type = 1, ylim = c(-20, 0))
plotbeta_n(six[[1]],six[[2]], plot_type = 1, ylim = c(-20, 0))
plotbeta_n(nine[[1]],nine[[2]], plot_type = 1, ylim = c(-20, 0))
plotbeta_n(ten[[1]],ten[[2]], plot_type = 1, ylim = c(-20, 0))


## Run the below to save plots into single pdf file
pdf('pigs_byanimal.pdf', width = 20, height = 30)
#jpeg('pigs_10_9_6_5.jpg')
par(mfrow=c(4,2), cex = 1.5)

plotbeta_n(five[[1]],five[[2]], plot_type = 1, ylim = c(-20, 0))
plotbeta_n(animal_reg[[5]][[1]], animal_reg[[5]][[2]], plot_type = 2, ylim = c(-6, 6))

plotbeta_n(six[[1]],six[[2]], plot_type = 1, ylim = c(-20, 0))
plotbeta_n(animal_reg[[6]][[1]], animal_reg[[6]][[2]], plot_type = 2, ylim = c(-6, 6))

plotbeta_n(nine[[1]],nine[[2]], plot_type = 1, ylim = c(-20, 0))
plotbeta_n(animal_reg[[9]][[1]], animal_reg[[9]][[2]], plot_type = 2, ylim = c(-6, 6))

plotbeta_n(ten[[1]],ten[[2]], plot_type = 1, ylim = c(-20, 0))
plotbeta_n(animal_reg[[10]][[1]], animal_reg[[10]][[2]], plot_type = 2, ylim = c(-6, 6))


dev.off()


