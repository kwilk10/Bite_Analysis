
library(readr); library(pastecs); library(stringr); library(data.table)
library(dplyr); library(fda)
#setwd('~/Desktop/FDA')
setwd('~/Desktop/FDA')

source('FDA_Common_funcs.R')

## So for interactions, we want to compare to the OVERALL mean
## So we need the full data!! 
## And then subset by animal and food
  ## Should end up with 12 interaction plots that show how each animal/food combo
    # differs from the overall mean 


## Read our interprolated data in
pigs = fread('pig_cycles_int_full.csv')


## Call pigs_working_data function from FDA_common_funcs.R
##  To get working data set 
##  Each column is now a single cycle associated with a pig and food
full_pigs = pigs_working_data(pigs, 'X.TMJdata.rz')

## Now we need to create our regression based off of this 

full_fda = fda_func(full_pigs, 20)

c5_stdr = readRDS('by_food/stderrList_fda_carrot5')
ap5_stdr = readRDS('by_food/stderrList_fda_apple5')
al5_stdr = readRDS('by_food/stderrList_fda_almond5')


apple_5 = int_plt(fda_data = full_fda, full_data = full_pigs, animal_num = 5, 
               food = 'apple', type = 'fda', 
                   stderrList = NULL)

carrot_5 = int_plt(fda_data = full_fda, full_data = full_pigs, animal_num = 5, 
        food = 'carrot', type = 'fda',
        stderrList = NULL)

almond_5 = int_plt(fda_data = full_fda, full_data = full_pigs, animal_num = 5, 
        food = 'almond', type = 'fda',
        stderrList = NULL)

plotbeta_n(apple_5[[1]], apple_5[[2]], plot_type = 2)
plotbeta_n(carrot_5[[1]], carrot_5[[2]], plot_type = 2)
plotbeta_n(almond_5[[1]], almond_5[[2]], plot_type = 2)


############################################################################################
c6_stdr = readRDS('by_food/stderrList_fda_carrot6')
ap6_stdr = readRDS('by_food/stderrList_fda_apple6')
al6_stdr = readRDS('by_food/stderrList_fda_almond6')


apple_6 = int_plt(fda_data = full_fda, full_data = full_pigs, animal_num = 6, 
        food = 'apple', type = 'fda', stderrList = NULL)

carrot_6 = int_plt(fda_data = full_fda, full_data = full_pigs, animal_num = 6, 
        food = 'carrot', type = 'fda', 
        stderrList = NULL)

almond_6 = int_plt(fda_data = full_fda, full_data = full_pigs, animal_num = 6, 
        food = 'almond', type = 'fda',
        stderrList = NULL)

plotbeta_n(apple_6[[1]], apple_6[[2]], plot_type = 2)
plotbeta_n(carrot_6[[1]],carrot_6[[2]], plot_type = 2)
plotbeta_n(almond_6[[1]], almond_6[[2]], plot_type = 2)

############################################################################################
c9_stdr = readRDS('by_food/stderrList_fda_carrot9')
ap9_stdr = readRDS('by_food/stderrList_fda_apple9')
al9_stdr = readRDS('by_food/stderrList_fda_almond9')


apple_9 = int_plt(fda_data = full_fda, full_data = full_pigs, animal_num = 9, 
        food = 'apple', type = 'fda', 
        stderrList = NULL)

carrot_9 = int_plt(fda_data = full_fda, full_data = full_pigs, animal_num = 9, 
        food = 'carrot', type = 'fda',
        stderrList = NULL)

almond_9 = int_plt(fda_data = full_fda, full_data = full_pigs, animal_num = 9, 
        food = 'almond', type = 'fda', 
        stderrList = NULL)

plotbeta_n(apple_9[[1]], apple_9[[2]], plot_type = 2)
plotbeta_n(carrot_9[[1]],carrot_9[[2]], plot_type = 2)

plotbeta_n(almond_9[[1]], almond_9[[2]], plot_type = 2)


############################################################################################
c10_stdr = readRDS('by_food/stderrList_fda_carrot10')
ap10_stdr = readRDS('by_food/stderrList_fda_apple10')
al10_stdr = readRDS('by_food/stderrList_fda_almond10')


apple_10 = int_plt(fda_data = full_fda, full_data = full_pigs, animal_num = 10, 
        food = 'apple', type = 'fda', 
        stderrList = NULL)

carrot_10 = int_plt(fda_data = full_fda, full_data = full_pigs, animal_num = 10, 
        food = 'carrot', type = 'fda', 
        stderrList = NULL)

almond_10 = int_plt(fda_data = full_fda, full_data = full_pigs, animal_num = 10, 
        food = 'almond', type = 'fda', 
        stderrList = NULL)

plotbeta_n(apple_10[[1]], apple_10[[2]], plot_type = 2)
plotbeta_n(carrot_10[[1]],carrot_10[[2]], plot_type = 2)
plotbeta_n(almond_10[[1]], almond_10[[2]], plot_type = 2)



############################################################################################
##################################################################################
##################################################################################

pdf('pigs_interactions.pdf', width = 20, height = 30)
#jpeg('pigs_10_9_6_5.jpg')
par(mfrow=c(4,3), cex = 1.5)

plotbeta_n(apple_5[[1]], apple_5[[2]], plot_type = 2, ylim = c(-8, 3) )
plotbeta_n(carrot_5[[1]], carrot_5[[2]], plot_type = 2, ylim = c(-8, 3))
plotbeta_n(almond_5[[1]], almond_5[[2]], plot_type = 2, ylim = c(-8, 3))



plotbeta_n(apple_6[[1]], apple_6[[2]], plot_type = 2, ylim = c(-8, 3))
plotbeta_n(carrot_6[[1]],carrot_6[[2]], plot_type = 2, ylim = c(-8, 3))
plotbeta_n(almond_6[[1]], almond_6[[2]], plot_type = 2, ylim = c(-8, 3))


plotbeta_n(apple_9[[1]], apple_9[[2]], plot_type = 2, ylim = c(-8, 3))
plotbeta_n(carrot_9[[1]],carrot_9[[2]], plot_type = 2, ylim = c(-8, 3))
plotbeta_n(almond_9[[1]], almond_9[[2]], plot_type = 2, ylim = c(-8, 3))

plotbeta_n(apple_10[[1]], apple_10[[2]], plot_type = 2, ylim = c(-8, 3))
plotbeta_n(carrot_10[[1]],carrot_10[[2]], plot_type = 2, ylim = c(-8, 3))
plotbeta_n(almond_10[[1]], almond_10[[2]], plot_type = 2, ylim = c(-8, 3))

dev.off()

##################################################################################
##################################################################################
##################################################################################

## Okay, now we need to run the above for each of the 12 different combinations and put them into a file..
pdf('interactions.pdf',width = 20, height = 30)

par(mfrow=c(4,3), cex = 1.5)
par(ask=F)
## pig 5
plt_func_f2(fda_data = apple_fda, full_data = apple_wide,
            animal_num = 5, type = 'fda')

par(ask=F)

plt_func_f2(fda_data = carrot_fda, full_data = carrot_wide,
            animal_num = 5, type = 'fda')
par(ask=F)
plt_func_f2(fda_data = almond_fda, full_data = almond_wide,
            animal_num = 5, type = 'fda')
par(ask=F)
## pig 6
plt_func_f2(fda_data = apple_fda, full_data = apple_wide,
            animal_num = 6, type = 'fda')
par(ask=F)
plt_func_f2(fda_data = carrot_fda, full_data = carrot_wide,
            animal_num = 6, type = 'fda')
par(ask=F)
plt_func_f2(fda_data = almond_fda, full_data = almond_wide,
            animal_num = 6, type = 'fda')
par(ask=F)

## pig 9

plt_func_f2(fda_data = apple_fda, full_data = apple_wide,
            animal_num = 9, type = 'fda')
par(ask=F)
plt_func_f2(fda_data = carrot_fda, full_data = carrot_wide,
            animal_num = 9, type = 'fda')
par(ask=F)
plt_func_f2(fda_data = almond_fda, full_data = almond_wide,
            animal_num = 9, type = 'fda')
par(ask=F)

## pig 10


plt_func_f2(fda_data = apple_fda, full_data = apple_wide,
            animal_num = 10, type = 'fda')
par(ask=F)
plt_func_f2(fda_data = carrot_fda, full_data = carrot_wide,
            animal_num = 10, type = 'fda')
par(ask=F)
plt_func_f2(fda_data = almond_fda, full_data = almond_wide,
            animal_num = 10, type = 'fda')

dev.off()



pdf('Almond_Carrot_apple.pdf',width = 20, height = 30)
par(mfrow=c(3,2), cex = 1.5)
par(ask=F)

par(ask=F)
plt_func_f2(fda_data = almond_fda, full_data = almond_wide,
            animal_num = 10, type = 'fda')
par(ask=F)
plt_func_f(full_fda, full_data,'full_fda','almond', al_fd)

plt_func_f2(fda_data = carrot_fda, full_data = carrot_wide,
            animal_num = 10, type = 'fda')
par(ask=F)
plt_func_f(full_fda, full_data,'full_fda','carrot', cr_fd)

par(ask=F)
plt_func_f2(fda_data = apple_fda, full_data = apple_wide,
            animal_num = 10, type = 'fda')
par(ask=F)
plt_func_f(full_fda, full_data,'full_fda','apple', ap_fd)

dev.off()

