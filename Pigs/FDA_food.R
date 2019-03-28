##############################
# Run some FDA by Food Type          #
#
# This script takes the interprolated data from the pigs_inter.R script
# And creates funcational data objects and runs fANOVA by food type
# 
#
#####################################


library(readr); library(pastecs); library(stringr); library(data.table)
library(dplyr); library(fda)
#setwd('~/Desktop/FDA')
setwd('~/Desktop/FDA')

source('FDA_Common_funcs.R')
## Read our interprolated data in
pigs = fread('pig_cycles_int_full.csv')

## Call pigs_working_data function from FDA_common_funcs.R
##  To get working data set 
##  Each column is now a single cycle associated with a pig and food
full_pigs = pigs_working_data(pigs, 'X.TMJdata.rz')

## First off, we are going to create our wide data and bind it all together into one master data.table
carrot_wide = full_pigs %>% dplyr::select(contains('carrot'))
almond_wide = full_pigs %>% dplyr::select(contains('almond'))
apple_wide = full_pigs %>% dplyr::select(contains('apple'))

food_reg = list()
for(f in c('almond','carrot','apple')){
  food_reg[[f]] = plt_func_f(fda_data = full_fda, full_data= full_pigs, food = f, type = 'fda', 
                               stderrList = readRDS(sprintf('by_animal/stderrList_%s_%s', 'fda', f) ) ,
                               ylim = NULL, plot_type = 1)
  
}



plotbeta_n(food_reg[['almond']][[1]],food_reg[['almond']][[2]], plot_type = 1, ylim = c(-20, 0))

carrot_fda = fda_func(carrot_wide, 20)
almond_fda = fda_func(almond_wide, 20)
apple_fda = fda_func(apple_wide, 20)

carrot = plt_func_a(fda_data = carrot_fda, full_data= carrot_wide, animal_num = 5, type = 'fda', 
           stderrList = NULL,
           ylim = NULL, plot_type = 1)
almond = plt_func_a(fda_data = almond_fda, full_data= almond_wide, animal_num = 5, type = 'fda', 
                    stderrList = NULL,
                    ylim = NULL, plot_type = 1)
apple = plt_func_a(fda_data = apple_fda, full_data= apple_wide, animal_num = 5, type = 'fda', 
                    stderrList = NULL,
                    ylim = NULL, plot_type = 1)



pdf('pigs_byfood.pdf', width = 20, height = 30)
#jpeg('pigs_10_9_6_5.jpg')
par(mfrow=c(4,2), cex = 1.5)


plotbeta_n(almond[[1]],almond[[2]], plot_type = 1, ylim = c(-20, 0))
abline(h = c(-10, -15,-10, -5, 0), v = c( 0, 0.2, 0.4, 0.6, 0.8, 1.0), lty = 2, lwd = 0.5)
plotbeta_n(food_reg[['almond']][[1]],food_reg[['almond']][[2]], plot_type = 2, ylim = c(-1.5, 1.5))

plotbeta_n(carrot[[1]],carrot[[2]], plot_type = 1, ylim = c(-20, 0))
abline(h = c(-10, -15,-10, -5, 0), v = c( 0, 0.2, 0.4, 0.6, 0.8, 1.0), lty = 2, lwd = 0.5)
plotbeta_n(food_reg[['carrot']][[1]],food_reg[['carrot']][[2]], plot_type = 2, ylim = c(-1.5, 1.5))

plotbeta_n(apple[[1]],apple[[2]], plot_type = 1, ylim = c(-20, 0))
abline(h = c(-10, -15,-10, -5, 0), v = c( 0, 0.2, 0.4, 0.6, 0.8, 1.0), lty = 2, lwd = 0.5)
plotbeta_n(food_reg[['apple']][[1]],food_reg[['apple']][[2]], plot_type = 2, ylim = c(-1.5, 1.5))

dev.off()
