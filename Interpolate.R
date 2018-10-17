##############################
# Interprolate!           #
#
# This script take our cycle flagged dataset 
# As created in our Cycle_Matching.R script and interprolates data for each cycle
# So that all the cycles are the same length
# 
# * Note: This is all done linearly which may need to be adjusted later down the road
# 
# This script outputs a csv file with the interprolated data
# ** Currently only done on .rz axis
#
#####################################

library(readr); library(pastecs); library(stringr); library(data.table)
library(dplyr); library(fda)

## Set this working directory to wherever you saved the full csv file
setwd('~/Desktop/FDA')

# Read in the data
pigs = fread('pig_cycles_full.csv')

## count := uniqueN(orig_cycle) adds a count by animal, animal_id, food, date, and ID
pigs[, count := uniqueN(orig_cycle), by = .(animal, animal_id, food, date, ID)]

## uniq_id := .GRP adds a unique identifier to each cycle
# This is used to pull each cycle individually and add in the necessary data
pigs_test = pigs[, uniq_id := .GRP, .(animal, animal_id, food, date, ID, orig_cycle)]
uniq_ids = unique(pigs_test[,.(animal, animal_id, food, date, ID, uniq_id)])
View(pigs_test)

## This is our interprolate function. Inputs our the data file and the cycle number
inter_func = function(data, cycle){
  # Max cycle for all the pigs. Will need to be adjusted for other animals
  max_cyc = 163
  data = as.data.table(data)
  # Select just the cycle we want to add data to 
  # Currently just keep the .rz axis. Will need to add other axis later
  data = data[orig_cycle == cycle][,.(ref_frame, X.TMJdata.rz)]
  
  n = max_cyc-dim(data)[1]
  
  ## This is the actual interprolation step! 
  interp = as.data.table(approx(data$ref_frame, data$X.TMJdata.rz, n = n+2 ))
  
  colnames(interp) = c('ref_frame', 'X.TMJdata.rz')
  data2 = rbind(interp, data)
  data2 = data2[order(ref_frame)][,`:=`(Cycle = cycle, uniq_id = i)]
  data2 = data2[2:164]
  ## Return our interprolated data
  return(data2)
  
}


# Create empty list to put all our data into
int_data = list()
##Unique groups are: animal, animal_id food, date, ID, orig_cycle

## Loop through all the cycles (hence why we needed a unique identifier for each cycle)
for(i in unique(pigs_test$uniq_id)){
  
  # Select just the animal/food/ID we want
  data = pigs_test[uniq_id == i]
  
  # pull which cycle it actually is
  cycle = data$orig_cycle

    # Honestly not sure why I have this loop here but it is working currently 
    # and I don't feel like messing with it too much currently
      for(j in unique(data$orig_cycle)){
      id = sprintf('%i-%i', i, j)
      
      # This is where we actually run our inter_func from above
      # on our data and puts it into our empty list
      int_data[[id]] = inter_func(data, j)
      
    }
}
  
  
## bind that list all together into a single dataframe
pigs_int = rbindlist(int_data)    
View(pigs_int)

## Here we now merge it with our all the unique id's so that we have all our identifing information
# That is which animal, food, id, etc that it is
pigs_fin = merge(pigs_int, uniq_ids, by = 'uniq_id', all.y = TRUE)

## Finally, we can out put our final dataset to a csv! Woohoo! 
write.csv(pigs_fin, file = "/Users/maraudersmap/Desktop/FDA/pig_cycles_int.csv")
