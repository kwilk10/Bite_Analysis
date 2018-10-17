##############################
# Flag Cycles as identified in parameters file (By Claire) #
#
# This script loops through the raw data files in your folder
# and flags the cycles as identified in a spearate parameters folder
# It makes a few assumptions about how the data files are stored
# on your local computer.
# First: Set your local working directory to a folder that contains
# a folder named each of the following:
#   1. params
#       - This contains the file that will be used to label the cycles
#   2. data
#       - This contains all the data files
#       - Here you should also have sub-folders for each animal and then food
#       - ex: data/pig_10/almond, data/pig_5/carrot, data/pig_6/apple, etc.
#   3. prepared
#       - This is just the folder that the final RData fill be outputted
#
# Finaly this script outputs an RData file (can also output a csv file) 
# That has all the animals and all the food types in a single file
#
# *Note: Frame work for all this work was pulled from Claytons batch_run.R script
# 
#
#####################################



# Libraries
library(zoo); library(data.table);library(lubridate) 


# This function will be used later to pull the animal type from the 
# parameter file
animal_from_filename <- function(name) {
  animal <- tolower(unlist(strsplit(name, c('_')))[1])
  # Strip plurality  TODO: make sure this works for all the aminals (sic).
  if (substr(animal, nchar(animal), nchar(animal)) == 's') {
    animal <- tolower(substr(animal, 1, nchar(animal) - 1))
  }
  print(animal)
}

# Function to:
# * read parameter files
# * generate filenames from data in parameter files
# * read data files and combine with classifications from parameters
# returns merged, classified data frame


read_data <- function(path, datapath) {
  # Here we pull in all of the descriptor files from which we shall pull
  # in all the data files.
  
  for (param in list.files(path)) {
    
    par_frame = fread(paste0(path,param))
    
    ## This is all just cleaning up the parameters data file
    par_frame[Food == '', Food := NA]
    par_frame[ID == '', ID := NA] 
    par_frame[Subject == '', Subject := NA]
    par_frame[,`:=`(Food = na.locf(Food), ID = na.locf(ID), Subject = na.locf(Subject))]
    par_frame[,Subject := as.character(Subject)]
    

    seq_only <- par_frame[as.character(par_frame$Date) != "",]
    par_frame[Date == '', Date := NA]
    par_frame[,Date:=na.locf(Date)]
    par_frame[,Date := ymd(Date)]
    seq_only[,Date := ymd(Date)]

    # Here we use our function from above to pull the animal type that we want
    animal <- animal_from_filename(param)
    
    #cyc = 1
    
    ## These create a few empty lists that we will fill later with our cycle flagged data
    full_data_list = list()
    cycframe_list_full = list()
    file_idx = 0
    
    # iterate over rows in parameter file in different data files
    for (rw in 1:nrow(seq_only)) {
      file_idx = file_idx + 1
      # construct filename to read:
      row <- seq_only[rw,]
      
      ## This formats the date column as we want it
      datestamp <- format(as.Date(row$Date, '%m/%d/%Y'), format = '%Y%m%d')
      
      # This pulls the animal subject number, food, and ID
      individual <- row$Subject
      food <- tolower(as.character(row$Food))
      id <- row$ID
      
      ## This was in the original parameters list but we didn't use it
      # will leave in case it become relevant
      #treated <- ifelse(row$Treatment == 0, "ctrl", "treated") #TODO: find out what "treated" actually is in the filenames
      
      ## This is just creating a string that will match our data file names
      regex <- paste0(animal, individual, "_", datestamp, "_", 'ctrl', "_", food, "_", id, "_",
                      'TMJBUTTER25_sm250', '.csv')
      # uncomment the following line to see what file name it is trying to find in your folder
      #print(regex)
      
      data = sprintf('%s/%s_%s/%s', datapath, animal, individual, food )
      
      # this now loops through your folder and puts all the files that match the regex into a list
      # There should only be one I believe
      datfile <- list.files(path = data, pattern = regex)
      
      # If the length of datfile == 0, then that file does not exist in your folder 
      # So we move on to the next argument and don't do anything with that particular regex string
      if (length(datfile) == 0) next
      
      # read it in
      
      ## We have now established that the file does exist in your folder so we can read it in! 
      newframe <- read.csv(paste0(data,'/',  datfile))
 
      
      # Flag the heck out of it!
      newframe$animal <- animal
      newframe$animal_id <- individual
      newframe$food <- food
      #newframe$sequence <- row$Sequence
      newframe$date <- datestamp
      newframe$ID <- id
      
      #cycle_frame <- subset(par_frame,
       #                     Date == row$Date & Food == row$Food & Subject == row$Subject,
        #                    c(`Frame Start`, `Frame end`))
      
      ## This subsets to just the date, food, subject, and id that we are working on (as chosen above)
      ## and selects just Frame Start and Frame End
      cycle_frame = par_frame[Date == row$Date & Food == row$Food & Subject == row$Subject & ID == row$ID,
                              .(start = `Frame Start`, end = `Frame end`)]
    
      ## A bunch of 0's, same number as number of rows as our read in file
      vec <- numeric(nrow(newframe))
      
      ## Add an empty column of NA's
      newframe$orig_cycle = NA
      
      # Start our cycle at 1 and create another empty list to fill
      cyc = 1
      cycframe_list = list()
      
      # Loop through the rows in our Frame Start and Frame End data set we created above
      for (i in 1:nrow(cycle_frame)) {
        
        # Making cycle numbers unique in case that makes something easier later
        #cyc <- ifelse(exists("ret_df"), max(ret_df$orig_cycle) + i, i)
      
        ## Here we pull the start and end of each cycle
        start = cycle_frame$start[i]
        end = cycle_frame$end[i]
        #Now select corresponding rows in full data and add cycle number
        newframe$orig_cycle[start:end] = cyc
        ## Add it to our empty list
        cycframe_list[[i]] = newframe[start:end,]
        # Advance cycle to next cycle number
        cyc = cyc + 1
        
        
        #rang <- start:end
        
        #vec[start:end] <- cyc
        
      }
      
      ## 
      #full_data_list[[file_idx]] = newframe
      
      # We want this one. Puts all the animals into their own part of the list
      cycframe_list_full[[file_idx]]=rbindlist(cycframe_list)
      
      #newframe$orig_cycle <- vec
      
      #if  (!exists("ret_df")) {
        #ret_df <- newframe
        #browser()
     # }
      #else {
       # ret_df <- rbind(ret_df, newframe)
      #}
    }
  }
  #full_data = rbindlist(full_data_list)
  #full_data[,ref_frame := round(frame)]
  
  ## Now we have looped through all the animals and food and there are all in our list
  ## So all we have to do is bind them together
  cycl_data = rbindlist(cycframe_list_full)
  # And finally, we add a column that is the frame number but rounded 
  # Just makes it easier to work with down the line
  cycl_data[,ref_frame := round(frame)]
  
}



# Main function
## This is the only function you actually run!! 
## The ONLY argument you need is your parent path 
## Aka the folder that contains your params, data, and prepared folder
run_all <- function(Parent_path) {
  # Settables
  PARENT = Parent_path
  PARAM_DIR <- paste0(PARENT, 'params/')
  print(PARAM_DIR)
  DATA_DIR <- paste0(PARENT, 'data/')
  SAVE_DIR <- paste0(PARENT, 'prepared/')
  
  df <- read_data(PARAM_DIR, DATA_DIR)
  #df <- real_minima(df, "TMJdata1.rz")
  
  #Here we save our file to our prepared folder as a RData file
  save(df, file = paste0(SAVE_DIR, "data.Rdata"))
}

## Running our function
Parent_path = '/Users/maraudersmap/Desktop/FDA/'
run_all(Parent_path)

## And now we can load our full, cycle labeled data! 
load('/Users/maraudersmap/Desktop/FDA/prepared/data.Rdata')
head(df)

## Here we can save it as a csv file to wherever we want so we can view it easier
write.csv(df, file = "/Users/maraudersmap/Desktop/FDA/pig_cycles_full.csv")


