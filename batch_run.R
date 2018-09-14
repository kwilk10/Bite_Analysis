# Settables -- Right up top

PARENT <- '/home/cgerstner/dad/'
PARAM_DIR <- paste0(PARENT, 'params/')
DATA_DIR <- paste0(PARENT, 'data/')
SAVE_DIR <- paste0(PARENT, 'prepared/')

# Libraries
library(zoo)

# Helper functions
animal_from_filename <- function(name) {
  animal <- strsplit(name, c('-'))[2]
  # Strip plurality  TODO: make sure this works for all the aminals (sic).
  if (substr(animal, nchar(animal), nchar(animal)) == 's') {
    animal <- substr(animal, 1, nchar(tst) - 1)
  }
  animal
}

# Function to:
# * read parameter files
# * generate filenames from data in parameter files
# * read data files and combine with classifications from parameters
# returns merged, classified data frame

# Can't shake the feeling I could have done this better, vectorized and more R-like, but it's fine.

# TODO: test at all
read_data <- function(path, datapath) {
  # Here we pull in all of the descriptor files from which we shall pull
  # in all the data files.
  for (param in list.files(path)) {
    par_frame <-read.csv(paste0(path, param))
    par_frame[par_Frame$Food == '',]$Food <- NA
    par_frame$Food <- na.locf(par_frame$Food)
    seq_only <- par[as.character(par$Date) != ""]
    
    # iterate over rows in parameter file in different data files
    for (row in 1:nrow(seq_only)) {
      # construct filename to read:
      
      # convert date format in parameter file to Ymd like filenames
      datestamp <- format(as.Date(newframe$Date, '%m/%d/%Y'), format = '%Y%m%d')
      # animal
      animal <- animal_from_filename(param)
      # subject #
      individual <- row$Individual
      # fewd
      food <- tolower(as.character(row$Food))
      # treatmentity
      treated <- ifelse(row$Treatment == 0, "ctrl", "treated") #TODO: find out what "treated" actually is in the filenames
      
      # regular expression match filename for this sequence
      regex <- paste0(animal, individual, "_", datestamp, "_", treated, "_", food, '.*\\.csv')
      datfile <- list.files(path = datapath, pattern = regex)
      
      # read it in
      newframe <- read.csv(paste0(datapath,  datfile))

      # Flag the heck out of it!
      newframe$animal <- animal
      newframe$food <- food
      newframe$sequence <- row$Sequence
      newframe$date <- datestamp
      
      cycle_frame <- subset(par_frame, 
                            Sequence == row$Sequence & Food == food & Individual == individual, 
                            c(start, end))
      for (i in nrow(cycle_frame)) {
        newframe[cycle_frame[i,"start"]:cycle_frame[i,"end"],]$orig_cycle <- i
      }
      
      # rbind, with minor trade of efficiency for readability
      df <- ifelse(exists(df), newframe, rbind(df, newframe))
    }
  }
  df
}



