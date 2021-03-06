# Settables are on the bottom in the main run_all function.

# Libraries
library(zoo) # na.locf, na.approx

# Helper functions
animal_from_filename <- function(name) {
  animal <- unlist(strsplit(name, c('-')))[2]
  # Strip plurality  TODO: make sure this works for all the aminals (sic).
  if (substr(animal, nchar(animal), nchar(animal)) == 's') {
    animal <- substr(animal, 1, nchar(animal) - 1)
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
    par_frame <- read.csv(paste0(path, param))
    par_frame[par_frame$Food == '',]$Food <- NA
    par_frame$Food <- na.locf(par_frame$Food)
    par_frame$ID <- na.locf(par_frame$ID)
    seq_only <- par_frame[as.character(par_frame$Date) != "",]

    animal <- animal_from_filename(param)
    
    # iterate over rows in parameter file in different data files
    for (rw in 1:nrow(seq_only)) {
      # construct filename to read:
      row <- seq_only[rw,]
      
      datestamp <- format(as.Date(row$Date, '%m/%d/%Y'), format = '%Y%m%d')
      
      individual <- row$Individual
      food <- tolower(as.character(row$Food))
      id <- row$ID
      treated <- ifelse(row$Treatment == 0, "ctrl", "treated") #TODO: find out what "treated" actually is in the filenames
      
      regex <- paste0(animal, individual, "_", datestamp, "_", treated, "_", food, "_", id, '.*\\.csv')
      datfile <- list.files(path = datapath, pattern = regex)
      
      if (length(datfile) == 0) next
      
      # read it in
      newframe <- read.csv(paste0(datapath,  datfile))

      # Flag the heck out of it!
      newframe$animal <- animal
      newframe$animal_id <- individual
      newframe$food <- food
      newframe$sequence <- row$Sequence
      newframe$date <- datestamp

      cycle_frame <- subset(par_frame,
                            Sequence == row$Sequence & Food == row$Food & Individual == row$Individual,
                            c(start, end))

      vec <- numeric(nrow(newframe))
      for (i in 1:nrow(cycle_frame)) {
        # Making cycle numbers unique in case that makes something easier later
        cyc <- ifelse(exists("ret_df"), max(ret_df$orig_cycle) + i, i)
        rang <- c(cycle_frame[i, "start"]:cycle_frame[i, "end"])
        vec[rang] <- cyc
      }
      newframe$orig_cycle <- vec

      if  (!exists("ret_df")) {
        ret_df <- newframe
      }
      else {
        ret_df <- rbind(ret_df, newframe)
      }
    }
  }
  ret_df
}

localextrema <- function(x, lookbehind = 1, threshold = 0, 
                         type = 'minima', scaled = FALSE, stddev = 0) {
  if (type == 'minima') {
    signdiff <- 2
  }
  else if (type == 'maxima') {
    signdiff <- -2
  }
  else {
    stop(c("Invalid type: ", type, " choose either 'maxima' or 'minima'"))
  }
  
  if (scaled == TRUE) {
    vec <- scale(x)
  }
  else {
    vec <- x
  }
  
  changes <- which(diff(sign(diff(vec))) == signdiff) + 1
  if (missing(stddev) | scaled == FALSE) {
    changes[vec[changes - lookbehind] - vec[changes] > threshold]
  }
  else {
    if (type == 'minima') {
      changes[vec[changes] < stddev]
    }
    else {
      changes[vec[changes] > stddev]
    }
  }
}

real_minima <- function(x, cname) {
  indices <- localextrema(x[,cname], scaled = TRUE, stddev = -1)
  x$real_cycle <- numeric(nrow(x))

  for (i in indices) {
    x[i, "real_cycle"] <- match(i, indices)
  }
  
  x[x$real_cycle == 0, "real_cycle"] <- NA
  x$real_cycle <- na.locf(x$real_cycle, na.rm = FALSE)
  x
}

# Main function
run_all <- function() {
  # Settables
  PARENT <- '/home/cgerstner/dad/'
  PARAM_DIR <- paste0(PARENT, 'params/')
  DATA_DIR <- paste0(PARENT, 'data/')
  SAVE_DIR <- paste0(PARENT, 'prepared/')

  df <- read_data(PARAM_DIR, DATA_DIR)
  df <- real_minima(df, "TMJdata1.rz")

  save(df, file = paste0(SAVE_DIR, "data.Rdata"))
  df
}


