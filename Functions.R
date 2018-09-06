library(readr); library(pastecs); library(stringr); library(data.table)
library(dplyr); library(fda); library(tidyr)

## Function to read in data and time normalize
identify_cycles = function(data, axis, cutoff = -10){
  
  data = fread(data)
  colnames(data) = c('frame','tx','ty','tz','rx','ry','rz')
  data = data[,.(frame, get(axis))]
  colnames(data) = c('frame', axis)
  
  
  turn_points = turnpoints(data[[axis]])
  peaks = extract(turn_points)
  
  data[, turnpoints := peaks]
  data[get(axis) > cutoff, turnpoints := 0]
  
  ## Add column for index to maintain order
  
  index = 1:nrow(data)
  data$Index = index
  
  
  data$Cycle = NA
  cycle = 0
  cycle_start = (data[turnpoints == -1])[1]
  
  for(i in 1:length(data[[axis]]) ){
    if(data$turnpoints[i] == -1){
      cycle = cycle + 1
      cycle_start = rbind(cycle_start, data[i])
    }
    data$Cycle[i] = cycle
  }
  
  ## Only select unique cycles
  cycle_start = unique(cycle_start)
  ## Add Cycles to cycles_start
  cycle_start$Cycle = (1:nrow(cycle_start)) - 1
  
  ## Bind data with cycle start, adding the duplicate rows
  ## That is, so cycle 2 starts at end point of cycle 1 etc. 
  #data = rbind(cycle_start, data)
  #data = data[order(Index)]
  return(data)
}

initial_plot = function(cycle_data, axis){
  
  min = cycle_data[turnpoints == -1]


  ggplot(cycle_data, aes(x = frame, y = cycle_data[[axis]]))+
    geom_line() +
    geom_vline(xintercept = min$frame, col = 'red') +
    geom_point(data = min, aes(x = frame, y = min[[axis]]), col = 'red')
  
}

data = "pig10_20150826_ctrl_carrot_1_TMJBUTTER25_sm250.csv"
axis = 'TMJdata1.rz'

cycle_data = identify_cycles(data, axis)
initial_plot(cycle_data, axis)

################################################
##
##        Time normalizing function           ##
##
#################################################

## Function to interprolate data to make all cycles the same length
## Will call later in our time normalizing function
interpolation_func = function(data, cycle, max_cycle, axis){
  max_cyc = max_cycle
  data = as.data.table(data)
  data = data[Cycle == cycle][,.(frame, get(axis))]
  colnames(data) = c('frame',axis)
  n = max_cyc-dim(data)[1]
  interp = as.data.table(approx(data$frame, data[[axis]], n = n+2 ))
  colnames(interp) = c('frame', axis)
  data2 = rbind(interp, data)
  data2 = data2[order(frame)][,Cycle := cycle]
  data2 = data2[2:(max_cyc+1)]
  return(data2)
}

time_normalize = function(data, axis, cycles){

  data = identify_cycles(data, axis)
  ## Keep only those cycles we deem good
  data = data[Cycle %in% cycles,]
  
  ##Calculate cycle lengths
  cycle_lengths = data[,.N, by = Cycle]
  max_cycle = max(cycle_lengths)
  
  ## Now we want to normalize so
  ## Each cycle to start at frame 0
  data = data[,frame := frame - min(frame), by = Cycle]
  ##Each cycle ends at frame 1
  data = data[, frame := frame/max(frame), by = Cycle]
  
  ##call Interprolation function on all cycles
  data_inter_list = list()
  for(i in cycles){
    data_inter_list[[i]] = interpolation_func(data, i, max_cycle, axis)
    
  }
  
  data = rbindlist(data_inter_list)
  
  

}

data = "pig10_20150826_ctrl_carrot_1_TMJBUTTER25_sm250.csv"
axis = 'TMJdata1.rz'

cycles = c(1,2,3,4,5,7,8,9,10,11)
test1 = time_normalize(data, axis, cycles)

saturated_fd = function(data, axis, cycles, lambda = 1e-12,
                        norder = 6){
  data = time_normalize(data, axis, cycles)
  cycle_lengths = data[,.N, by = Cycle]
  max_cycle = max(cycle_lengths)
  n_cycles = length(cycles)
  
  id = rep(1:max_cycle, n_cycles )
  data[,id := id]
  data = data[,.(id, Cycle, get(axis))]
  colnames(data) = c('id', 'Cycle',axis)
  data_wide = dcast(data, id ~ Cycle, value.var = axis )
  
  colnames(data_wide) = c('id', paste0('C', cycles))
  
  data_wide = data_wide[,!c('id'), with = FALSE]
  
  #split 0 to 1 into 84 values (length of our cycles )
  samples = seq(0,1,length = 84)
  # Define the number of basis functions we want
  nbasis = length(samples) + norder -2
  # Create basic basis functions with our parameters
  mybasis = create.bspline.basis(c(0,1), nbasis, norder, samples)
  myfdPar = fdPar(mybasis, 4, lambda)
  # Dataset has to be in matrix form
  dataset = as.matrix(data_wide)
  ## Lets make it a fd object based off of our basis functions! woohoo! 
  myfd = smooth.basis(samples, dataset, myfdPar)$fd
  
  
  return(myfd)
}

data = "pig10_20150826_ctrl_carrot_1_TMJBUTTER25_sm250.csv"
axis = 'TMJdata1.rz'
cycles = c(1,2,3,4,5,7,8,9,10,11)

myfd1 = saturated_fd(data, axis, cycles)

landmark_reg = function(raw_data, axis, cycles){
  fd_obj = saturated_fd(raw_data, axis, cycles)
  raw_data = time_normalize(raw_data, axis, cycles)
  
  #fd_obj = deriv.fd(fd_obj, 2)
  
  ## Landmark Registration
  
  #First take the splines created from our fd object 
  #and put into separate data.table
  splines = as.data.table(fd_obj$coefs)
  
  #Create an empty list to put things in later
  peak_x = list()
  #These are the x values on the plots
  x_vals = as.matrix(fd_obj$basis$params)
  
  #Loop through all columns in our splines list
  for(col in colnames(splines)){
    print(col)
    c = splines %>% select(col)
    
    #split colname string so we can select appropriate cycle
    cyc = as.numeric(strsplit(col, 'C')[[1]][2])
    #Select frames for corresponding cycle
    print(cyc)
    frames = raw_data[Cycle == cyc,][,.(frame)]
    print(frames)
    #Add column with spline points
    frames = frames[,spline := c]
    #Now we just want the x value (frame) corresponding to maximum spline point
    peak_x[[cyc]] = frames[which.max(frames$spline)]$frame
  }
  #Unlist our peak_x to use late
  peak_x = unlist(peak_x)
  print(peak_x)
  #Take the mean of it
  mean_peak = mean(peak_x)
  
  #Create some W functions and stuff
  #Still not 100% certain what these do
  wbasisLM = create.bspline.basis(c(0,1), 4, 3, c(0, mean_peak,1))
  WfdLM    = fd(matrix(0,4,1),wbasisLM)
  WfdParLM = fdPar(WfdLM,1,1e-12)
  
  #Throw it all into the landmarkreg function
  #peak_x are the points we are landmarking (max vals)
  regListLM = landmarkreg(fd_obj, peak_x,mean_peak, WfdParLM, TRUE) 
}

data = "pig10_20150826_ctrl_carrot_1_TMJBUTTER25_sm250.csv"
axis = 'TMJdata1.rz'
cycles = c(1,2,3,4,5,7,8,9,10,11)

test = landmark_reg(data, axis, cycles)
