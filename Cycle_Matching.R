#install.packages('lubridate')
#install.packages('zoo')

library(lubridate); library(zoo); library(stringr); library(data.table)

## Read in master cycle data sheet
keep = c('Individual','Date','Food','ID','Cycle','start','end')
cycle_id = fread("NIH-pigs-5-6-9-10_Ctrl_2018-08-28_for-Geoff.csv", select = keep)


##First we are going to just look at pig 10 so lets filter out to that
pig10_cycles = cycle_id[Individual == 10]

##Read in pig 10 data and run identify_cycles function on it
data = "pig10_20150827_ctrl_Apple_1_TMJBUTTER25_sm250.csv"
axis = 'rz'
pig10_apple_1 = identify_cycles(data, axis)

## Ignore this ##
date = str_extract(data,"[0-9]{8}")
food = str_extract(data,"[A-Za-z]{5,6}_\\d{1}")

## Want all blank spaces to be NA
pig10_cycles$Date[pig10_cycles$Date == ''] <- NA
pig10_cycles$Food[pig10_cycles$Food == ''] <- NA
pig10_cycles$ID[pig10_cycles$ID == ''] <- NA

## This adds data into the NA slots 
pig10_cycles[,`:=`(Date=na.locf(Date), Food = na.locf(Food), ID = na.locf(ID))]
## Change the Date column to a date. Not super necessary step yet but will be nice to later match
pig10_cycles[,Date := mdy(Date)]


## We just want to start with matching Pig 10 Apple 1 on 08/27/15
pig10_cycles_ap1 = pig10_cycles[Food == 'Apple' & ID ==1 & Date == '2015-08-27']


## Loop through and add a new Cycle column that has the cycles as labled in the NIH file
for(i in 1:nrow(pig10_cycles_ap1)){
  start = pig10_cycles_ap1$start[i]
  end = pig10_cycles_ap1$end[i]
  cyc = pig10_cycles_ap1$Cycle[i]
  
  pig10_apple_1[Index %in% start:end, Cycle_2 := cyc]
}


## Check it ouuut
View(pig10_apple_1)
