## Analysis of pigs chewing cycles

### Scripts

Scripts should be run in the following order (FDA_food.R and FDA_by_animal.R are interchangeable)
 - [Cycle_Matching.R](./Cycle_Matching.R)
 - [Interprolate.R](./Interpolate.R)
 - [FDA_food.R](./FDA_food.R) or [FDA_by_animal.R](./FDA_by_animal.R)

#### [Cycle_Matching.R](./Cycle_Matching.R)

Takes raw data files found in the following format in working directory *data/pig#_ymd_ctrl_food_id_TMJBUTTER25_sm250.csv*. Labels cycles in raw data using the parameters file (labeled by Claire) under the following format, */params/pig_params.csv*. Combines all animals and foods into one master data file. Cycles are labeled by each unique animal_id, food, and date combination. That is, animal_id == 5, date = 20150319, food = apple has 14 labeled cycles while animal_id == 5, date = 20150319, food = almond has 6 labeled cycles. Outputs all labeled cycles and corresponding data into one .Rdata and .csv file.  

Outputs: 

    data.Rdata in prepared folder
    pig_cycles_full.csv in working directory
    
    
#### [Interprolate.R](./Interpolate.R)

This script uses the raw, labeled data, *pig_cycles_full.csv* as outputted by [Cycle_Matching.R](./Cycle_Matching.R) script. It interprolates data to each cycle so that each cycle is 163 frames long. It is currently set up just to work with X.TMJdata.rz axis (as of 12/28/2018). It again combines all food, animal, and data combinations into a single master file but only keeps the X.TMJdata.rz axis data (as of 12/28/2019). 

Outputs:
    
    pig_cycles_int.csv in working directory
 
#### [FDA_food.R](./FDA_food.R)

This script uses the raw, interprolated data outputed by [Interprolate.R](./Interprolate.R) script and runs some FDA by food type. Currently built only to work with X.TMJdata.rz axis. Runs some continuous registration and regular functional data analysis. 

#### [FDA_by_animal.R](./FDA_by_animal.R)

This script uses the raw, interprolated data outputed by [Interprolate.R](./Interprolate.R) script and runs some FDA by animal ID. Currently built only to work with X.TMJdata.rz axis. Runs some continuous registration and regular functional data analysis.  
