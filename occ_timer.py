#####################################################################################################################################################
#########################################################  OCCULTATION TIMING PROGRAM  ##############################################################
#####################################################################################################################################################


#####################################################################################################################################################
###############################################################  INPUTS & OUTPUTS  ##################################################################
#####################################################################################################################################################


#########################################################################################################################################################################################################################  IMPORTS  ######################################################################
#####################################################################################################################################################

import numpy as np
import tables as tb
import os
from tqdm import trange, tqdm_notebook
from time import sleep

########################################################################################################################################################################################################################   CODE   ########################################################################
#####################################################################################################################################################


def occ_timer(path, n_source, thresh_factor):

    with open(path+str('timingData.txt'), 'r') as file:
        lines = [line.strip() for line in file.readlines()]
        numbers = [float(num) for num in lines]
        timingData = np.array(numbers, dtype=np.float64)

    time_array = timingData-timingData[0]
    
    o = open(path + str('occ_timer_results.txt'), "a")
    with open(path+str('photData/phot_source_')+str(n_source)+str('.txt'), 'r') as file:
        lines = [line.strip() for line in file.readlines()]
        numbers = [float(num) for num in lines]
        source_array = np.array(numbers, dtype=np.float64)
    mean = np.mean(source_array)
    threshold = thresh_factor*mean
        
    o.write(str('source_')+str(n_source)+str(' mean: ')+str(mean)+'\n')
    o.write(str('source_')+str(n_source)+str(' threshold: ')+str(threshold)+'\n'+'\n')
    o.close()
        
    for j in range(len(source_array)):
        if source_array[j] <= threshold and source_array[j-1] >= threshold:
            start_time = timingData[j]
        elif source_array[j] >= threshold and source_array[j-1] <= threshold:
            end_time = timingData[j]
        else:
            continue
                
    o = open(path + str('occ_timer_results.txt'), "a")
    o.write("Occultation Source: "+str(n_source)+'\n')
    o.write("Occultation Duration =  " +str(end_time-start_time)+str('s'))
    o.close()
    
    print('Occultation Source : ' +str(n_source))
    print("Occultation Duration =  " +str(end_time-start_time)+str('s'))
    