#####################################################################################################################################################
###########################################################  DATA REDUCTION PROGRAM  ################################################################
#####################################################################################################################################################

''' PYTHON VERSION: 3.12.1'''
'''Code to take a .fits image, remove the overscan region and reduce it by flat-fielding it, then save as a new fits file. '''

#####################################################################################################################################################
###############################################################  INPUTS & OUTPUTS  ##################################################################
#####################################################################################################################################################


########################################################################################################################################################################################################################  IMPORTS  #######################################################################
#####################################################################################################################################################

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import tables as tb
import os
import re

########################################################################################################################################################################################################################   CODE   ########################################################################
#####################################################################################################################################################

def reducto(im_file, flat_file, directory, make_dir):

    if make_dir == True:       # Make a directory called 'reduced for reduced images if required, if False, then code assumes that a directory with                                   that name exists already.
        os.mkdir(directory.replace('images/', 'reduced'))
        print(f"Directory 'reduced' created successfully.") 
        
    
    #Find the index of the filename.
    pattern = r'_imageData_(\d+)\.fits$'
    match = re.search(pattern, im_file)
    file_num = int(match.group(1))
    
    # Open the image file.
    hdul_im = fits.open(directory + im_file)
    hdul_flat = fits.open(directory + flat_file)

    #Remove the overscan regions from the image.
    crop_im = np.delete(np.delete(hdul_im[0].data, np.s_[530:538:1],1), np.s_[0:25:1],1)
    crop_flat = np.delete(np.delete(hdul_flat[0].data, np.s_[530:538:1],1), np.s_[0:25:1],1)

    #Perform the reduction routine: (image - bias - bias offset)/flat
    red_data = crop_im/crop_flat  
    
    
    #Create our basic fits header.
    
    red_header = fits.Header()
    red_header['SIMPLE'] = 'T'
    red_header['BITPIX'] = -64  # 64-bit floating point
    red_header['NAXIS'] = 2
    red_header['NAXIS1'] = crop_im.shape[1]
    red_header['NAXIS2'] = crop_im.shape[0]


    red_hdu = fits.PrimaryHDU(red_data, header=red_header)

    # Write the bias FITS file to 'reduced' directory.
    red_hdu.writeto(directory.replace('images/', 'reduced/')+ im_file.replace('.fits','_reduced.fits'), overwrite=True)
    #print(f"FITS file '{im_file.replace('.fits','_reduced.fits')}' created successfully.")

