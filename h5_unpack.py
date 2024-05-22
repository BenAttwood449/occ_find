###################################################################################################################################################################################################################  HDF5 FILE UNPACKER  ###############################################################
####################################################################################################################################################

''' PYTHON VERSION: 3.12.1'''
''' Code to open a hdf5 spool file from the 1.54m Danish telescope and convert to a directory containing the data as .fits and .txt files.'''

####################################################################################################################################################
#############################################################  INPUTS & OUTPUTS  ###################################################################
####################################################################################################################################################

'''INPUTS: 

1. file_in - Input hdf5 file.
2. out_dir - Directory for output.


   
   OUTPUTS:
   
Directory containing:
1. 'flat' - flat frame image, .fits file
2. 'bias' - bias frame image, .fits file
3. 'imageData' - directory containing raw, unreduced .fits files of observations
4. 'cosmicsData' - data on CCD cosmic ray strikes, .txt file.
5. 'reducedData' - corrections required to reduce the 'imageData' files, .txt file.
6. 'timingData' - Timestamp data for images, .txt file.

For further information on the 'cosmicsData' and 'reducedData' datasets, or the structure of the input spool files, see Skotfelt et.al., A&A 574, A54, 2015 for further details.
'''
######################################################################################################################################################################################################################  IMPORTS  #######################################################################
####################################################################################################################################################

import numpy as np
import os
import tables as tb
from astropy.io import fits
from datetime import datetime
from tqdm import trange, tqdm_notebook
from time import sleep

######################################################################################################################################################################################################################   CODE   ########################################################################
####################################################################################################################################################


def h5_unpack(file, out_dir):


    # Make directory for files.
    filename = file.replace('.h5', '')
    parent_dir = out_dir
    im_dir = 'images'
    file_path = os.path.join(parent_dir, filename)
    im_path = os.path.join(parent_dir, filename, im_dir)
    os.mkdir(file_path)
    os.mkdir(im_path)
    
    print("Directory '% s' created" % filename)     
    print("Directory '% s' created" % im_dir) 

    
    t = tb.open_file(file)
    tabs = ['bias','flat','imageData','reducedData','cosmicsData','spurious','timingData']



    
    # Create the bias header
    bias_header = fits.Header()
    bias_header['SIMPLE'] = 'T'
    bias_header['BITPIX'] = -64  # 64-bit floating point
    bias_header['NAXIS'] = 2
    bias_header['NAXIS1'] = t.root.bias.shape[1]
    bias_header['NAXIS2'] = t.root.bias.shape[0]

    # Create a bias HDU (Header Data Unit)
    bias_hdu = fits.PrimaryHDU(t.root.bias[:,:], header=bias_header)

    # Write the bias FITS file
    bias_hdu.writeto(out_dir + filename + str('/') + im_dir + str('/') + filename + str('_bias.fits'), overwrite=True)
    print(f"FITS file '{filename + str('_bias.fits')}' created successfully.")

            
    # Create the flat header
    flat_header = fits.Header()
    flat_header['SIMPLE'] = 'T'
    flat_header['BITPIX'] = -64  # 64-bit floating point
    flat_header['NAXIS'] = 2
    flat_header['NAXIS1'] = t.root.flat.shape[1]
    flat_header['NAXIS2'] = t.root.flat.shape[0]

    # Create the flat HDU
    flat_hdu = fits.PrimaryHDU(t.root.flat[:,:], header=flat_header)

    # Write the flat FITS file
    flat_hdu.writeto(out_dir +filename + str('/') + im_dir + str('/') + filename + str('_flat.fits'), overwrite=True)
    print(f"FITS file '{filename + str('_flat.fits')}' created successfully.")



    
    times = []
    for i in range(len(t.root.imageData)):
        (d,m) = divmod(i,100)
        _t = t.root.timingData[d][0] + t.root.imageData.attrs['ACC_TIME']*m
        times.append(_t)
    times_utc = [datetime.utcfromtimestamp(_t).isoformat() for _t in times]
    

    a = open(out_dir + str('/') + filename + str('/') + str('cosmicsData.txt'), 'w')
    b = open(out_dir + str('/') + filename + str('/') + str('reducedData.txt'), 'w')
    c = open(out_dir + str('/') + filename + str('/') + str('timingData.txt'), 'w')
    d = open(out_dir + str('/') + filename + str('/') + str('timingData_UTC.txt'), 'w')
   # e = open(out_dir + file_dir + str('imageData_header.txt'), 'w')

    with open(out_dir + str('/') + filename + str('/') + str('cosmicsData.txt'), 'w') as a:
        for i in range(len(t.root.cosmicsData)):
            a.write(str(t.root.cosmicsData[i])+'\n')
    print(f"txt file 'cosmicsData.txt' created successfully.")
    
    with open(out_dir + str('/') + filename + str('/') + str('reducedData.txt'), 'w') as b:
        for i in range(len(t.root.reducedData)):
            b.write(str(t.root.reducedData[i])+'\n')
    print(f"txt file 'reducedData.txt' created successfully.")
    
    with open(out_dir + str('/') + filename + str('/') + str('timingData.txt'), 'w') as c:
        for i in range(len(times)):
            c.write(str(times[i])+'\n')
    print(f"txt file 'timingData.txt' created successfully.")
    
    with open(out_dir + str('/') + filename + str('/') + str('timingData_UTC.txt'), 'w') as d:
        for i in range(len(times_utc)):
            d.write(str(times_utc[i])+'\n')
    print(f"txt file 'timingData_UTC.txt' created successfully.")

    with open(out_dir + str('/') + filename + str('/') + str('headerData.txt'), 'w') as e:
        e.write(str('ACC_TIME: ')+ str(t.root.imageData.attrs.configuration['ACC_TIME'])+'\n')
        e.write(str('AIRMASS: ')+ str(t.root.imageData.attrs.configuration['AIRMASS'])+'\n')
        e.write(str('DATE_OBS: ')+ str(t.root.imageData.attrs.configuration['DATE_OBS'])+'\n')
        e.write(str('OBJECT: ')+ str(t.root.imageData.attrs.configuration['OBJECT'])+'\n')
        e.write(str('TEL_RA: ')+ str(t.root.imageData.attrs.configuration['TEL_RA'])+'\n')
        e.write(str('TEL_DEC: ')+ str(t.root.imageData.attrs.configuration['TEL_DEC'])+'\n')
        e.write(str('EXP_TIME: ')+ str(t.root.imageData.attrs.configuration['EXP_TIME'])+'\n')
        e.write(str('EM_GAIN: ')+ str(t.root.imageData.attrs.configuration['EM_GAIN'])+'\n')
    print(f"txt file 'headerData.txt' created successfully.")


    for k in trange(len(t.root.imageData)):
           
        # Create a FITS header
        im_header = fits.Header()
        im_header['SIMPLE'] = 'T'
        im_header['BITPIX'] = -64                                                     # 64-bit floating point
        im_header['NAXIS'] = 2
        im_header['NAXIS1'] = t.root.imageData[k][0].shape[1]
        im_header['NAXIS2'] = t.root.imageData[k][0].shape[0]

        # Create a FITS HDU (Header Data Unit)
        im_hdu = fits.PrimaryHDU(t.root.imageData[k][0], header=im_header)

        # Write the FITS file
        im_hdu.writeto(out_dir + filename + str('/') + im_dir + str('/') + filename + str('_imageData_') + str(k) +str('.fits'), overwrite=True)
        sleep(0.01)


