####################################################################################################################################################
#########################################################  OCCULTATION FINDING PROGRAM  ############################################################
####################################################################################################################################################

''' PYTHON VERSION: 3.12.1'''
''' Code to go through a folder of .fits files, detect the sources in the images, perform aperture photometry on each source in every image and return a lightcurve as a plot and .txt files with lightcurve data for each source.''' 

####################################################################################################################################################
###############################################################  INPUTS & OUTPUTS  #################################################################
####################################################################################################################################################

''' INPUTS:
1. 'path' - str - directory in which the extracted spool data is located.
2. 'spool_file' - str - spool file from which the data originated. The spool file will be read in to determine the number of images in the spool                            file.
3. 'make_dir' - bool - Creates a 'PhotData' directory if marked 'True'.
4. 'reduced' - bool - Specifies whether or not the images have been reduced.
5. 'fwhm' - int - Full-Width at Half Maximum to specfy for PSF fitting by IRAFStarFinder, refer to photutils documentation for further details.
6. 'thresh_factor' - float - Specifies the percentage above the median pixel value at which IRAFStarFinder will look for sources.

IMPORTANT NOTE ON INPUT FILES: IF NO TIMESTAMP DATA IS SPECIFIED IN THE FITS FILES, THE CODE WILL PLOT/WRITE TO FILE, THE IMAGE NUMBER INSTEAD OF THE TIME.

    OUTPUTS:

Directory called photData containing:
    1. phot_source_i.txt - txt file containg background subtracted photometry data for source 'i'
    2. source_list.txt - txt file containing the list of sources detected by IRAFStarFinder on which, photometry has been performed. This is               printed when occ_find is run.
    3. spool_file_lightcurve_source_i.png - png file showing the lightcurve plot for source 'i'. These are displayed when occ_find is run.
    4. sources_plot.png - Plot of photometry apertures overlaid onto detected sources.
'''

#######################################################################################################################################################################################################################  IMPORTS  ######################################################################
####################################################################################################################################################

import numpy as np
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy.stats import sigma_clipped_stats   
from astropy import visualization as aviz
from photutils.aperture import CircularAnnulus, CircularAperture, ApertureStats, aperture_photometry
from photutils.detection import IRAFStarFinder

import tables as tb
import os
from tqdm import trange, tqdm_notebook
from time import sleep

######################################################################################################################################################################################################################   CODE   ########################################################################
####################################################################################################################################################


def occ_find(path, spool_file, make_dir, reduced, fwhm, thresh_factor):

    out_dir = path + spool_file.replace('.h5','/photData')
    
    if make_dir == True:
        os.mkdir(out_dir)
        print(f"Directory 'photData' created successfully.") 
        
    t = tb.open_file(spool_file)  #Open spool file with PyTables, this is for easy reading of the lengths of datasets.

    #Specifying different paths to images depending on whether or not they are reduced.
    if reduced == True:
        hdul_start = fits.open(path+spool_file.replace('.h5','/reduced/')+spool_file.replace('.h5','_imageData_0_reduced.fits'))
    else:
        hdul_start = fits.open(path+spool_file.replace('.h5','/images/')+spool_file.replace('.h5','_imageData_0.fits'))
        
    hdul_start.info()
    mean = np.mean(hdul_start[0].data)
    med = np.median(hdul_start[0].data)
    std = np.std(hdul_start[0].data)
    print(mean, med, std)                                   
    threshold = thresh_factor*med
    
    #Source finding function.
    find = IRAFStarFinder(fwhm, threshold, sharplo=0, sharphi=1, roundlo = -1, roundhi=1)

    sources = find(hdul_start[0].data - mean)  
    for col in sources.colnames:  
        if col not in ('id', 'npix'):
            sources[col].info.format = '%.2f'  # for consistent table output

    source_doc = open(out_dir+str('/source_list.txt'), "a")
    source_doc.write(str(sources))
    source_doc.close()
    sources.pprint(max_width=100)
    
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=fwhm)
    annulus_aperture = CircularAnnulus(positions, r_in=2*fwhm, r_out=3*fwhm)

    norm = aviz.ImageNormalize(hdul_start[0].data, interval=aviz.ZScaleInterval())
    plt.imshow(hdul_start[0].data,origin='lower', norm=norm, interpolation='nearest')
    plt.colorbar(label='pixel value')
    apertures.plot(color='red', lw=1.5)
    annulus_aperture.plot(color='cyan',lw=1.5)
    plt.savefig(out_dir+str('/sources_plot.png'))
    plt.show()

    # Looping through all images to peform aperture photometry at the positions of the sources found by IRAF Star Finder in the first image of the        spool. After performing aperture photometry on each image, the photometry results are written to a .txt file for each source.
    bkg_counts = []

    if reduced == True:
        for i in trange(len(t.root.imageData)):
            hdul_im = fits.open(path+spool_file.replace('.h5','/reduced/')+spool_file.replace('.h5','_imageData_')+str(i)+str('_reduced.fits'))
            aperstats = ApertureStats(hdul_im[0].data, annulus_aperture)
            bkg_mean = aperstats.mean
            phot_table = aperture_photometry(hdul_im[0].data, apertures)
    
            apertures.area 
            aperture_area = apertures.area_overlap(hdul_im[0].data)
    
            total_bkg = bkg_mean * aperture_area
    
            phot_bkgsub = phot_table['aperture_sum'] - total_bkg
            phot_table['total_bkg'] = total_bkg
            phot_table['aperture_sum_bkgsub'] = phot_bkgsub
    
            for col in phot_table.colnames:
                phot_table[col].info.format = '%.8g'  # for consistent table output
    
            for j in range(len(sources)):
                o = open(out_dir + str('/phot_source_') + str(j+1) + str('.txt'), "a")
                if phot_table[j][5] < 0:
                    o.write(str(0)+'\n')
                else:
                    o.write(str(phot_table[j][5])+'\n')
                o.close()
            sleep(0.01)
    else:
        
        for i in trange(len(t.root.imageData)):
            hdul_im = fits.open(path+spool_file.replace('.h5','/images/')+spool_file.replace('.h5','_imageData_')+str(i)+str('.fits'))
            aperstats = ApertureStats(hdul_im[0].data, annulus_aperture)
            bkg_mean = aperstats.mean
            phot_table = aperture_photometry(hdul_im[0].data, apertures)
    
            apertures.area 
            aperture_area = apertures.area_overlap(hdul_im[0].data)
    
            total_bkg = bkg_mean * aperture_area
    
            phot_bkgsub = phot_table['aperture_sum'] - total_bkg
            phot_table['total_bkg'] = total_bkg
            phot_table['aperture_sum_bkgsub'] = phot_bkgsub
    
            for col in phot_table.colnames:
                phot_table[col].info.format = '%.8g'  # for consistent table output
    
            for j in range(len(sources)):
                o = open(out_dir + str('/phot_source_') + str(j+1) + str('.txt'), "a")
                if phot_table[j][5] < 0:
                    o.write(str(0)+'\n')
                else:
                    o.write(str(phot_table[j][5])+'\n')
                o.close()
            sleep(0.01)

    with open(path+spool_file.replace('.h5','/')+str('timingData.txt'), 'r') as file:
        lines = [line.strip() for line in file.readlines()]
        numbers = [float(num) for num in lines]
        timingData = np.array(numbers, dtype=np.float64)

    time_array = timingData-timingData[0]

    
    
    plt.title(spool_file.replace('.h5','_lightcurves'))
    plt.xlabel('time (s)')
    plt.ylabel('normalised pixel counts')
    
    for i in range(len(sources)):
        with open(out_dir+str('/phot_source_')+str(i+1)+str('.txt'), 'r') as file:
            lines = [line.strip() for line in file.readlines()]
            numbers = [float(num) for num in lines]
            source_array = np.array(numbers, dtype=np.float64)

        norm_source = source_array/np.max(source_array) + i
        
        plt.plot(time_array, norm_source, label=str('source ')+str(i+1))
    plt.legend(loc=1)
    plt.savefig(out_dir+str('/')+spool_file.replace('.h5','_lightcurve.png'))    
    plt.show()    
