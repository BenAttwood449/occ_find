from astropy.io import fits
from astropy import visualization as aviz
import numpy as np
import matplotlib.pyplot as plt

def fits_to_png(fits_filename, im_dir, save_dir):
    # Open the FITS file
    hdul = fits.open(im_dir + fits_filename)
    # Remove the overscan region.
    crop_im = np.delete(np.delete(hdul[0].data, np.s_[530:538:1],1), np.s_[0:25:1],1)
    norm = aviz.ImageNormalize(hdul[0].data, interval=aviz.ZScaleInterval())
    # Display the image using matplotlib
    plt.imshow(crop_im, origin='lower', norm=norm)
    plt.colorbar()
    plt.title(fits_filename.replace('.fits',''))
    # Save the image as PNG
    plt.savefig(save_dir + fits_filename.replace('.fits','.png'))
    plt.close()
    hdul.close()
    #print(f"Saved {fits_filename.replace('.fits','.png')}")