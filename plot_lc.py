import numpy as np
from matplotlib import pyplot as plt
import os
from tqdm import trange, tqdm_notebook
from time import sleep

def plot_lc(spool_file, path, n_sources, v, t_occ, delta_t_occ):

    out_dir = path + spool_file.replace('.h5','/photData')
    
    with open(path+spool_file.replace('.h5','/')+str('timingData.txt'), 'r') as file:
        lines = [line.strip() for line in file.readlines()]
        numbers = [float(num) for num in lines]
        timingData = np.array(numbers, dtype=np.float64)

    time_array = timingData-timingData[0]

    
    fig, ax = plt.subplots()
    #plt.title(spool_file.replace('.h5','_lightcurves'))
    plt.xlabel('time (s)')
    plt.ylabel('normalised pixel counts')

    textstr = '\n'.join((
    r'$t_{occ}=%.2f$' % (t_occ, ) +str('$s$'),
    r'$v=%.2f$' % (v, ) +str('$kms^{-1}$'),
    r'$l_{chord}=%.2f$' % (v*t_occ, ) +str('$km$')))

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=1)

    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
    
    for i in range(n_sources):
        with open(out_dir+str('/phot_source_')+str(i+1)+str('.txt'), 'r') as file:
            lines = [line.strip() for line in file.readlines()]
            numbers = [float(num) for num in lines]
            source_array = np.array(numbers, dtype=np.float64)

        norm_source = source_array/np.max(source_array) + i
        
        ax.plot(time_array, norm_source, label=str('source ')+str(i+1), color='blue')
    plt.savefig(out_dir+str('/')+spool_file.replace('.h5','_lightcurves.png'))    
    plt.show()


