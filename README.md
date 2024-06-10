
README FOR OCC_FIND SOFTWARE PACKAGE


Author: Ben Attwood

occ_find is a Python-based software package for finding occultation events from HDF5 data files from TCI on the 1.54m Danish telescope.

Before running, ensure the following python packages are installed:

Python 3.12.1
numpy 1.2.24
astropy 5.2.2
photutils 1.8.0
matplotlib 3.7.4
pandas 2.0.3
imageio 2.33.1
PyTables 3.8.0
tqdm 4.66.1

The core package consists of 2 scripts:

1. h5_unpack.py  - For extracting image data and metadata from HDF5 file to a directory.
3. occ_find.py - For performing automated, high speed aperture photometry to find stellar occultation events.

There is an additional 3 scripts in the package which are not essential:

1. reducto.py (optional, see below) - For reduction of images by flat-fielding.
2. occ_timer - For timing duration of any detected occultations, this script is rather buggy and not currently reccommended for use.
3. plot_lc - For re-plotting and re-formatting any obtained lightcurves.

Further information on each script is given as a commented block at the top of each source script.


HOW TO USE OCC_FIND

occ_find is optimised for use in the jupyter notebook platform, and an example notebook is included called 'occ_find tutorial.ipynb'. Please study this example notebook before running on any spool files. It is reccommended that a directory 'occ_find' (or whatever you want to call it) is created, and all code is stored in there, along with any spool files to be processed.
