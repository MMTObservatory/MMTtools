"""
mmirs_qrp
=========

Python script that quickly reduces longslit and MOS spectra from MMT/MMIRS
"""

import sys, os

from chun_codes import systime

from os.path import exists
from astropy.io import ascii as asc
from astropy.io import fits

from scipy.ndimage.interpolation import shift

import numpy as np

import matplotlib.pyplot as plt
import glob

from astropy.table import Table
from astropy import log

def main(rawdir, prefix, bright=False, silent=False, verbose=True):

    '''
    Main function of mmirs_qrp

    Parameters
    ----------

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 13 October 2017
     - Handle gzip files as inputs
     - Change data cube dimensions
     - Handle FITS overwrite
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    if rawdir[-1] != '/': rawdir = rawdir + '/'

    dcorr_files = glob.glob(rawdir+prefix+'*dcorr.fits*')
    n_files     = len(dcorr_files)
    
    hdu0 = fits.getheader(dcorr_files[0])
    naxis1, naxis2 = hdu0['NAXIS1'], hdu0['NAXIS2']
    
    d_cube0     = np.zeros((n_files, naxis1, naxis2))
    peak_val    = np.zeros(n_files)
    shift_cube0 = np.zeros((n_files, naxis1, naxis2))

    for ii in range(n_files):
        d_cube0[ii] = fits.getdata(dcorr_files[ii])

        if bright == True:
            im_test  = d_cube0[ii].copy()
            med0_row = np.median(im_test, axis=1)

            resize   = np.repeat(med0_row, naxis1).reshape((naxis2,naxis1))
            im_test  = im_test - resize

            med0_col     = np.median(im_test, axis=0)
            peak_val[ii] = np.argmax(med0_col)
        else:
            log.info('### Not supported yet')

    shift_val = peak_val[0] - peak_val

    log.info('### Shift values for spectra : ')
    log.info('### '+" ".join(shift_val))

    for ii in range(n_files):
        shift_cube0[ii] = shift(d_cube0[ii], shift_val[ii])

    fits.writeto(rawdir+prefix+'_stack.fits', shift_cube0, overwrite=True)

    if silent == False: log.info('### End main : '+systime())
#enddef

