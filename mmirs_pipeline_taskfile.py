"""
mmirs_pipeline_taskfile
====

Python code to create task files for MMIRS IDL pipeline:
http://tdc-www.harvard.edu/software/mmirs_pipeline.html

Operates in a given path containing raw files, find common files and create
task files to execute with mmirs_pipeline
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import glob

from astropy.table import Table
from astropy import log

def get_header_info(files0):
    '''
    Get information from FITS header

    Parameters
    ----------
    files0 : list
      List containing full path to MMIRS files

    Returns
    -------
    tab0: astropy.table.table
      Astropy Table containing FITS header

    Notes
    -----
    Created by Chun Ly, 28 November 2017
     - Added documentation
     - Return Astropy Table
    '''

    n_files0 = len(files0)
    filename = np.array(['']*n_files0)
    exptime  = np.zeros(n_files0)
    object0  = np.array(['']*n_files0)
    imagetyp = np.array(['']*n_files0)
    aptype   = np.array(['']*n_files0)
    aperture = np.array(['']*n_files0)
    filter0  = np.array(['']*n_files0)
    disperse = np.array(['']*n_files0)
    
    for ii in range(n_files0):
        hdr = fits.getheader(files0[ii])

        filename[ii] = hdr['FILENAME'].split('/')[-1]
        exptime[ii]  = hdr['EXPTIME']
        object0[ii]  = hdr['OBJECT']
        imagetyp[ii] = hdr['IMAGETYP']
        aptype[ii]   = hdr['APTYPE']
        aperture[ii] = hdr['APERTURE']
        filter0[ii]  = hdr['FILTER']
        disperse[ii] = hdr['DISPERSE']

    #endfor

    arr0   = [filename, exptime, object0, imagetyp, aptype, aperture, filter0, disperse]
    names0 = ('filename','exptime','object','imagetype','aptype','aperture','filter',
              'disperse')
    tab0 = Table(arr0, names=names0)
    return tab0
#enddef

def create(rawdir, silent=False, verbose=True):

    '''
    Main function to create task files to execute

    Parameters
    ----------
    rawdir : str
      Full path to MMIRS files

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 28 November 2017
    '''
    
    if silent == False: log.info('### Begin create : '+systime())

    if rawdir[-1] != '/': rawdir = rawdir + '/'

    files0   = glob.glob(rawdir+'*.????.fits')
    n_files0 = len(files0)
    if silent == False: log.info('### Number of FITS files found : '+str(n_files0))

    
    if silent == False: log.info('### End create : '+systime())
#enddef

