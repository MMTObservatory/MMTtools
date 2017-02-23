"""
mmtcam
======

Provide description for code here.
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits
from astropy import log

import numpy as np

import matplotlib.pyplot as plt
import glob

from astropy.table import Table

from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder

from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture
norm = ImageNormalize(stretch=SqrtStretch())

from matplotlib.backends.backend_pdf import PdfPages

def get_seqno(files):
    t_files = [os.path.basename(file) for file in files]
    seqno = [file.replace('.fits.gz','').replace('.fits','').split('.')[1] for
             file in t_files]
    return np.array(seqno)

def get_files(path0):
    # Later Mod on 23/02/2017
    files = glob.glob(path0+'/*fits*')
    seqno = get_seqno(files) # Later + on 23/02/2017

    s_idx = np.argsort(np.array(seqno))

    files = np.array(files)[s_idx]
    return files, seqno[s_idx]
#enddef

def find_stars(files=None, path0=None, plot=False, out_pdf_plot=None,
               silent=False, verbose=True):
    '''
    Find stars in an image and return a catalog of positions

    Parameters
    ----------
    files : list
      List of files

    path0 : string
      Path to files. If not provided it is assumed that [files] has the full
      path name

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 23 February 2017
    '''

    if silent == False: log.info('### Begin: '+systime())

    if files == None and path0 == None:
        log.error('files and path0 keywords not provided')
        log.error('Exiting!!!')
        return

    if files == None and path0 != None:
        files, seqno = get_files(path0)
        path0 = None # Reset since files will have full path
    else:
        seqno = get_seqno(files)

    if path0 == None:
        s_date = os.path.dirname(files[0]).split('/')[-1]

    if plot == True:
        if out_pdf_plot == None:
            out_pdf_plot = os.path.dirname(files[0])+'/find_stars.pdf'
        pp = PdfPages(out_pdf_plot)

    for ff in xrange(len(files)):
        image = fits.getdata(files[ff])
        mean, median, std = sigma_clipped_stats(image, sigma=3.0, iters=5)    
        if verbose == True:
            log.info('%s : %f %f %f' % (seqno[ff], mean, median, std))
        daofind = DAOStarFinder(fwhm=8.0, threshold=10.*std)    
        s_cat = daofind(image - median)

        # Exclude saturated objects
        unsat = np.where(s_cat['peak'] <= 60000.0)[0]
        sat   = np.where(s_cat['peak'] > 60000.0)[0]
        cat_sat = s_cat[sat]
        s_cat = s_cat[unsat]
        s_cat.sort(['peak'])
        s_cat.reverse()
        if ff == 0: s_cat.pprint()

        if plot == True:
            pos0  = (s_cat['xcentroid'], s_cat['ycentroid'])
            aper0 = CircularAperture(pos0, r=8.)

            pos0      = (cat_sat['xcentroid'], cat_sat['ycentroid'])
            sat_aper0 = CircularAperture(pos0, r=8.)

            fig, ax = plt.subplots()
            ax.imshow(image-median, cmap='Greys', origin='lower', norm=norm)
            aper0.plot(color='blue', lw=1.5, alpha=0.5)
            sat_aper0.plot(color='red', lw=1.5, alpha=0.5)

            t_ann = s_date+'/'+os.path.basename(files[ff])
            ax.annotate(t_ann, [0.025,0.975], xycoords='axes fraction',
                         ha='left', va='top')
            fig.savefig(pp, format='pdf', bbox_inches='tight')

    #endfor
    if plot == True:
        if silent == False:
            log.info('## Writing : '+out_pdf_plot+' | '+systime())
        pp.close()

    if silent == False: log.info('### End: '+systime())
#enddef

def run_all(files=None, path0=None, silent=False, verbose=True):
    '''
    Run all functions related to MMTCam analysis

    Parameters
    ----------
    files : list
      List of files

    path0 : string
      Path to files. If not provided it is assumed that [files] has the full
      path name

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 23 February 2017
    '''

    if files == None and path0 == None:
        log.error('files and path0 keywords not provided')
        log.error('Exiting!!!')
        return

    if silent == False: log.info('### Begin: '+systime())

    if files == None and path0 != None:
        files, seqno = get_files(path0)
        path0 = None # Reset since files will have full path

    if silent == False:
        log.info('The following files will be analyzed from : ')
        log.info(path0)
        for file in files: log.info(os.path.basename(file))

    find_stars(files=files, path0=path0, plot=True)

    if silent == False: log.info('### End: '+systime())
