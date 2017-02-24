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
import astropy.units as u

import numpy as np

import matplotlib.pyplot as plt
import glob

from astropy.table import Table

from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder

# Mod on 23/02/2017 to use IRAF's zscale intervals
from astropy.visualization import ZScaleInterval
zscale = ZScaleInterval()
#from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture

from astropy.nddata import Cutout2D

from matplotlib.backends.backend_pdf import PdfPages

out_cat_dir = 'daofind/' # + on 23/02/2017


def get_seqno(files):
    # + on 23/02/2017
    t_files = [os.path.basename(file) for file in files]
    seqno = [file.replace('.fits.gz','').replace('.fits','').split('.')[1] for
             file in t_files]
    return np.array(seqno)

def get_files(path0):
    # + on 23/02/2017
    # Later Mod on 23/02/2017 for seqno
    files = glob.glob(path0+'*fits*') #Mod on 23/02/2017
    seqno = get_seqno(files) # Later + on 23/02/2017

    s_idx = np.argsort(np.array(seqno))

    files = np.array(files)[s_idx]
    return files, seqno[s_idx]
#enddef

def remove_dup_sources(s_cat):
    # + on 23/02/2017
    n_sources = len(s_cat)

    bad = []
    x0 = s_cat['xcentroid']
    y0 = s_cat['ycentroid']

    for nn in xrange(n_sources):
        t_idx  = np.arange(nn+1,n_sources)
        x_diff, y_diff = x0[t_idx]-x0[nn], y0[t_idx]-y0[nn]
        i_match = np.where(np.sqrt(x_diff**2 + y_diff**2) <= 15.0)[0]
        if len(i_match) > 0:
            i_match = t_idx[i_match]
            if len(i_match) == 1:
                if s_cat[i_match]['peak'] < s_cat[nn]['peak']:
                    bad += [i_match.tolist()[0]]
                else:
                    bad += [nn]
            else:
                print nn, 'too many'
                bad += i_match.tolist()
    return bad

def find_stars(files=None, path0=None, plot=False, out_pdf_plot=None,
               silent=False, verbose=True):
    '''
    Find stars in an image and return a catalog of positions

    Parameters
    ----------
    files : list
      List of files

    path0 : string
      Directory path to files.

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 23 February 2017
     - Later modified to plot images and overlay sources
     - Adjust scale to using IRAF's zscale
     - Call remove_dup_sources()
    '''

    if silent == False: log.info('### Begin: '+systime())

    if files == None and path0 == None:
        log.error('files and path0 keywords not provided')
        log.error('Exiting!!!')
        return

    if files == None and path0 != None:
        files, seqno = get_files(path0)
        # path0 = None # Reset since files will have full path
    else:
        if files != None: seqno = get_seqno(files)

    # Later + on 23/02/2017
    out_cat_dir0 = path0+out_cat_dir
    if not exists(out_cat_dir0):
        if silent == False: log.info('Creating : '+out_cat_dir0)
        os.mkdir(out_cat_dir0)

    s_date = path0.split('/')[-1]

    if plot == True:
        if out_pdf_plot == None:
            out_pdf_plot = path0+'find_stars.pdf'
        pp = PdfPages(out_pdf_plot)

    for ff in xrange(len(files)):
        basename = os.path.basename(files[ff])
        image = fits.getdata(files[ff])
        mean, median, std = sigma_clipped_stats(image, sigma=2.0, iters=5)
        image_sub = image - median

        if verbose == True:
            log.info('%s mean/med/sig: %f %f %f' %
                     (seqno[ff], mean, median, std))
        #Later Mod on 23/02/2017 to lower threshold
        daofind = DAOStarFinder(fwhm=8.0, threshold=5.*std)
        s_cat = daofind(image_sub)

        # Exclude saturated objects
        unsat = np.where(s_cat['peak'] <= 60000.0)[0]
        sat   = np.where(s_cat['peak'] >  60000.0)[0]
        cat_sat = s_cat[sat]
        s_cat = s_cat[unsat]
        s_cat.sort(['peak'])
        s_cat.reverse()
        # s_cat.pprint()

        # Later + on 23/02/2017
        bad = remove_dup_sources(s_cat)
        print bad
        if len(bad) > 0:
            cat_bad = s_cat[bad]
            s_cat.remove_rows(bad)

        if ff == 0 and silent == False: s_cat.pprint()

        # + on 23/02/2017
        out_cat = out_cat_dir0+basename.replace('.fits.gz','.tbl').replace('.fits','.tbl')
        s_cat.write(out_cat, format='ascii.fixed_width_two_line',
                    overwrite=True)

        # Later + on 23/02/2017
        if len(bad) >0 and verbose == True:
            log.info('The following will be removed : ')
            cat_bad.pprint()
            out_cat_bad = out_cat.replace('.tbl','.bad.tbl')
            cat_bad.write(out_cat_bad, format='ascii.fixed_width_two_line',
                          overwrite=True)

        if plot == True:
            pos0  = (s_cat['xcentroid'], s_cat['ycentroid'])
            aper0 = CircularAperture(pos0, r=8.)

            pos0      = (cat_sat['xcentroid'], cat_sat['ycentroid'])
            sat_aper0 = CircularAperture(pos0, r=8.)

            pos0      = (cat_bad['xcentroid'], cat_bad['ycentroid'])
            bad_aper0 = CircularAperture(pos0, r=8.)

            fig, ax = plt.subplots()
            z1, z2 = zscale.get_limits(image_sub)
            print z1, z2
            norm = ImageNormalize(vmin=z1, vmax=z2) #stretch=SqrtStretch())
            ax.imshow(image_sub, cmap='Greys', origin='lower', norm=norm)
            aper0.plot(color='blue', lw=1.5, alpha=0.5)
            sat_aper0.plot(color='red', lw=1.5, alpha=0.5)
            bad_aper0.plot(color='magenta', lw=1.5, alpha=0.5)

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

def make_postage(files=None, path0=None, n_stack=5, size=50,
                 silent=False, verbose=True):
    '''
    Create cut-outs and median stack to produce image of the
    point-spread function

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
        # path0 = None # Reset since files will have full path
    else:
        if files != None: seqno = get_seqno(files)

    post_dir0 = path0 + 'post/'
    if not exists(post_dir0):
        if silent == False: log.info('Creating : '+post_dir0)
        os.mkdir(post_dir0)

    out_cat_dir0 = path0+out_cat_dir
    for ff in xrange(len(files)):
        basename = os.path.basename(files[ff])
        image = fits.getdata(files[ff])
        mean, median, std = sigma_clipped_stats(image, sigma=2.0, iters=5)
        image_sub = image - median

        in_cat = out_cat_dir0+basename.replace('.fits.gz','.tbl').\
                 replace('.fits','.tbl')
        s_cat  = asc.read(in_cat, format='fixed_width_two_line')

        n_bright = np.min([n_stack,len(s_cat)])
        bright   = range(n_bright)
        s_cat    = s_cat[bright]

        x0 = np.round_(s_cat['xcentroid'])
        y0 = np.round_(s_cat['ycentroid'])

        im0 = np.zeros( (len(bright), size, size))
        size2d = u.Quantity((size, size), u.pixel)
        for ii in range(n_bright):
            pos0 = (x0[ii], y0[ii])
            cutout = Cutout2D(image_sub, pos0, size2d).data
            im0[ii] = cutout/np.max(cutout)

        out_fits = post_dir0+seqno[ff]+'.fits'
        psf_im = np.median(im0, axis=0)
        fits.writeto(out_fits, psf_im, overwrite=True)

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
        # path0 = None # Reset since files will have full path

    if silent == False:
        log.info('The following files will be analyzed from : ')
        log.info(path0)
        for file in files: log.info(os.path.basename(file))

    find_stars(files=files, path0=path0, plot=True, verbose=False)
    make_postage(files=files, path0=path0, verbose=False)

    if silent == False: log.info('### End: '+systime())
