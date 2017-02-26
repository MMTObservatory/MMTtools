"""
mmtcam
======

Set of functions to identify stars and stack them to construct the PSF for
MMTCam for each image. Intended to understand the cause of unusual PSFs
(e.g., double, oscillations, elongated profiles)
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
from pylab import subplots_adjust # + on 24/02/2017

import scipy.optimize as opt # + on 25/02/2017

# + on 26/02/2017
from datetime import datetime
from astropy.time import Time, TimezoneInfo

out_cat_dir = 'daofind/' # + on 23/02/2017

# + on 24/02/2017
bbox_props = dict(boxstyle="square,pad=0.3", fc="white", alpha=0.75, ec="none")

c_levels = 0.2+0.1*np.arange(9)

f_s = 2*np.sqrt(2*np.log(2)) # sigma-FWHM conversion | + on 25/02/2017
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
#enddef

def hdr_annotate(h0, ax):
    '''
    Define string for annotation in psf_contours()

    Parameters
    ----------
    h0 : astropy.io.fits.header.Header
      FITS header file

    ax : matplotlib.axes._subplots.AxesSubplot
      Axes to use for annotation

    Returns
    -------
    None.

    Notes
    -----
    Created by Chun Ly, 24 February 2017
    Modified by Chun Ly, 26 February 2017
     - Add commmanded RA,Dec offsets to annotated text
    '''

    txt0 = r'$t_{\rm exp}$=%.1f, sz=%.2f, HA=%s' % (h0['EXPTIME'],h0['AIRMASS'],
                                                    h0['HA'])
    txt0 += '\n'

    txt0 += r'$\alpha$=%s, $\delta$=%s' % (h0['RA'], h0['DEC'])
    txt0 += '\n'

    # Include commanded offsets | + on 26/02/2017
    rao, deo = np.float(h0['RAOFF']), np.float(h0['DECOFF'])
    txt0 += r'Offsets: $\alpha$=%.2f", $\delta$=%.2f"' % (rao, deo)
    txt0 += '\n'

    txt0 += 'Alt=%s, Az=%s' % (h0['OBJCTALT'], h0['OBJCTAZ'])
    txt0 += '\n'

    txt0 += r'$\theta_{\rm rot}$=%.2f, $\theta_{\rm para}$=%.2f, ' % \
            (h0['ROTANGLE'], h0['PARANG'])
    txt0 += r'$\theta_{\rm pos}$=%.2f' % h0['POSANG']
    ax.annotate(txt0, [0.025,0.025], ha='left', va='bottom',
                xycoords='axes fraction', fontsize=8, bbox=bbox_props)

def get_mst(h0):
    '''
    Get string-formatted MST time from UTC datetime

    Parameters
    ----------
    h0 : astropy.io.fits.header.Header
      FITS header file

    Returns
    -------
    t.strftime : string
      Time in HH:MM:SS format

    Notes
    -----
    Created by Chun Ly, 26 February 2017
    '''

    utc_mst = TimezoneInfo(utc_offset=-7*u.hour)
    t = Time(h0['DATE-OBS']).to_datetime(timezone=utc_mst)
    return t.strftime('%H:%M:%S')
#enddef

def fwhm_fwqm_size(post, pscale):
    '''
    Computes FWHM based on a provided cutout of the PSF. Code uses the number
    of pixels. Handles non-Gaussian PSFs

    Parameters
    ----------
    image : numpy.array
      2-D image

    Returns
    -------
    fwhm0 : float
      Full-width at half maximum based on number of pixels (area)

    fwqm0 : float
      Full-width at quarter maximum based on number of pixels (area)

    Notes
    -----
    Created by Chun Ly, 25 February 2017
    '''

    max0 = np.max(post)

    i_fwhm    = np.where(post >= max0/2.0)
    i_fwqm    = np.where(post >= max0/4.0)

    # Area in arcsec^2
    area_fwhm = len(i_fwhm[0]) * (pscale.to(u.arcsec).value)**2
    area_fwqm = len(i_fwqm[0]) * (pscale.to(u.arcsec).value)**2

    fwhm0     = 2.0*np.sqrt(area_fwhm/np.pi)
    fwqm0     = 2.0*np.sqrt(area_fwqm/np.pi)

    return fwhm0, fwqm0
#enddef

def gauss2d((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    '''
    2-D Gaussian for opt.curve_fit()

    Parameters
    ----------
    (x,y) : numpy.ndarray
      x,y grid from numpy.meshgrid()

    amplitude : float
      Peak of Gaussian

    xo : float
      Gaussian center value along x

    yo : float
      Gaussian center value along y

    sigma_x : float
      Gaussian sigma along x

    sigma_y : float
      Gaussian sigma along y

    theta : float
      Orientation along major axis of Gaussian. Positive is clock-wise.

    offset : float
      Level of continuum

    Returns
    -------
    g.ravel() : numpy.ndarray
      Contiguous flattened array

    Notes
    -----
    Created by Chun Ly, 25 February 2017
    '''

    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g.ravel()

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

    s_date = path0.split('/')[-2] # Mod on 25/02/2017 for minor bug

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
            fig.savefig(pp, format='pdf') #, bbox_inches='tight')

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
    Modified by Chun Ly, 24 February 2017
     - Include FITS header in cutout images
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
        image, hdr = fits.getdata(files[ff], header=True)
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
            cutout = Cutout2D(image_sub, pos0, size2d, mode='partial',
                              fill_value=np.nan)
            im0[ii] = cutout.data/np.max(cutout.data)

        out_fits = post_dir0+seqno[ff]+'.fits'
        psf_im = np.nanmedian(im0, axis=0)
        fits.writeto(out_fits, psf_im, hdr, overwrite=True)

    if silent == False: log.info('### End: '+systime())
#enddef

def psf_contours(files=None, path0=None, out_pdf_plot=None, silent=False,
                 verbose=True):
    '''
    Generate contour plots for the MMTCam PSF

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
    Created by Chun Ly, 24 February 2017
     - Later mod to handle plotting styles
     - Later Mod to include header info in annotation
     - Use filled contours with plasma cmap
     - Add colorbar
    Modified by Chun Ly, 25 February 2017
     - Add colorbar for last subplot
     - Get FWHM and FWQM from fwhm_fwqm_image()
     - Call opt.curve_fit() to fit 2-D Gaussians
     - Overlay cyan contours for best 2-D fit
    '''

    if files == None and path0 == None:
        log.error('files and path0 keywords not provided')
        log.error('Exiting!!!')
        return

    if silent == False: log.info('### Begin psf_contours: '+systime())

    if files == None and path0 != None:
        files, seqno = get_files(path0)
    else:
        if files != None: seqno = get_seqno(files)

    post_dir0 = path0 + 'post/'

    if out_pdf_plot == None: out_pdf_plot = path0+'psf_contours.pdf'

    pp = PdfPages(out_pdf_plot)

    pscale = 0.16 * u.arcsec
    v_pscale = pscale.to(u.arcsec).value
    ncols, nrows = 3, 3

    n_files = len(files)
    for ff in xrange(n_files):
        psf_file = post_dir0+seqno[ff]+'.fits'
        psf_im, h0 = fits.getdata(psf_file, header=True)
        psf_im  /= np.max(psf_im)

        if ff == 0:
            shape0 = psf_im.shape
            x0 = pscale*np.arange(-1*shape0[0]/2.0,shape0[0]/2.0)
            y0 = pscale*np.arange(-1*shape0[1]/2.0,shape0[1]/2.0)

        if ff % (ncols*nrows) == 0:
            fig, ax = plt.subplots(nrows, ncols)

        row, col = ff / ncols % nrows, ff % ncols

        # Later mod on 24/02/2017
        cf = ax[row,col].contourf(x0,y0,psf_im, levels=c_levels,
                                  cmap=plt.cm.plasma)

        # Mod on 25/02/2017 to include colorbar for last subplot
        # Mod on 26/02/2017 to shrink height
        if col == ncols-1:
            cax = fig.add_axes([0.925, 0.76-0.32*row, 0.01, 0.14])
        if ff == n_files-1:
            cax = fig.add_axes([0.605, 0.76-0.32*row, 0.01, 0.14])
        if col == ncols-1 or ff == n_files-1:
            cbar = fig.colorbar(cf, ax=ax[row,col], cax=cax)
            cbar.ax.tick_params(labelsize=8)

        if row == nrows-1:
            ax[row,col].set_xlabel('X [arcsec]')
        else:
            if ((n_files-1)-ff) > ncols-1:
                ax[row,col].set_xticklabels([])

        if ff == n_files-1:
            for cc in range(ncols): ax[row,cc].set_xlabel('X [arcsec]')

        if col == 0:
            ax[row,col].set_ylabel('Y [arcsec]')
        else: ax[row,col].set_yticklabels([])

        # Mod on 26/02/2017
        t_label = seqno[ff]+'.'+h0['FILTER']
        ax[row,col].annotate(t_label, [0.025,0.975], weight='bold', ha='left',
                             va='top', xycoords='axes fraction', fontsize=10)

        # Compute image quality | + on 25/02/2017
        f_annot = 'UTC='+h0['UT']+'  MST='+get_mst(h0)+'\n' # + on 26/02/2017
        fwhm0, fwqm0 = fwhm_fwqm_size(psf_im, pscale)
        f_annot += 'Area: FWHM=%.2f", FWQM=%.2f"\n' % (fwhm0, fwqm0)

        sigG = fwhm0/f_s/pscale.to(u.arcsec).value
        ini_guess = (1.0, 25, 25, sigG, sigG, 0.0, 0.0)
        gx = np.linspace(0,shape0[0]-1,shape0[0])
        gy = np.linspace(0,shape0[1]-1,shape0[1])
        gx, gy = np.meshgrid(gx, gy)

        psf_im_re = psf_im.reshape(shape0[0]*shape0[1])
        popt, pcov = opt.curve_fit(gauss2d, (gx, gy), psf_im_re, p0=ini_guess)
        FWHMx = popt[3] * f_s * pscale.to(u.arcsec).value
        FWHMy = popt[4] * f_s * pscale.to(u.arcsec).value
        f_annot += r'2DFit: FW$_1$=%.2f", FW$_2$=%.2f", ' % (FWHMx,FWHMy)
        f_annot += r'$\theta$=%.2f' % np.degrees(popt[5])
        ax[row,col].annotate(f_annot, [0.025,0.915], xycoords='axes fraction',
                             ha='left', va='top', fontsize=8, zorder=10)

        # Overlay 0.25, 0.50, and 0.75 quartile contours for 2-D Gaussian fit
        # + on 25/02/2017
        f_data = gauss2d((gx, gy), *popt)
        levels = np.array([0.25,0.5,0.75])*np.max(f_data)
        ax[row,col].contour(x0, y0, f_data.reshape(shape0[0],shape0[1]),
                            colors='c', linewidth=2, cmap=None,
                            levels=levels.tolist())
        # Plot center of fit | + on 26/02/2017
        xcen = v_pscale*(-1*shape0[0]/2.0+popt[1])
        ycen = v_pscale*(-1*shape0[1]/2.0+popt[2])
        ax[row,col].plot(xcen, ycen, 'o', mfc='c', mec='none', alpha=0.5)

        hdr_annotate(h0, ax[row,col]) # + on 24/02/2017

        if ff == n_files-1:
            for cc in range(col+1,ncols): ax[row,cc].axis('off')
            for rr in range(row+1,nrows):
                for cc in range(ncols): ax[rr,cc].axis('off')

        if ff % (ncols*nrows) == ncols*nrows-1 or ff == n_files-1:
            ax[0,1].set_title('MMTCam : '+path0.split('/')[-2], loc=u'center',
                              fontsize=14, weight='bold')
            subplots_adjust(left=0.025, bottom=0.025, top=0.975, right=0.975,
                            wspace=0.02, hspace=0.02)
            fig.set_size_inches(8,8)
            fig.savefig(pp, format='pdf', bbox_inches='tight')
    #endfor

    pp.close()
    if silent == False: log.info('### End psf_contours: '+systime())
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
    psf_contours(files=files, path0=path0, verbose=False)

    if silent == False: log.info('### End: '+systime())
