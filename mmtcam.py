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
from datetime import datetime, timedelta # Mod on 28/02/2017
from astropy.time import Time, TimezoneInfo

from ccdproc import cosmicray_median # + on 26/02/2017

from scipy.ndimage import uniform_filter # + on 26/02/2017

from astroquery.irsa import Irsa as IRSA # + on 27/02/2017
import astropy.coordinates as coords # + on 27/02/2017
from astropy.wcs import WCS # + on 28/02/2017

import pymysql # + on 28/02/2017

out_cat_dir = 'daofind/' # + on 23/02/2017

# + on 24/02/2017 | Mod on 27/02/2017 for less padding
bbox_props = dict(boxstyle="square,pad=0.15", fc="w", alpha=0.75, ec="none")

c_levels = 0.2+0.1*np.arange(9)
f_s      = 2*np.sqrt(2*np.log(2)) # sigma-FWHM conversion | + on 25/02/2017

#Moved up on 27/02/2017
pscale   = 0.16 * u.arcsec
v_pscale = pscale.to(u.arcsec).value

# Moved up on 01/03/2017
utc_mst = TimezoneInfo(utc_offset=-7*u.hour)
mst_utc = TimezoneInfo(utc_offset=+7*u.hour)

# Convert from m/s to mph | + on 01/03/2017
mph = u.imperial.mile/u.hour
mph_conv = (1 * u.m/u.s).to(mph).value

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

def remove_dup_sources(s_cat, verbose=False):
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
                if verbose == True: print nn, 'too many : ', len(i_match)
                bad += i_match.tolist()
    return bad
#enddef

def hdr_annotate(h0, ax):
    '''
    Define string for annotation in psf_contours()

    Parameters
    ----------
    h0 : astropy.io.fits.header.Header
      FITS header

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
      FITS header

    Returns
    -------
    t.strftime : string
      Time in HH:MM:SS format

    Notes
    -----
    Created by Chun Ly, 26 February 2017
    '''

    t = Time(h0['DATE-OBS']).to_datetime(timezone=utc_mst)
    return t.strftime('%H:%M:%S')
#enddef

def draw_NE_vector(h0, ax0):
    '''
    Draw N and E vectors based on CD matrix

    Parameters
    ----------
    h0 : astropy.io.fits.header.Header
      FITS header

    Returns
    -------
    None

    Notes
    -----
    Created by Chun Ly, 26 February 2017
    '''

    cd = [h0[key] for key in ['CD1_1','CD1_2','CD2_1','CD2_2']]

    # Draw N vector
    mN = np.max(np.abs(cd[2:]))
    dxN, dyN = cd[2]/mN, cd[3]/mN

    mE = np.max(np.abs(cd[:2]))
    dxE, dyE = cd[0]/mE, cd[1]/mE

    if dxN != 1 or dxE != 1: rx, ry = -2.40, 0.0
    if dxN == 1 or dxE == 1: rx, ry = -3.00, 0.0
    if dyN == 1 and dxE == 1: rx, ry = -3.50, 0.0

    ax0.arrow(rx, ry, dxN, dyN, head_width=0.05, head_length=0.1, fc='k', ec='k')
    if np.abs(dyN/dxN) > 1:
        if dyN == -1: haN, vaN = 'left', 'top'
        if dyN == +1: haN, vaN = 'left', 'bottom'
    if np.abs(dxN/dyN) > 1:
        if dxN == -1: haN, vaN = 'right', 'bottom'
        if dxN == +1: haN, vaN = 'left', 'bottom'
    ax0.annotate('N', [rx+dxN, ry+dyN], xycoords='data', ha=haN, va=vaN)

    # Draw E vector
    ax0.arrow(rx, ry, dxE, dyE, head_width=0.05, head_length=0.1, fc='k', ec='k')
    if np.abs(dyE/dxE) > 1:
        if dyE == -1: haE, vaE = 'left', 'top'
        if dyE == +1: haE, vaE = 'left', 'bottom'
    if np.abs(dxE/dyE) > 1:
        if dxE == -1: haE, vaE = 'right', 'bottom'
        if dxE == +1: haE, vaE = 'left', 'bottom'
    ax0.annotate('E', [rx+dxE, ry+dyE], xycoords='data', ha=haE, va=vaE)

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

    max0 = np.nanmax(post) # Mod on 28/02/2017

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
#enddef

def query_mmtlog_wind(u_start, u_stop, user='webuser', passwd='', path0='',
                      silent=False, verbose=True):
    '''
    Query ops.mmto.arizona.edu's log for wind data.
    Note: This code requires specifying the password

    Parameters
    ----------
    u_start : string
     UTC start time. Formatted as 'YYYY-MM-DD HH:MM:SS'

    u_stop : string
     UTC stop time. Formatted as 'YYYY-MM-DD HH:MM:SS'

    user : string
     Username to login. Default: 'webuser'

    passwd : string
     Password for user. Default: ''

    path0 : string
      Directory path to files.

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
    tab0 : astropy.table.Table
     Astropy table containing young and young2 wind data

    Notes
    -----
    Created by Chun Ly, 28 February 2017
    Modified by Chun Ly, 01 March 2017
     - Return tab0
    '''

    if passwd == '':
        log.error('Must specify password!')
        log.error('Exiting!!!')
        return

    if silent == False: log.info('### Begin query_mmtlog_wind: '+systime())

    m_start = Time(u_start).to_datetime(timezone=utc_mst)
    m_start = m_start.strftime('%Y-%m-%d %H:%M%:%S')
    m_stop  = Time(u_stop).to_datetime(timezone=utc_mst) + \
              timedelta(seconds=10*60.0) # Add 10 min. to have enough buffer
    m_stop  = m_stop.strftime('%Y-%m-%d %H:%M%:%S')

    conn = pymysql.connect(host='ops.mmto.arizona.edu', user=user,
                           passwd=passwd, db='mmtlogs')
    cur = conn.cursor()

    sql1 = "SELECT timestamp,young_wind_speed,young_wind_direction FROM "+\
           "young_background_log where timestamp >= '"+m_start+"' AND "+\
           "timestamp < '"+m_stop+"'"
    # print sql1
    n_entries = cur.execute(sql1)
    results1 = cur.fetchall()

    sql2 = sql1.replace('young', 'young2')
    cur.execute(sql2)
    results2 = cur.fetchall()

    time0   = np.repeat('XXXX-XX-XX XX:XX:XX', n_entries)
    speed1  = np.zeros(n_entries)
    direct1 = np.zeros(n_entries)
    speed2  = np.zeros(n_entries)
    direct2 = np.zeros(n_entries)

    for nn in xrange(n_entries):
        time0[nn]   = results1[nn][0].isoformat()
        speed1[nn]  = results1[nn][1]
        direct1[nn] = results1[nn][2]

        speed2[nn]  = results2[nn][1]
        direct2[nn] = results2[nn][2]

    outfile = path0+'wind_data.tbl'
    vec0   = [time0, speed1, direct1, speed2, direct2]
    names0 = ('MST_time','speed1', 'direct1', 'speed2', 'direct2')
    tab0 = Table(vec0, names=names0)
    if silent == False: log.info('## Writing : '+outfile)
    asc.write(tab0, outfile, format='fixed_width_two_line', overwrite=True)
    if silent == False: log.info('### End query_mmtlog_wind: '+systime())
    return tab0
#enddef

def wind_avg_max(wind_tab0, h0):
    '''
    Compute average, maximum and avg direction of wind from wind data for
    each observation period

    Parameters
    ----------
    wind_tab0 : astropy.table.Table
     Astropy table containing young and young2 wind data

    h0 : astropy.io.fits.header.Header
      FITS header

    Returns
    -------
    h0 : astropy.io.fits.header.Header
      Updated FITS header

    Notes
    -----
    Created by Chun Ly, 1 March 2017
    '''
    t_start = Time(h0['DATE-OBS']).to_datetime(timezone=utc_mst)
    t_stop  = t_start + timedelta(seconds=h0['EXPTIME'])

    # Note: mjd_start/mjd_stop is MJD associated with MST time not UTC
    mjd_start = Time(t_start.strftime('%Y-%m-%d %H:%M%:%S')).mjd
    mjd_stop  = Time(t_stop.strftime('%Y-%m-%d %H:%M%:%S')).mjd

    time0 = Time(wind_tab0['MST_time'])
    mjd0  = time0.mjd
    t_idx = np.where((mjd0 >= mjd_start) & (mjd0 <= mjd_stop))[0]

    avg1 = np.average(wind_tab0['speed1'][t_idx]*mph_conv)
    max1 = np.max(wind_tab0['speed1'][t_idx]*mph_conv)
    dir1 = np.average(wind_tab0['direct1'][t_idx])

    avg2 = np.average(wind_tab0['speed2'][t_idx]*mph_conv)
    max2 = np.max(wind_tab0['speed2'][t_idx]*mph_conv)
    dir2 = np.average(wind_tab0['direct2'][t_idx])

    #young1 = {'avg':avg1, 'max':max1, 'dir':dir1}
    #young2 = {'avg':avg2, 'max':max2, 'dir':dir2}
    #return young1, young2

    # Update header with wind data
    h0.set('Y1_AVG', avg1, 'YOUNG1 avg wind speed [mph]')
    h0.set('Y1_MAX', max1, 'YOUNG1 max wind speed [mph]')
    h0.set('Y1_DIR', dir1, 'YOUNG1 avg wind direction [deg]')

    h0.set('Y2_AVG', avg2, 'YOUNG2 avg wind speed [mph]')
    h0.set('Y2_MAX', max2, 'YOUNG2 max wind speed [mph]')
    h0.set('Y2_DIR', dir2, 'YOUNG1 avg wind direction [deg]')

    return h0
#enddef

def check_extended(h0, s_cat, seqno, return_irsa_cat=False, silent=False,
                   verbose=True):
    '''
    Query the 2MASS extended source catalog (XSC) to identify contamination
    from extended sources

    Parameters
    ----------
    s_cat : astropy.table.Table
     Astropy-formatted table

    hd0 : FITS header
     FITS header containing WCS to determine RA/Dec

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 27 February 2017
    Modified by Chun Ly, 28 February 2017
     - Include extended galaxies slightly outside of MMTCam FoV
     - Use elliptical formula to determine if inside extended source
    '''

    if silent == False: log.info('### Begin check_extended: '+systime())

    size0 = 1.25*h0['NAXIS1'] * pscale # Size of region to search
    c0 = coords.SkyCoord(ra=h0['CRVAL1'], dec=h0['CRVAL2'], unit=u.deg)
    i_cat0 = IRSA.query_region(c0, catalog='fp_xsc', spatial='Box',
                             width=size0)

    flag_ext = np.zeros(len(s_cat)) # + on 28/02/2017

    # Check if daofind sources are within elliptical region of extended sources
    # + on 28/02/2017
    if len(i_cat0) == 0:
        log.info('No extended source found for : '+seqno)
    else:
        log.info('Extended sources found for : '+seqno)
        if verbose == True: print i_cat0

        w0 = WCS(h0)
        sRA, sDec = w0.wcs_pix2world(s_cat['xcentroid'], s_cat['ycentroid'], 1)

        sc = coords.SkyCoord(ra=sRA, dec=sDec, unit=u.deg)

        for cc in range(len(i_cat0)):
            ic = coords.SkyCoord(ra=i_cat0['clon'][cc], dec=i_cat0['clat'][cc],
                                 unit=(u.hour, u.deg))
            #dist0 = ic.separation(sc).to(u.arcsec).value
            dra  = (ic.ra.deg - sRA)*3600.0 * np.cos(np.radians(ic.dec.deg))
            ddec = (ic.dec.deg - sDec)*3600.0

            ang0 = np.radians(90.0-i_cat0['sup_phi'])
            maj0 = i_cat0['r_k20fe'][cc]
            min0 = maj0*i_cat0['sup_ba'][cc]
            dist0 = ((dra*np.cos(ang0) + ddec*np.sin(ang0))/maj0)**2 + \
                    ((dra*np.sin(ang0) - ddec*np.cos(ang0))/min0)**2
            ext0 = np.where(dist0 <= 1.0)[0]
            flag_ext[ext0] = 1
            # print s_cat[ext0]
    if silent == False: log.info('### End check_extended: '+systime())

    if return_irsa_cat == False:
        return flag_ext
    else: flag_ext, i_cat0
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
    Modified by Chun Ly, 28 February 2017
     - Call check_extended() to get flag indicating extended, flag_ext
    '''

    if silent == False: log.info('### Begin find_stars: '+systime())

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

    for ff in xrange(len(files)): #[34,35,36,37,38,39]: #xrange(len(files)):
        basename = os.path.basename(files[ff])
        image, hdr = fits.getdata(files[ff], header=True)
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

        # + on 28/02/2017
        flag_ext = check_extended(hdr, s_cat, seqno[ff], verbose=False)
        i_extend = np.where(flag_ext == 1)[0]
        i_point  = np.where(flag_ext == 0)[0]
        if len(i_extend) > 0:
            cat_ext = s_cat[i_extend]
            s_cat   = s_cat[i_point]
        #endif

        # Later + on 23/02/2017
        bad = remove_dup_sources(s_cat)
        if len(bad) > 0:
            cat_bad = s_cat[bad]
            s_cat.remove_rows(bad)

        if ff == 0 and silent == False: s_cat.pprint()

        # + on 23/02/2017
        out_cat = out_cat_dir0+basename.replace('.fits.gz','.tbl')
        out_cat = out_cat.replace('.fits','.tbl')
        s_cat.write(out_cat, format='ascii.fixed_width_two_line',
                    overwrite=True)

        # Write extended catalog | + on 28/02/2017
        if len(i_extend) >0:
            out_cat_ext = out_cat.replace('.tbl','.ext.tbl')
            cat_ext.write(out_cat_ext, format='ascii.fixed_width_two_line',
                          overwrite=True)

        # Later + on 23/02/2017
        if len(bad) >0 and verbose == True:
            log.info('The following will be removed : ')
            cat_bad.pprint()
        if len(bad) >0:
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
            # print z1, z2
            norm = ImageNormalize(vmin=z1, vmax=z2) #stretch=SqrtStretch())
            ax.imshow(image_sub, cmap='Greys', origin='lower', norm=norm)
            aper0.plot(color='blue', lw=1.5, alpha=0.5)
            sat_aper0.plot(color='red', lw=1.5, alpha=0.5)
            bad_aper0.plot(color='magenta', lw=1.5, alpha=0.5)

            # Label bright sources | + on 28/02/2017
            bright = np.where(s_cat['peak'] >= 0.33*max(s_cat['peak']))[0]
            for nn in bright:
                t_pos = [s_cat['xcentroid'][nn], s_cat['ycentroid'][nn]+10]
                ax.annotate(str(nn+1), t_pos, xycoords='data', ha='center',
                            va='bottom', color='b', weight='medium')

            # Mark sources excluded by extended criteria | + on 28/02/2017
            if len(i_extend) > 0:
                ax.plot(cat_ext['xcentroid'], cat_ext['ycentroid'], 'rx', linewidth=2)

            t_ann = s_date+'/'+os.path.basename(files[ff])
            ax.set_title(t_ann, loc=u'center', fontsize=14, weight='bold')
            #ax.annotate(t_ann, [0.025,0.975], xycoords='axes fraction',
            #             ha='left', va='top', bbox=bbox_props)
            ax.set_xlim([0,hdr['NAXIS1']])
            ax.set_ylim([0,hdr['NAXIS2']])
            ax.set_xlabel('X [pixels]')
            ax.set_ylabel('Y [pixels]')

            fig.set_size_inches(8,8)
            fig.savefig(pp, format='pdf', bbox_inches='tight')
    #endfor
    if plot == True:
        if silent == False:
            log.info('## Writing : '+out_pdf_plot+' | '+systime())
        pp.close()

    if silent == False: log.info('### End find_stars: '+systime())
#enddef

def make_postage(files=None, path0=None, n_stack=5, size=50,
                 user='webuser', passwd='', silent=False, verbose=True):
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
    Modified by Chun Ly, 26 February 2017
     - Use cosmicray_median() to interpolate over CRs
     - Include number of stack sources in FITS header
    Modified by Chun Ly, 28 February 2017
     - Call query_mmtlog_wind() function
     - Add user and passwd keyword to pass on
    Modified by Chun Ly, 1 March 2017
     - Check if wind data table is available before running query_mmtlog_wind()
     - Call wind_avg_max() function
    '''

    if files == None and path0 == None:
        log.error('files and path0 keywords not provided')
        log.error('Exiting!!!')
        return

    if silent == False: log.info('### Begin make_postage: '+systime())

    if files == None and path0 != None:
        files, seqno = get_files(path0)
        # path0 = None # Reset since files will have full path
    else:
        if files != None: seqno = get_seqno(files)

    # Query for wind data | + on 28/02/2017
    # Mod on 01/03/2017 to check if file exists
    wind_file = path0+'wind_data.tbl'
    if not exists(wind_file):
        u_start = fits.getheader(files[0])['DATE-OBS']
        u_stop  = fits.getheader(files[-1])['DATE-OBS']
        wind_tab0 = query_mmtlog_wind(u_start, u_stop, user=user,
                                      passwd=passwd, path0=path0)
    else:
        if silent == False: log.info('### File found! Reading : '+wind_file)
        wind_tab0 = asc.read(wind_file, format='fixed_width_two_line')

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

        # Handle failure if only one source is available and is near edge
        # + on 28/02/2017
        not_edge = np.where((s_cat['xcentroid'] > 50.0) &
                            (s_cat['xcentroid'] <= hdr['NAXIS1']-50) &
                            (s_cat['ycentroid'] > 50.0) &
                            (s_cat['ycentroid'] <= hdr['NAXIS2']-50))[0]
        s_cat = s_cat[not_edge]

        # + on 28/02/2017
        good  = np.where(s_cat['peak'] >= 0.33*max(s_cat['peak']))[0]
        s_cat = s_cat[good]

        n_bright = np.min([n_stack,len(s_cat)])
        bright   = range(n_bright)
        s_cat    = s_cat[bright]

        x0 = np.round_(s_cat['xcentroid'])
        y0 = np.round_(s_cat['ycentroid'])

        im0 = np.zeros( (len(bright), size, size) )
        size2d = u.Quantity((size, size), u.pixel)
        for ii in range(n_bright):
            pos0 = (x0[ii], y0[ii])
            cutout = Cutout2D(image_sub, pos0, size2d, mode='partial',
                              fill_value=np.nan)
            # Identify and interpolate over CRs
            cutout_cr, crmask = cosmicray_median(cutout.data, thresh=5, rbox=11)
            im0[ii] = cutout_cr/np.max(cutout_cr)

        out_fits = post_dir0+seqno[ff]+'.fits'
        psf_im = np.nanmedian(im0, axis=0)
        hdr.set('NBRIGHT', n_bright) # + on 26/02/2017

        wind_avg_max(wind_tab0, hdr) # + on 01/03/2017

        fits.writeto(out_fits, psf_im, hdr, overwrite=True)

    if silent == False: log.info('### End make_postage: '+systime())
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
    Modified by Chun Ly, 26 February 2017
     - Minor stylistic plotting changes
     - Call get_mst() to get MST time
     - Draw center of best fit
     - Call draw_NE_vector() function
    Modified by Chun Ly, 26 February 2017
     - Use cosmicray_median() to interpolate over CRs
     - Use psf_im_cr over psf_im
     - Use uniform_filter() to smooth data with a size of 3 pixels
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

    ncols, nrows = 3, 3

    n_files = len(files)
    for ff in xrange(n_files):
        psf_file = post_dir0+seqno[ff]+'.fits'
        psf_im, h0 = fits.getdata(psf_file, header=True)

        # Identify and interpolate over any extraneous CRs | + on 26/02/2017
        psf_im_cr, mask = cosmicray_median(psf_im, thresh=5, rbox=11)
        psf_im_cr /= np.max(psf_im_cr)
        psf_im_sm = uniform_filter(psf_im_cr, size=3) # + on 26/02/2017

        if ff == 0:
            shape0 = psf_im_cr.shape
            x0 = pscale*np.arange(-1*shape0[0]/2.0,shape0[0]/2.0)
            y0 = pscale*np.arange(-1*shape0[1]/2.0,shape0[1]/2.0)

        if ff % (ncols*nrows) == 0:
            fig, ax = plt.subplots(nrows, ncols)

        row, col = ff / ncols % nrows, ff % ncols

        # Later mod on 24/02/2017, 26/02/2017
        cf = ax[row,col].contourf(x0, y0, psf_im_sm, levels=c_levels,
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
        t_nstack = r'N$_{\rm stack}$=%i' % h0['NBRIGHT']
        ax[row,col].annotate(t_nstack, [0.975,0.975], weight='bold', ha='right',
                             va='top', xycoords='axes fraction', fontsize=10)

        # Compute image quality | + on 25/02/2017
        f_annot = 'UTC='+h0['UT']+'  MST='+get_mst(h0)+'\n' # + on 26/02/2017
        fwhm0, fwqm0 = fwhm_fwqm_size(psf_im_cr, pscale)
        f_annot += 'Area: FWHM=%.2f", FWQM=%.2f"\n' % (fwhm0, fwqm0)

        sigG = fwhm0/f_s/pscale.to(u.arcsec).value
        ini_guess = (1.0, 25, 25, sigG, sigG, 0.0, 0.0)
        gx = np.linspace(0,shape0[0]-1,shape0[0])
        gy = np.linspace(0,shape0[1]-1,shape0[1])
        gx, gy = np.meshgrid(gx, gy)

        psf_im_re = psf_im_cr.reshape(shape0[0]*shape0[1])
        popt, pcov = opt.curve_fit(gauss2d, (gx, gy), psf_im_re, p0=ini_guess)
        FWHMx = popt[3] * f_s * pscale.to(u.arcsec).value
        FWHMy = popt[4] * f_s * pscale.to(u.arcsec).value
        f_annot += r'2DFit: FW$_1$=%.2f", FW$_2$=%.2f", ' % (FWHMx,FWHMy)
        f_annot += r'$\theta$=%.2f' % np.degrees(popt[5])
        # Mod on 27/02/2017 to have a fill color
        ax[row,col].annotate(f_annot, [0.025,0.915], xycoords='axes fraction',
                             ha='left', va='top', fontsize=8, zorder=10,
                             bbox=bbox_props)

        # Overlay 0.25, 0.50, and 0.75 quartile contours for 2-D Gaussian fit
        # + on 25/02/2017
        f_data  = gauss2d((gx, gy), *popt)
        f_data /= np.max(f_data)
        levels  = np.array([0.25,0.5,0.75])
        CS = ax[row,col].contour(x0, y0, f_data.reshape(shape0[0],shape0[1]),
                                 colors='c', linewidth=2, cmap=None,
                                 levels=levels.tolist())
        if col==0: ax[row,col].clabel(CS, CS.levels, fmt='%.2f', inline=1,
                                      inline_spacing=0.25, fontsize=8)

        # Plot center of fit | + on 26/02/2017
        xcen = v_pscale*(-1*shape0[0]/2.0+popt[1])
        ycen = v_pscale*(-1*shape0[1]/2.0+popt[2])
        ax[row,col].plot(xcen, ycen, 'o', mfc='c', mec='none', alpha=0.5)

        hdr_annotate(h0, ax[row,col]) # + on 24/02/2017

        draw_NE_vector(h0, ax[row,col]) # + on 26/02/2017

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

def run_all(files=None, path0=None, user='webuser', passwd='',
            silent=False, verbose=True):
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
    Modified by Chun Ly, 28 February 2017
     - Add user and passwd keyword to pass on to make_postage()
    '''

    if files == None and path0 == None:
        log.error('files and path0 keywords not provided')
        log.error('Exiting!!!')
        return

    if silent == False: log.info('### Begin run_all: '+systime())

    if files == None and path0 != None:
        files, seqno = get_files(path0)
        # path0 = None # Reset since files will have full path

    if silent == False:
        log.info('The following files will be analyzed from : ')
        log.info(path0)
        for file in files: log.info(os.path.basename(file))

    find_stars(files=files, path0=path0, plot=True, verbose=False)
    make_postage(files=files, path0=path0, user=user, passwd=passwd,
                 verbose=False)
    psf_contours(files=files, path0=path0, verbose=False)

    if silent == False: log.info('### End run_all: '+systime())
