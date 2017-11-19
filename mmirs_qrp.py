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

# Mod on 07/11/2017
from astropy.stats import sigma_clip
#from ccdproc import cosmicray_median

from scipy.optimize import curve_fit

# + on 12/11/2017
import astropy.units as u
pscale = 0.2012008872545049 # arcsec/pix

def gauss1d(x, a0, a, x0, sigma):
    return a0 + a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def main(rawdir, prefix, bright=False, dither='ABApBp', silent=False,
         verbose=True):

    '''
    Main function of mmirs_qrp

    Parameters
    ----------
    rawdir : str
      Full path to MMIRS files

    prefix : str
      Prefix for files to process.  Code will search for dcorr files
      in rawdir + prefix + "*dcorr.fits*"

    bright : boolean
      Indicate whether a bright star/target is in slit.  If True,
      code will use bright star to determine offsets to shift spectra.
      If False, will use FITS header to determine dithered size
      (to be implemented). Default: False

    dither : str
      Dither sequence type.  Accepts 'ABApBp', 'ABAB' (to be implemented),
      and 'ABBA' (to be implemented). Default: 'ABApBp'

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
     - Added dithering keyword
     - Handle ABA'B' dithering to get "sky" frame
     - AB subtraction for background removal
     - Add difference data cube for combine
     - Compute average from background-subtracted, image-shifted images
    Modified by Chun Ly, 14 October 2017
     - Bug in shift call. Need to specify as row,column shift values
    Modified by Chun Ly, 16 October 2017
     - Documentation added
    Modified by Chun Ly, 19 October 2017
     - Attempt CR rejection with ccdproc.cosmicray_median
    Modified by Chun Ly, 20 October 2017
     - Handle ABBA dithering to get "sky" frame
    Modified by Chun Ly,  7 November 2017
     - Use astropy.stats.sigma_clip to mask CRs
    Modified by Chun Ly, 12 November 2017
     - Get FITS header dither values, compute differences from first frame
     - Write ASCII table with dither offsets
     - Handle dithering for bright == False using FITS dither info
     - Handle output ASCII table for bright == False

     - Quality Assurance: Compute and plot FWHM
     - Plot seqno on x-axis for FWHM plot
     - Aesthetics for plots
     - Require bright=True for FWHM calculations
    Modified by Chun Ly, 17 November 2017
     - Write npz file of arrays to expedite analysis
     - Write npz file in compressed form
     - Handle masked arrays
     - Simplify npz to limit to just masked arrays
     - Read in npz file, Use npz file for sigma_clip mask
     - Bug fix: seqno handling of .gz files
    Modified by Chun Ly, 18 November 2017
     - Minor bug fix: NAXIS2 -> NAXIS1 and vice versa
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    if rawdir[-1] != '/': rawdir = rawdir + '/'

    dcorr_files = glob.glob(rawdir+prefix+'*dcorr.fits*')
    n_files     = len(dcorr_files)
    
    hdu0 = fits.getheader(dcorr_files[0])
    naxis1, naxis2 = hdu0['NAXIS1'], hdu0['NAXIS2']

    npz_file = rawdir+prefix+'.npz'
    if exists(npz_file):
        if silent == False: log.info('## Reading : '+npz_file)
        npz0 = np.load(npz_file)

    d_cube0     = np.zeros((n_files, naxis2, naxis1))
    peak_val    = np.zeros(n_files)
    shift_cube0 = np.zeros((n_files, naxis2, naxis1))

    # Mod on 12/11/2017
    dither_az   = np.zeros(n_files)
    dither_el   = np.zeros(n_files)
    diff_cube0  = np.zeros((n_files, naxis2, naxis1))

    if dither == 'ABApBp' or dither == 'ABBA':
        i_off = [1, -1] * (n_files/2)
        if n_files % 2 == 1: i_off.append(-1) # Odd number correction
        i_sky = np.arange(n_files)+np.array(i_off)
    print i_sky

    for ii in range(n_files):
        d_cube0[ii], t_hdr = fits.getdata(dcorr_files[ii], header=True)
        # Mod on 12/11/2017
        dither_az[ii] = t_hdr['INSTAZ']
        dither_el[ii] = t_hdr['INSTEL']

        ## + on 19/10/2017. Mod on 07/11/2017
        #data_crfree, crmask = cosmicray_median(d_cube0[ii], thresh=5, rbox=11)
        #d_cube0[ii] = data_crfree

    #Flag pixels that are outliers | Mod on 07/11/2017
    #arr_med = np.median(d_cube0, axis=0)
    #arr_std = np.std(d_cube0, axis=0)
    #
    #d_cube_mask = np.zeros_like(d_cube0)
    #
    #for ii in range(n_files):
    #    mask = np.where(np.absolute((d_cube0[ii] - arr_med)/arr_std) >= 4)
    #    d_cube_mask[ii][mask] = 1

    # Mod on 19/10/2017
    for ii in range(n_files):
        t_sky  = d_cube0[i_sky[ii]]
        t_diff = d_cube0[ii] - t_sky

        # Mod on 07/11/2017
        #t_sky  = np.ma.masked_array(d_cube0[i_sky[ii]],
        #                            mask=d_cube_mask[i_sky[ii]])
        #t_diff = np.ma.masked_array(d_cube0[ii], mask=d_cube_mask[ii]) - t_sky
        diff_cube0[ii] = t_diff

        # Compute offsets using bright source
        if bright == True:
            med0_row = np.median(t_diff, axis=1)

            resize   = np.repeat(med0_row, naxis1).reshape((naxis2,naxis1))
            im_test  = t_diff - resize

            med0_col     = np.median(im_test, axis=0)
            peak_val[ii] = np.argmax(med0_col)
        else:
            if ii == 0: log.info('### Using FITS dither information')

    dither_diff = (dither_el[0] - dither_el) / pscale * u.pix

    # Mod on 12/11/2017
    if bright == True:
        shift_val = peak_val[0] - peak_val
    else:
        shift_val = dither_diff.value

    log.info('### Shift values for spectra : ')
    # Mod on 12/11/2017
    files_short = [str0.replace(rawdir,'') for str0 in dcorr_files]
    if bright == True:
        diff0  = dither_diff - shift_val * u.pix
        arr0   = [files_short, dither_az * u.arcsec, dither_el * u.arcsec,
                  dither_diff, shift_val * u.pix, diff0]
        names0 = ('file','dither_az','dither_el','dither_diff','shift_val',
                  'difference')
    else:
        arr0   = [files_short, dither_az * u.arcsec, dither_el * u.arcsec,
                  dither_diff]
        names0 = ('file','dither_az','dither_el','dither_diff')

    dither_tab = Table(arr0, names=names0)
    dither_tab.pprint(max_lines=-1)

    # + on 12/11/2017
    out_dither_file1 = rawdir+prefix+'_dither_ecsv.cat'
    if silent == False: log.info('## Writing : '+out_dither_file1)
    dither_tab.write(out_dither_file1, format='ascii.ecsv')

    # + on 12/11/2017
    out_dither_file2 = rawdir+prefix+'_dither.cat'
    if silent == False: log.info('## Writing : '+out_dither_file2)
    dither_tab.write(out_dither_file2, format='ascii.fixed_width_two_line')
    #log.info('### '+" ".join([str(a) for a in shift_val]))

    for ii in range(n_files):
        # Bug fix - Mod on 14/10/2017
        shift_cube0[ii] = shift(diff_cube0[ii], [0,shift_val[ii]])

    # + on 07/11/2017
    if not exists(npz_file):
        shift_cube0_mask = sigma_clip(shift_cube0, sigma=3., iters=5, axis=0)
    else:
        shift_cube0_mask = np.ma.array(shift_cube0, mask=npz0['shift_cube0_mask'])

    stack0 = np.ma.average(shift_cube0_mask, axis=0)

    # + on 17/11/2017
    if silent == False: log.info('## Writing : '+npz_file)
    np.savez_compressed(npz_file, dither_tab=dither_tab,
                        shift_cube0_mask=shift_cube0_mask.mask)
    #stack0d=stack0.data, stack0m=stack0.mask,

    # Mod on 07/11/2017
    fits.writeto(rawdir+prefix+'_stack.fits', stack0.data, overwrite=True)
    #fits.writeto(rawdir+prefix+'_stack.fits', stack0, overwrite=True)

    # Compute FWHM and plot | + on 12/11/2017
    if bright == True:
        out_fwhm_pdf = rawdir+prefix+'_stack_FWHM.pdf'
        fig, ax = plt.subplots()

        FWHM0 = np.zeros(n_files)
        seqno = [str0.replace('.gz','').replace('_dcorr.fits','')[-4:] for \
                 str0 in dcorr_files] # Mod on 17/11/2017

        for ii in range(n_files):
            im0    = shift_cube0_mask[ii]
            med0   = np.ma.median(im0, axis=0)
            x0     = np.arange(len(med0))
            x0_max = np.argmax(med0)

            p0         = [0.0, max(med0), x0_max, 2.0]
            popt, pcov = curve_fit(gauss1d, x0, med0, p0=p0)
            FWHM0[ii]  = popt[3]*2*np.sqrt(2*np.log(2)) * pscale

        ax.plot(seqno, FWHM0, marker='o', color='b', alpha=0.5)
        ax.set_xlabel('Image Frame No.')
        ax.set_ylabel('FWHM [arcsec]')
        ax.minorticks_on()
        fig.set_size_inches(8,8)
        fig.savefig(out_fwhm_pdf, bbox_inches='tight')

    if silent == False: log.info('### End main : '+systime())
#enddef

