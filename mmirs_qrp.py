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

# + on 20/11/2017
from matplotlib.backends.backend_pdf import PdfPages
from pylab import subplots_adjust

# + on 12/11/2017
import astropy.units as u

pscale = 0.2012008872545049 # arcsec/pix

bbox_props = dict(boxstyle="square,pad=0.15", fc="w", alpha=0.75, ec="none")

def gauss1d(x, a0, a, x0, sigma):
    return a0 + a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def main(rawdir, prefix, bright=False, dither='ABApBp', flats=[],
         silent=False, verbose=True):

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

    flats: list
      List of files or seqno for flats.
      For example, flats=[1100,1101] or flats=['flat.1100','flat.1101']
      If not provided, flat fielding will not be performed

      NOTE: THERE ARE SOME BUGS WITH FLATFIELDING

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
    Modified by Chun Ly, 20 November 2017
     - Plot fit to PSF profiles
    Modified by Chun Ly, 25 November 2017
     - PSF profiles showed double peak with a lower peak that is 20% of max,
       but visual inspection does not indicate double-peaked PSF.
     - Unclear of the full cause, but ultimately computed FWHM in the central
       200 pix without masking (masking seems to reject pixels from bright objects)
     - Masking was part of the cause as it would produce a double peak
       distribution
    Modified by Chun Ly, 26 November 2017
     - Run curve_fit to get accurate line center for shifting
     - Set default shifting as integer pixels since CRs are present
       in stack (may need to grow mask)
    Modified by Chun Ly, 9 December 2017
     - Compute transparency: Integrate flux for bright source and plot
     - Normalize spectrum by exposure time for proper transparency computation
     - Plotting aesthetic improvements: legend, labels
     - Normalize transparency value to best, plotting aesthetics
     - Plotting aesthetic improvements: different linestyle and widths, smaller
       legend

     - Incorporate flat fielding (still some unexpected problems)

    Modified by Chun Ly, 11 December 2017
     - Expand dither_tab stdout,
     - Add annotation of avg/median/sigma for FWHM
     - Handle unused FWHM subplots
     - Minor stdout for avg/median/sigma for FWHM
     - Annotation for source info in FWHM and transparency plots
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
    peak_val0   = np.zeros(n_files) # + on 26/11/2017
    shift_cube0 = np.zeros((n_files, naxis2, naxis1))

    exptime0 = np.zeros(n_files) # + on 09/12/2017

    # Define and combine flat files | + on 09/12/2017
    do_flat = 0
    if len(flats) != 0:
        if 'flat' in str(flats[0]):
            flats = [flat+'_dcorr.fits' for flat in flats]
        else:
            flats = ['flat.%04i_dcorr.fits' % t_seq for t_seq in flats]

        flat_arr = np.zeros((len(flats), naxis2, naxis1))
        for ff in range(len(flats)):
            t_flat, flat_hdr = fits.getdata(rawdir+flats[ff], header=True)

            flat_arr[ff] = t_flat / np.max(t_flat)

        flat_arr_mask = sigma_clip(flat_arr, sigma=3., iters=3, axis=0)

        flat_avg = np.ma.average(flat_arr_mask, axis=0)
        flat0 = flat_avg.data

        bad = np.where((np.isfinite(flat0) == False) | (flat0 == 0) |
                       (flat0 == np.min(flat0)))
        if len(bad) > 0: flat0[bad] = 1.0
        # print np.min(flat0), np.max(flat0)

        out_fits_file = rawdir+prefix+'_flat.fits'
        log.info('### Writing : '+out_fits_file)
        fits.writeto(out_fits_file, flat0, flat_hdr, overwrite=True)
        do_flat = 1

    # Set this to use curvefit to compute fractional offset
    do_curvefit_center = 0 # + on 26/11/2017

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
        d_data, t_hdr = fits.getdata(dcorr_files[ii], header=True)
        # print ii, np.min(d_data), np.max(d_data)

        if do_flat: d_data = d_data / flat0
        # print ii, np.min(d_data), np.max(d_data)

        d_cube0[ii] = d_data

        # Mod on 12/11/2017
        dither_az[ii] = t_hdr['INSTAZ']
        dither_el[ii] = t_hdr['INSTEL']

        exptime0[ii] = t_hdr['EXPTIME'] # + on 09/12/2017

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

            # + on 26/11/2017
            if do_curvefit_center:
                p0 = [0.0, np.max(med0_col), peak_val[ii], 2.0]
                x0 = np.arange(len(med0_col))
                popt, pcov    = curve_fit(gauss1d, x0, med0_col, p0=p0)
                peak_val0[ii] = popt[2]
        else:
            if ii == 0: log.info('### Using FITS dither information')

    dither_diff = (dither_el[0] - dither_el) / pscale * u.pix

    # Mod on 12/11/2017
    if bright == True:
        if do_curvefit_center: # Mod on 26/11/2017
            shift_val = peak_val0[0] - peak_val0 # Mod on 26/11/2017
        else:
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
    dither_tab.pprint(max_lines=-1, max_width=-1)

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

    if bright == True:
        # Compute FWHM and plot | + on 12/11/2017
        out_fwhm_pdf = rawdir+prefix+'_stack_FWHM.pdf'

        # + on 20/11/2017
        pp = PdfPages(out_fwhm_pdf)
        ncols, nrows = 3, 2

        FWHM0 = np.zeros(n_files)
        seqno = [str0.replace('.gz','').replace('_dcorr.fits','')[-4:] for \
                 str0 in dcorr_files] # Mod on 17/11/2017

        for ii in range(n_files):
            if ii % (ncols*nrows) == 0: # + on 20/11/2017
                fig, ax = plt.subplots(nrows, ncols)

            row, col = ii / ncols % nrows, ii % ncols # + on 20/11/2017

            im0    = shift_cube0_mask[ii].data
            med0   = np.median(im0[1024-100:1024+100], axis=0) # Mod on 25/11/2017
            x0     = np.arange(len(med0))
            x0_max = np.argmax(med0)
            y0_max = np.max(med0)

            p0         = [0.0, y0_max, x0_max, 2.0]

            popt, pcov = curve_fit(gauss1d, x0, med0, p0=p0)
            FWHM0[ii]  = popt[3]*2*np.sqrt(2*np.log(2)) * pscale

            # + on 20/11/2017
            ax[row,col].plot((x0-popt[2])*pscale, med0/y0_max, color='black')
            ax[row,col].plot((x0-popt[2])*pscale, gauss1d(x0, *popt)/y0_max,
                             color='blue', alpha=0.5)

            # + on 20/11/2017
            if row == nrows-1:
                ax[row,col].set_xlabel('X [arcsec]')
            else:
                if ((n_files-1)-ii) > ncols-1:
                    ax[row,col].set_xticklabels([])

            # + on 11/12/2017
            if ii == n_files-1:
                for cc in range(ncols): ax[row,cc].set_xlabel('X [arcsec]')

            # + on 20/11/2017
            if col == 0:
                ax[row,col].set_ylabel('Normalized Flux')
            else: ax[row,col].set_yticklabels([])

            # + on 20/11/2017
            ax[row,col].annotate(seqno[ii], [0.025,0.975], ha='left',
                                 va='top', xycoords='axes fraction',
                                 weight='bold', fontsize=10)
            ax[row,col].annotate('FWHM = %.2f"' % FWHM0[ii], [0.975,0.975],
                                 ha='right', va='top', xycoords='axes fraction',
                                 weight='bold', fontsize=10)

            # + on 20/11/2017
            ax[row,col].set_xlim([-2.5,2.5])
            ax[row,col].set_ylim([-0.05,1.1])

            # + on 11/12/2017
            if ii == n_files-1:
                for cc in range(col+1,ncols): ax[row,cc].axis('off')
                for rr in range(row+1,nrows):
                    for cc in range(ncols): ax[rr,cc].axis('off')

            # + on 20/11/2017
            if ii % (ncols*nrows) == ncols*nrows-1 or ii == n_files-1:
                subplots_adjust(left=0.025, bottom=0.025, top=0.975,
                                right=0.975, wspace=0.02, hspace=0.02)
                fig.set_size_inches(8,6)
                fig.savefig(pp, format='pdf', bbox_inches='tight')
        #endfor

        fig, ax = plt.subplots() # + on 20/11/2017
        ax.plot(seqno, FWHM0, marker='o', color='b', alpha=0.5)
        ax.set_xlabel('Image Frame No.')
        ax.set_ylabel('FWHM [arcsec]')
        ax.minorticks_on()

        filt0, disp0, ap0 = t_hdr['FILTER'], t_hdr['DISPERSE'], \
                            t_hdr['APERTURE']

        label0 = '%s   %s   %s   %s' % (prefix, filt0, disp0, ap0)
        ax.annotate(label0, [0.025,0.975], ha='left', va='top', weight='bold',
                    xycoords='axes fraction', fontsize=11, bbox=bbox_props)

        # + on 11/12/2017
        avg_FWHM0 = np.average(FWHM0)
        med_FWHM0 = np.median(FWHM0)
        sig_FWHM0 = np.std(FWHM0)
        txt0  = 'Average : %.2f"\n Median : %.2f"\n' % (avg_FWHM0, med_FWHM0)
        txt0 += r'$\sigma$ : %.2f"' % sig_FWHM0
        if silent == False:
            log.info('Average : %.2f"' % avg_FWHM0)
            log.info('Median  : %.2f"' % med_FWHM0)
            log.info('Sigma   : %.2f"' % sig_FWHM0)

        ax.annotate(txt0, [0.975,0.975], ha='right', va='top',
                    xycoords='axes fraction', fontsize=11, bbox=bbox_props)

        # Mod on 20/11/2017
        fig.set_size_inches(8,6)
        fig.savefig(pp, format='pdf', bbox_inches='tight')

        # + on 20/11/2017
        if silent == False:
            log.info('## Writing : '+out_fwhm_pdf+' | '+systime())
            pp.close()

        # Compute transparency and plot | + on 09/12/2017
        out_trans_pdf = rawdir+prefix+'_stack_trans.pdf'

        pp = PdfPages(out_trans_pdf)

        trans0 = np.zeros(n_files)
        spec0  = np.zeros((n_files, naxis2))

        fig, ax = plt.subplots()

        lstyle0 = ['solid','dashed','dashdot','dotted',
                   'solid','dashed','dashdot','dotted']
        lw0 = [0.25, 0.33, 0.50, 0.75, 0.25, 0.33, 0.50, 0.75]

        for ii in range(n_files):
            x0    = np.arange(naxis1)
            t_sig = FWHM0[ii] / (2*np.sqrt(2*np.log(2))) / pscale
            idx   = np.where(np.abs(x0 - peak_val[ii])/t_sig <= 2.5)[0]

            im0        = diff_cube0[ii] / exptime0[ii]
            spec0[ii]  = np.sum(im0[:,idx], axis=1)
            trans0[ii] = np.sum(im0[:,idx])

            lstyle, lw = lstyle0[ii / 7], lw0[ii / 7]
            ax.plot(x0, spec0[ii], linewidth=lw, label=seqno[ii], alpha=0.5,
                    linestyle=lstyle)

        ax.set_xlim([0,2100])
        ax.set_xlabel('X [pixels]')
        ax.set_ylabel('Flux [ADU/s]')
        ax.legend(loc='upper left', fontsize='9', ncol=3, framealpha=0.5,
                  columnspacing=0.5)

        fig.set_size_inches(8,6)
        fig.savefig(pp, format='pdf', bbox_inches='tight')

        fig, ax = plt.subplots() # + on 20/11/2017
        i_max = np.argmax(trans0)
        ax.plot(seqno, trans0/trans0[i_max], marker='o', color='b', alpha=0.5)
        ax.set_xlabel('Image Frame No.')
        ax.set_ylabel('Relative Transparency')
        ax.set_ylim([0,1.1])
        ax.annotate(label0, [0.025,0.975], ha='left', va='top', weight='bold',
                    xycoords='axes fraction', fontsize=11, bbox=bbox_props)
        ax.minorticks_on()

        # Mod on 20/11/2017
        fig.set_size_inches(8,6)
        fig.savefig(pp, format='pdf', bbox_inches='tight')

        # + on 20/11/2017
        if silent == False:
            log.info('## Writing : '+out_trans_pdf+' | '+systime())
            pp.close()

    if silent == False: log.info('### End main : '+systime())
#enddef

