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

import collections

# + on 30/11/2017
co_filename = __file__
co_path     = os.path.dirname(co_filename) + '/'

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
     - List for variables with strings
     - Get proper FITS extension
     - Add [seqno] for sorting purposes
    Modified by Chun Ly, 30 November 2017
     - Add airmass info
     - Add date-obs info
     - Re-order table columns
    Modified by Chun Ly, 8 December 2017
     - Add PI and PropID info
    '''

    n_files0 = len(files0)
    filename = [] #np.array(['']*n_files0)
    seqno    = []
    exptime  = np.zeros(n_files0)
    airmass  = np.zeros(n_files0)
    dateobs  = [] #np.array(['']*n_files0)
    object0  = [] #np.array(['']*n_files0)
    imagetyp = [] #np.array(['']*n_files0)
    aptype   = [] #np.array(['']*n_files0)
    aperture = [] #np.array(['']*n_files0)
    filter0  = [] #np.array(['']*n_files0)
    disperse = [] #np.array(['']*n_files0)
    pi       = [] # + on 08/12/2017
    propid   = [] # + on 08/12/2017

    for ii in range(n_files0):
        hdr = fits.getheader(files0[ii], ext=1)

        exptime[ii] = hdr['EXPTIME']
        airmass[ii] = hdr['AIRMASS']

        t_filename = hdr['FILENAME'].split('/')[-1]
        seqno.append(t_filename.split('.')[-1])
        filename.append(t_filename)
        dateobs.append(hdr['DATE-OBS'])
        object0.append(hdr['OBJECT'])
        imagetyp.append(hdr['IMAGETYP'])
        aptype.append(hdr['APTYPE'])
        aperture.append(hdr['APERTURE'])
        filter0.append(hdr['FILTER'])
        disperse.append(hdr['DISPERSE'])
        pi.append(hdr['PI']) # + on 08/12/2017
        propid.append(hdr['PROPID']) # + on 08/12/2017
    #endfor

    arr0   = [filename, seqno, dateobs, object0, pi, propid, imagetyp,
              aptype, exptime, airmass, aperture, filter0, disperse]
    names0 = ('filename','seqno','dateobs','object','PI','PropID','imagetype',
              'aptype','exptime','airmass','aperture','filter','disperse')
    tab0 = Table(arr0, names=names0)
    tab0.sort('seqno')
    return tab0
#enddef

def read_template(longslit=False, mos=False):
    '''
    Read in MMIRS templates to populate with information

    Parameters
    ----------
    longslit : bool
      Indicate whether to use longslit template. Default: False

    mos : bool
      Indicate whether to use MOS template. Default: False

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 29 November 2017

    Modified by Chun Ly, 30 November 2017
     - Output ordered dictionary
    '''

    if longslit == False and mos == False:
        log.warn('## Must specify longslit or mos keyword!!!')
        log.warn('## Exiting!')
        return

    if longslit == True:
        temp_file = co_path + 'mmirs_longslit_template.txt'

    if mos == True:
        temp_file = co_path + 'mmirs_mos_template.txt'

    #log.info('## Reading : '+temp_file)
    f = open(temp_file)

    f0 = f.readlines()

    keyword = [str0.split('=')[0] for str0 in f0]
    text0   = [str0.split('=')[-1] for str0 in f0]

    temp_dict0 = collections.OrderedDict()
    temp_dict0['keyword'] = keyword
    temp_dict0['text']    = text0

    return temp_dict0 #{'keyword': keyword, 'text': text0}
#enddef

def get_calib_files(name, tab0):
    '''
    Get appropriate calibration files for each dataset

    Parameters
    ----------
    name : str
      object + aperture + filter + disperse name from organize_targets()

    tab0: astropy.table.table
      Astropy Table containing FITS header info

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 1 December 2017
    Modified by Chun Ly, 5 December 2017
     - Determine and return dark frames for science exposures
    '''
    len0 = len(tab0)

    itype0 = tab0['imagetype']
    aper0  = tab0['aperture']
    filt0  = tab0['filter']
    disp0  = tab0['disperse']

    t_str  = name.split('_')
    t_obj  = t_str[0]
    t_ap   = t_str[1]
    t_filt = t_str[2]
    t_disp = t_str[3]

    # + on 05/12/2017
    exptime    = tab0['exptime']
    i_obj      = [ii for ii in range(len0) if tab0['object'][ii] == t_obj]
    dark_etime = list(set(np.array(exptime)[i_obj]))
    dark_str0  = []
    for etime in dark_etime:
        i_dark = [ii for ii in range(len0) if
                  (itype0[ii] == 'dark' and exptime[ii] == etime)]
        print i_dark
        t_txt = ",".join(tab0['filename'][i_dark])
        dark_str0.append(t_txt)
        log.info("## List of science dark files for %.3fs : %s" % (etime, t_txt))

    ## COMPS
    i_comp = [ii for ii in range(len0) if
              (itype0[ii] == 'comp' and aper0[ii] == t_ap and \
               filt0[ii] == t_filt and disp0[ii] == t_disp)]

    comp_str0  = ",".join(tab0['filename'][i_comp])
    comp_itime = tab0['exptime'][i_comp[0]]

    # Darks for comps
    id_comp = [ii for ii in range(len0) if
               (itype0[ii] == 'dark' and tab0['exptime'][ii] == comp_itime)]
    comp_dark = ",".join(tab0['filename'][id_comp])


    ## FLATS
    i_flat = [ii for ii in range(len0) if
              (itype0[ii] == 'flat' and aper0[ii] == t_ap and \
               filt0[ii] == t_filt and disp0[ii] == t_disp)]


    flat_str0 = ",".join(tab0['filename'][i_flat])
    flat_itime = tab0['exptime'][i_flat[0]]

    # Darks for flats
    id_flat = [ii for ii in range(len0) if
               (itype0[ii] == 'dark' and tab0['exptime'][ii] == flat_itime)]
    flat_dark = ",".join(tab0['filename'][id_flat])

    log.info("## List of comp files : "+comp_str0)
    log.info("## List of comp dark files : "+comp_dark)

    log.info("## List of flat files : "+flat_str0)
    log.info("## List of flat dark files : "+flat_dark)

    return comp_str0, flat_str0, dark_etime, dark_str0
#enddef

def generate_taskfile(temp0, tab0):
    '''
    Modify the default task file template for each science exposure

    Parameters
    ----------
    temp0 : dict
      Dictionary containing task file info

    Returns
    -------
    tab0: astropy.table.table
      Astropy Table containing FITS header

    Notes
    -----
    Created by Chun Ly, 29 November 2017

    Modified by Chun Ly, 1 December 2017
     - Change input to accept temp0 dict
    '''

    col1 = ['RAW_DIR', 'R_DIR', 'W_DIR', 'RAWEXT', 'SLIT', 'GRISM', 'FILTER',
            'BRIGHT', 'SCI', 'SCI2', 'DITHPOS', 'DITHPOS2', 'DARKSCI', 'ARC',
            'DARKARC', 'FLAT', 'DARKFLAT']

    common0 = ['STAR', 'DARKST', 'STTYPE']
    for ss in range(1,6):
        col1 += [t0 + ('%02i' % ss) for t0 in common0]

    # + on 01/12/2017
    t_keyword = temp0['keyword']
    t_text    = temp0['text']


#enddef

def organize_targets(tab0):
    '''
    Use FITS header information to organize targets based on name, aperture,
    filter, and disperse

    Parameters
    ----------
    tab0: astropy.table.table
      Astropy Table containing FITS header info

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 28 November 2017
     - Get unique lists of combinations
    '''

    len0 = len(tab0)

    itype = tab0['imagetype']
    nondark = [ii for ii in range(len0) if itype[ii] == 'dark']
    obj     = [ii for ii in range(len0) if itype[ii] == 'object']
    comp    = [ii for ii in range(len0) if itype[ii] == 'comp']
    flat    = [ii for ii in range(len0) if itype[ii] == 'flat']

    comb0 = ['N/A'] * len0

    for kk in [obj, comp, flat]:
        for ii in kk:
            tab0_o = tab0[ii]
            comb0[ii] = tab0_o['object'] + '_' + tab0_o['aperture'] + '_' + \
                        tab0_o['filter'] + '_' + tab0_o['disperse']

    obj_comb0 = list(set(np.array(comb0)[obj]))

    n_obj_comb0 = len(obj_comb0)
    log.info('## Total number of combinations found : '+str(n_obj_comb0))
    for oo in range(n_obj_comb0):
        print '## '+obj_comb0[oo]

    return comb0, obj_comb0
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
     - Include .gz files
     - Call get_header_info
     - Bug fix: Ignore align_ and QRP _stack files
     - Call organize_targets
    Modified by Chun Ly, 29 November 2017
     - Write ASCII table containing header info -> 'obs_summary.tbl'
     - Begin for loop for writing MMIRS pipeline control task files
     - Call read_template function
    Modified by Chun Ly, 30 November 2017
     - Get template dict from read_template()
    '''
    
    if silent == False: log.info('### Begin create : '+systime())

    if rawdir[-1] != '/': rawdir = rawdir + '/'

    files0 = glob.glob(rawdir+'*.????.fits*')
    files0 = [t_file for t_file in files0 if 'align_' not in t_file]
    files0 = [t_file for t_file in files0 if '_stack' not in t_file]

    n_files0 = len(files0)
    if silent == False:
        log.info('### Number of FITS files found : '+str(n_files0))

    # Get header information
    tab0 = get_header_info(files0)
    tab0.pprint(max_lines=-1, max_width=-1)

    # + on 30/11/2017
    tab0_outfile = rawdir + 'obs_summary.tbl'
    if silent == False:
        if exists(tab0_outfile):
            log.info('## Overwriting : '+tab0_outfile)
        else:
            log.info('## Writing : '+tab0_outfile)
    #endif
    tab0.write(tab0_outfile, format='ascii.fixed_width_two_line',
               overwrite=True)

    comb0, obj_comb0 = organize_targets(tab0)

    # Get default templates | + on 30/11/2017
    longslit_temp0 = read_template(longslit=True)
    mos_temp0      = read_template(mos=True)

    # Create task files | + on 30/11/2017
    for name in obj_comb0:
        if 'HD' not in name and 'HIP' not in name and 'BD' not in name:
            if verbose == True: log.info('## Working on : '+name)

            idx   = [ii for ii in range(len(comb0)) if comb0[ii] == name]
            n_idx = len(idx)

            # Mod on 30/11/2017
            if '-long' in name: temp0 = longslit_temp0.copy()
            if 'mos' in name: temp0 = mos_temp0.copy()


    if silent == False: log.info('### End create : '+systime())
#enddef

