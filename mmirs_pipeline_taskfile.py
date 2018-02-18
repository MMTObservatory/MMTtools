"""
mmirs_pipeline_taskfile
====

Python code to create task files for MMIRS IDL pipeline:
http://tdc-www.harvard.edu/software/mmirs_pipeline.html

Operates in a given path containing raw files, find common files, and create
task files to execute with mmirs_pipeline

This code will create several files:
(1) 'obs_summary.tbl'
     - A summary of observations in a specified 'rawdir' with FITS header
       information

(2) MMIRS input taskfiles
     - FITS header compliant ASCII files that are provided to the
       MMIRS IDL pipeline

(3) 'IDL_input_[name].lis'
     - A list of files to run the MMIRS IDL pre-processing (i.e., non-linearity
      correction)

   Here, [name] is a naming convention that contains the name of the
   target, LS or MOS mode, and grism and filter combinations. It is defined
   from mmirs_pipeline_taskfile.

TO EXECUTE:
1. import this code:
     from MMTtools import mmirs_pipeline_taskfile


2. Call code for specified path:
     rawdir = '/path/to/raw/files/' <- Should end with forward slash
     mmirs_pipeline_taskfile.create(rawdir, w_dir='', dither='ABApBp',
                                    bright=True)

   Notes:
    1. w_dir can be changed. Default is to create a 'reduced' folder in rawdir
    2. dither: If NOT specified, code will determine dither pattern based on
       FITS header
    3. Set bright to True of False if there is a bright object in slit or MOS


3. Run the IDL script run_mmirs_pipeline_nonlin_script.idl that is automatically
   generated from step 2 in the rawdir path
     idl run_mmirs_pipeline_nonlin_script.idl

   Note: All the pre-processed files will be placed in the 'preproc' folder.


4. After creating pre-processed files, you can now run the MMIRS pipeline by
   creating an IDL script in rawdir to run run_pipeline.pro. The script should
   look like the following:
     .run run_pipeline

     run_pipeline, rawdir + 'reduced/[target_folder]/'

     exit

   Here, [target_folder] is a naming convention that contains the name of the
   target, LS or MOS mode, and grism and filter combinations. It is defined
   from mmirs_pipeline_taskfile.

   If your rawdir contains multiple targets, mmirs_pipeline_taskfile _should_
   separate out the targets in a respective manner.  Thus, the above IDL script
   will need to be created and run for each target combination.
"""

__version__ = '0.1' # Set on 16/02/2018

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

from astropy.time import Time # + on 31/01/2018

# + on 29/01/2018
from astroquery.simbad import Simbad
if not any('sptype' in sfield for sfield in Simbad.get_votable_fields()):
    Simbad.add_votable_fields('sptype')

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
    Modified by Chun Ly, 11 December 2017
     - Add instrument elevation offset (dithers along slit)
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
    instel   = np.zeros(n_files0) # + on 11/12/2017

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

        instel[ii] = hdr['INSTEL'] # + on 11/12/2017
    #endfor

    arr0   = [filename, seqno, dateobs, object0, pi, propid, imagetyp,
              aptype, exptime, airmass, aperture, filter0, disperse, instel]
    names0 = ('filename','seqno','dateobs','object','PI','PropID','imagetype',
              'aptype','exptime','airmass','aperture','filter','disperse','instel')
    tab0 = Table(arr0, names=names0)
    tab0.sort('seqno')
    return tab0
#enddef

def get_header_comments(f0):
    '''
    Extract FITS header comments. This function is to resolve issue with FITS block size,
    which prevents fits.Header.tofile() from working

    Parameters
    ----------
    f0 : list
      List containing full header as strings


    Returns
    -------
    comm0 : list
      List containing FITS header comments

    Notes
    -----
    Created by Chun Ly, 24 January 2018
    Modified by Chun Ly, 25 January 2018
     Moved up to work with read_template()
    '''

    split0 = [str0.split(' / ')[-1].replace('\n','') for str0 in f0]

    comm0 = [split0[xx] if split0[xx] != f0[xx].replace('\n','') else ''
             for xx in range(len(split0))]
    return comm0
#enddef

def read_template(longslit=False, mos=False):
    '''
    Read in MMIRS templates to populate with information

    Parameters
    ----------
    longslit : bool
      Indicate whether to use longslit template. Default: False
      Either longslit=True or mos=True

    mos : bool
      Indicate whether to use MOS template. Default: False
      Either mos=True or longslit=True

    Returns
    -------
    temp_hdr : astropy.io.fits.header.Header
      Astropy FITS-formatted header class

    Notes
    -----
    Created by Chun Ly, 29 November 2017

    Modified by Chun Ly, 30 November 2017
     - Output ordered dictionary

    Modified by Chun Ly, 11 December 2017
     - Bug fix for splitting with '='
     - Bug fix for COMMENT entries

    Modified by Chun Ly, 18 December 2017
     - Switch from dict to FITS header for simplification
     - Update documentation

    Modified by Chun Ly, 24 January 2018
     - Return string list version of template

    Modified by Chun Ly, 25 January 2018
     - Call get_header_comments() to get FITS keywords' comments
     - Return hdr0_comm, FITS keywords' comments

    Modified by Chun Ly, 17 February 2018
     - Remove return of f0 and hdr0_comm (obsolete)
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

    # + on 18/12/2017
    temp_hdr = fits.Header.fromstring("".join(f0), sep='\n')

    # hdr0_comm = get_header_comments(f0) # + on 25/01/2018

    return temp_hdr #, f0, hdr0_comm


    #keyword = np.array([str0.split('= ')[0] for str0 in f0])
    #text0   = np.array([str0.split('= ')[-1] for str0 in f0])
    #
    #i_comm0 = [ii for ii in range(len(text0)) if
    #           ('COMMENT' in text0[ii] or 'END' in text0[ii])]
    #keyword[i_comm0] = ''
    #
    #temp_dict0 = collections.OrderedDict()
    #temp_dict0['keyword'] = keyword
    #temp_dict0['text']    = text0
    #
    #return temp_dict0 #{'keyword': keyword, 'text': text0}
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
    calib_dict0 : dict
      Ordered dictionary containing information of calibration files

    Notes
    -----
    Created by Chun Ly, 1 December 2017
    Modified by Chun Ly, 5 December 2017
     - Determine and return dark frames for science exposures
    Modified by Chun Ly, 8 December 2017
     - Return ordered dict instead of individual variables
    Modified by Chun Ly, 11 December 2017
     - Use PI and PropID to further filter out comps and flats. Need a time constraint
     - Include darks for comps, flats
    Modified by Chun Ly, 26 January 2018
     - Minor code documentation
    '''
    len0 = len(tab0)

    itype0 = tab0['imagetype']
    aper0  = tab0['aperture']
    filt0  = tab0['filter']
    disp0  = tab0['disperse']

    # + on 11/12/2017
    pi     = tab0['PI']
    propid = tab0['PropID']

    t_str  = name.split('_')
    t_obj  = t_str[0]
    t_ap   = t_str[1]
    t_filt = t_str[2]
    t_disp = t_str[3]

    i_obj      = [ii for ii in range(len0) if tab0['object'][ii] == t_obj]

    # + on 11/12/2017
    t_pi     = pi[i_obj[0]]
    t_propid = propid[i_obj[0]]

    # + on 05/12/2017
    exptime    = tab0['exptime']
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
    # Mod on 11/12/2017
    i_comp = [ii for ii in range(len0) if
              (itype0[ii] == 'comp' and aper0[ii] == t_ap and
               filt0[ii] == t_filt and disp0[ii] == t_disp and \
               pi[ii] == t_pi and propid[ii] == t_propid)]
    if len(i_comp) > 2:
        log.warn('### Too many comps!!!')

    comp_str0  = ",".join(tab0['filename'][i_comp])
    comp_itime = tab0['exptime'][i_comp[0]]

    # Darks for comps
    id_comp = [ii for ii in range(len0) if
               (itype0[ii] == 'dark' and tab0['exptime'][ii] == comp_itime)]
    comp_dark = ",".join(tab0['filename'][id_comp])


    ## FLATS
    # Mod on 11/12/2017
    i_flat = [ii for ii in range(len0) if
              (itype0[ii] == 'flat' and aper0[ii] == t_ap and \
               filt0[ii] == t_filt and disp0[ii] == t_disp and \
               pi[ii] == t_pi and propid[ii] == t_propid)]
    if len(i_flat) > 2:
        log.warn('### Too many flats!!!')

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

    calib_dict0 = collections.OrderedDict()
    calib_dict0['comp_str']   = comp_str0
    calib_dict0['comp_dark']  = comp_dark # + on 11/12/2017
    calib_dict0['flat_str']   = flat_str0
    calib_dict0['flat_dark']  = flat_dark # + on 11/12/2017
    calib_dict0['dark_etime'] = dark_etime
    calib_dict0['dark_str']   = dark_str0

    return calib_dict0
#enddef

def get_tellurics(tab0, idx, comb0):
    '''
    Determining tellurics to use

    Parameters
    ----------
    tab0: astropy.table.table
      Astropy Table containing FITS header

    idx : list or np.array
      Index of entries for a given target

    comb0 : list
      List of strings that combines name, aperture, filter, and disperse

    Returns
    -------
    tell_dict0 : dict
      Dictionary containing telluric filenames, exptime, and associated darks

    Notes
    -----
    Created by Chun Ly, 25 January 2018
     - Require that str_tell be a list
    Modified by Chun Ly, 26 January 2018
     - Minor code documentation
     - Bug fix: incorrect indexing, tmp -> i_tell
    Modified by Chun Ly, 28 January 2018
     - Include exptime in separating telluric datasets
     - Get darks for each telluric datasets, str_dark
     - Return full telluric information through dictionary, tell_dict0
    Modified by Chun Ly, 29 January 2018
     - Bug fixes: len(n_tell) -> n_tell, str(etime), specify size when calling range
     - Import astroquery.simbad to get spectral type for telluric star; export
       in tell_dict0
    '''

    obj   = tab0['object']
    etime = tab0['exptime'] # + on 28/01/2018

    mmirs_setup = [b.replace(a+'_','') for a,b in zip(obj, comb0)]

    target_setup = mmirs_setup[idx[0]]

    i_tell = [xx for xx in range(len(obj)) if
              (('HD' in obj[xx] or 'HIP' in obj[xx]) and
               (mmirs_setup[xx] == target_setup))]

    # Include exptime should data with multiple exptime for same target is taken
    # Mod on 28/01/2018
    obj_etime  = [a+'_'+str(b) for a,b in zip(obj, etime)]
    tell_comb0 = list(set(np.array(obj_etime)[i_tell]))

    n_tell = len(tell_comb0)

    if n_tell == 0:
        log.warn('### No telluric data found!!!')
        str_tell = []

    str_tell   = ['']  * n_tell
    tell_etime = [0.0] * n_tell
    str_dark   = ['']  * n_tell # + on 28/01/2018
    tell_stype = ['']  * n_tell # + on 29/01/2018

    if n_tell == 1:
        log.info('### Only one telluric star is found!!!')
        log.info('### '+obj[i_tell[0]]+' '+str(etime[i_tell[0]]))

    if n_tell >= 1:
        for tt in range(n_tell):
            tmp = [xx for xx in range(len(obj)) if
                   obj_etime[xx] == tell_comb0[tt]]

            # tell_time[tt] = tab0['dateobs'][tmp[0]]
            str_tell[tt]   = ",".join(tab0['filename'][tmp])
            tell_etime[tt] = etime[tmp[0]] # + on 28/01/2018

            # + on 28/01/2018
            i_dark = [xx for xx in range(len(obj)) if
                      (tab0['imagetype'][xx] == 'dark' and
                       etime[xx] == tell_etime[tt])]
            str_dark[tt] = ",".join(tab0['filename'][i_dark])

            # + on 29/01/2018
            t_simbad = Simbad.query_object(tab0['object'][tmp[0]])
            tell_stype[tt] = t_simbad['SP_TYPE'][0].lower()
        #endfor

    # Mod on 28/01/2018, 29/01/2018
    tell_dict0 = {'name': str_tell, 'etime': tell_etime, 'dark': str_dark,
                  'stype': tell_stype}
    return tell_dict0
#enddef

def get_diff_images(tab0, idx, dither=None):
    '''
    Determining dither pattern and identify sky frame for image differencing

    Parameters
    ----------
    tab0: astropy.table.table
      Astropy Table containing FITS header

    idx : list or np.array
      Index of entries for a given target

    dither : str
      Dithering style. Either: 'ABApBp', 'ABAB', or 'ABBA'.
      If not given, determines based on dither pattern. Default: None

    Returns
    -------
    im_dict : dict
      Ordered dictionary with science and sky frames and dither positions

    Notes
    -----
    Created by Chun Ly, 24 January 2018
    Modified by Chun Ly, 26 January 2018
     - Bug fix: Incorrectly used offset_val over instel
    '''

    tab0 = tab0[idx]
    n_files = len(tab0)

    instel = tab0['instel']
    offset_val = list(set(instel))

    if dither == None:
        log.info('### Determining dither sequence...')

        if len(offset_val) == 4: dither="ABApBp"
        if len(offset_val) == 2:
            # Mod on 26/01/2018
            if instel[1] == instel[2]: dither="ABBA"
            if instel[1] != instel[2]: dither="ABAB"

        log.info('## Dither sequence is : '+dither)

    i_off = [1, -1] * (n_files/2)
    if dither != 'ABBA':
        if n_files % 2 == 1: i_off.append(-1) # Odd number correction
    i_sky = np.arange(n_files)+np.array(i_off)

    im_dict = collections.OrderedDict()
    im_dict['sci']      = tab0['filename']
    im_dict['sci2']     = tab0['filename'][i_sky] # This is the sky frame
    im_dict['dithpos']  = instel
    im_dict['dithpos2'] = instel[i_sky]

    return im_dict
#enddef

def generate_taskfile(hdr0, rawdir, w_dir, name, c_dict0, tell_dict0, tab0,
                      idx, dither=None):
    '''
    Modify the default task file template for each science exposure

    Parameters
    ----------
    hdr0 : astropy.io.fits.header.Header
      Astropy FITS-formatted header class containing task file info

    rawdir : str
      Full path for where raw files are

    w_dir : str
      Full path for where to place reduction data

    name : str
      object + aperture + filter + disperse name from organize_targets()

    c_dict0 : dict
      Ordered dictionary containing information of calibration files

    tell_dict0 : dict
      Dictionary containing telluric filenames, exptime, and associated darks

    tab0: astropy.table.table
      Astropy Table containing FITS header

    idx : list or np.array
      Index of entries for a given target

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 29 November 2017

    Modified by Chun Ly, 1 December 2017
     - Change input to accept temp0 dict

    Modified by Chun Ly, 11 December 2017
     - Pass in rawdir, w_dir, name
     - Define val0 list to update template
     - Update taskfile template and return
     - Strip spaces in keywords for exact match. Deals with W_DIR and RAW_DIR
       confusion
     - Pass in calibration dict, c_dict0
     - Simplify col1 to include common changes
     - Bug fix for dark_str
     - Bug fix for COMMENT and END entries

    Modified by Chun Ly, 18 December 2017
     - Switch modifications from ASCII to FITS header

    Modified by Chun Ly, 24 January 2018
     - Add idx input
     - Call get_diff_images()
     - Add dither keyword input
     - Add hdr0_comm input
     - Write ASCII taskfile for each science exposure

    Modified by Chun Ly, 25 January 2018
     - Remove handling BRIGHT FITS keyword (handled prior to generate_taskfile()
       call in create())

    Modified by Chun Ly, 26 January 2018
     - Minor code documentation
     - Update col1 for only tellurics that is needed
     - Update val0 for str_tell

    Modified by Chun Ly, 28 January 2018
     - Change str_tell to tell_dict0
     - Update val0 with telluric filenames and darks for telluric datasets

    Modified by Chun Ly, 29 January 2018
     - Bug fix: Incorrect name, tell_dict -> tell_dict0
     - Update spectral type for telluric star from tell_dict0['stype']

    Modified by Chun Ly, 30 January 2018
     - Remove extra keywords for telluric datasets

    Modified by Chun Ly, 31 January 2018
     - Remove RAW_DIR keyword handling. Done in create()
     - Change w_dir for each science exposure
     - Change output file for zero padding
     - Change outfile location from rawdir to w_dir to place taskfiles

    Modified by Chun Ly, 15 February 2018
     - Add preproc in rawdir path

    Modified by Chun Ly, 16 February 2018
     - Define c_hdr0 to avoid always changing hdr0
     - Bug fix: Missing ']' closure
     - Resolve issue with writing FITS headers to ASCII files

    Modified by Chun Ly, 17 February 2018
     - Remove hdr_comm0 input
    '''

    c_hdr0 = hdr0.copy()

    col1 = ['R_DIR', 'W_DIR', 'RAWEXT', 'SLIT', 'GRISM', 'FILTER',
            'DARKSCI', 'ARC', 'DARKARC', 'FLAT', 'DARKFLAT']

    n_tell = len(tell_dict0['name']) # + on 26/01/2018, Mod on 28/01/2018

    common0 = ['STAR', 'DARKST', 'STTYPE']
    for ss in range(1,n_tell+1): # Mod on 26/01/2018
        col1 += [t0 + ('%02i' % ss) for t0 in common0]

    # + on 30/01/2018
    hdr_star = [key0 for key0 in c_hdr0 if 'STAR' in key0]
    if len(hdr_star) > n_tell:
        del_idx = range(n_tell+1,len(hdr_star)+1)
        keys0 = sum([[key1+('%02i' % ii) for key1 in common0] for
                     ii in del_idx], [])
        for key in keys0:
            del c_hdr0[key]

    ## + on 01/12/2017
    #t_keyword0 = temp0['keyword']
    #t_keyword  = [str0.replace(' ','') for str0 in t_keyword0]
    #t_text    = temp0['text']

    # + on 11/12/2017
    t_str = name.split('_')
    t_obj, t_ap, t_filt, t_disp  = t_str[0], t_str[1], t_str[2], t_str[3]

    # + on 11/12/2017
    if '-long' in name:
        slit = t_ap.replace('pixel','_pixel').replace('-long','')
    else: slit = 'mos'

    # Mod on 15/02/2018
    val0 = [rawdir+'preproc/', w_dir, '.gz', slit, t_disp, t_filt,
            c_dict0['dark_str'][0], c_dict0['comp_str'], c_dict0['comp_dark'],
            c_dict0['flat_str'], c_dict0['flat_dark']]
    # Note: need to handle dark_str for different exposure time

    # + on 26/01/2018, Mod on 28/01/2018, 29/01/2018
    for ss in range(n_tell):
        val0 += [tell_dict0['name'][ss], tell_dict0['dark'][ss],
                 tell_dict0['stype'][ss]]

    # + on 11/12/2017
    for vv in range(len(val0)):
        c_hdr0[col1[vv]] = val0[vv]

    col2 = ['SCI', 'SCI2', 'DITHPOS', 'DITHPOS2']

    # + on 24/01/2018
    im_dict = get_diff_images(tab0, idx, dither=dither)

    # Write ASCII taskfiles | + on 24/01/2018
    keys1 = c_hdr0.keys()

    for ii in range(len(im_dict['sci'])):
        print log.info('### Writing taskfile for : '+im_dict['sci'][ii])
        for t_key in col2:
            c_hdr0[t_key] = im_dict[t_key.lower()][ii]

        c_hdr0['W_DIR'] = w_dir + format(ii+1, '02') + '/'

        outfile = w_dir+name+'_'+format(ii+1, '02')+'.txt' # Mod on 31/01/2018
        log.info('## Writing : '+outfile)

        # Mod on 16/02/2018
        c_hdr0.tofile(outfile, sep='\n', padding=False, overwrite=True)
        # False padding to avoid IOError: Header size (5725) is not a multiple
        # of block size (2880)
    #endfor

    return c_hdr0
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

def create(rawdir, w_dir='', dither=None, bright=False, silent=False,
           verbose=True):

    '''
    Main function to create task files to execute

    Parameters
    ----------
    rawdir : str
      Full path to MMIRS files

    w_dir : str
      Full path for where to place reduction data
      Default: based on rawdir path with 'reduced' appended

    dither : str
      Dithering style. Either: 'ABApBp', 'ABAB', or 'ABBA'.
      If not given, determines based on dither pattern. Default: None

    bright: boolean
      Indicate whether spectra has a bright target. Default: False

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
    Modified by Chun Ly, 8 December 2017
     - Call get_calib_files()
    Modified by Chun Ly, 8 December 2017
     - Add w_dir keyword input
    Modified by Chun Ly, 11 December 2017
     - Call generate_taskfile()
     - Pass calibration dict into generate_taskfile
    Modified by Chun Ly, 18 December 2017
     - Rename variables from template dict to FITS header
    Modified by Chun Ly, 22 January 2018
     - Bug fix: Write out FITS header file via fits.Header.tofile()
    Modified by Chun Ly, 23 January 2018
     - Get list containing FITS header in string (handle larger than 80 characters)
    Modified by Chun Ly, 24 January 2018
     - Get/save string list of template from read_template()
     - Call get_header_comments() and add FITS comments to list string [str0]
     - Bug fix: Missing '\n' for str_hdr
     - Pass idx to generate_taskfile()
     - Add dither keyword input
     - Pass hdr0_comm to generate_taskfile()
    Modified by Chun Ly, 25 January 2018
     - Get FITS keywords' comments from read_template() to improve efficiency
     - Remove obsolete code (incorporated into generate_taskfile())
     - Minor code documentation
     - Add bright keyword input, Update BRIGHT keyword for hdr0
     - Call get_tellurics()
     - Pass str_tell to generate_taskfile()
    Modified by Chun Ly, 28 January 2018
     - Get tell_dict0 from get_tellurics()
     - Pass tell_dict0 to generate_taskfile()
    Modified by Chun Ly, 31 January 2018
     - Get and set appropriate RAW_DIR in hdr0
     - Change w_dir to include object+aperture+filter+disperse in path
     - Create w_dir if not present
    Modified by Chun Ly, 12 February 2018
     - Generate ASCII input file for IDL pre-processing
     - Include comps, flats and associated darks for ASCII IDL pre-processing file
     - Include telluric data and associated darks for ASCII IDL pre-processing file
    Modified by Chun Ly, 15 February 2018
     - Speed improvement: Check size of tab0_outfile first before calling
       get_header_info()
    Modified by Chun Ly, 16 February 2018
     - Bug fix: Call get_header_info if tab0_outfile does not exists
    Modified by Chun Ly, 17 February 2018
     - Write IDL script for non-linear pre-processing with
       mmirs_pipeline_nonlin_script.pro
     - Remove hdr_comm0 input to generate_taskfile() call
     - Call read_template() appropriately (now only FITS header)
     - Remove hdr_str0 and hdr0_commm0 definitions
     - Remove column containing single-frame darks for IDL_input.lis
     - Change IDL_input.lis to generate for each target
     - Write IDL scripts for each target for main mmirs_pipeline
    '''

    if silent == False: log.info('### Begin create : '+systime())

    if rawdir[-1] != '/': rawdir = rawdir + '/'

    files0 = glob.glob(rawdir+'*.????.fits*')
    files0 = [t_file for t_file in files0 if 'align_' not in t_file]
    files0 = [t_file for t_file in files0 if '_stack' not in t_file]

    n_files0 = len(files0)
    if silent == False:
        log.info('### Number of FITS files found : '+str(n_files0))

    tab0_outfile = rawdir + 'obs_summary.tbl' # Moved up on 15/02/2018

    # Use ASCII catalog if available and size matches | + on 15/02/2018
    if exists(tab0_outfile):
        if silent == False: log.info('### Reading : '+tab0_outfile)
        tab0 = asc.read(tab0_outfile, format='fixed_width_two_line')
        if n_files0 != len(tab0):
            # Get header information
            tab0 = get_header_info(files0)

            # + on 30/11/2017.
            if silent == False:
                if exists(tab0_outfile):
                    log.info('## Overwriting : '+tab0_outfile)
                else:
                    log.info('## Writing : '+tab0_outfile)

            tab0.write(tab0_outfile, format='ascii.fixed_width_two_line',
                       overwrite=True)
        else:
            if silent == False: log.info('### Using existing tab0_outfile')
    else: # + on 16/02/2018
        tab0 = get_header_info(files0)

    tab0.pprint(max_lines=-1, max_width=-1)


    comb0, obj_comb0 = organize_targets(tab0)

    # Get default FITS template headers
    # + on 30/11/2017, Mod on 18/12/2017, 24/01/2018, 17/02/2018
    LS_hdr0  = read_template(longslit=True)
    mos_hdr0 = read_template(mos=True)

    # Create task files | + on 30/11/2017
    for name in obj_comb0:
        if 'HD' not in name and 'HIP' not in name and 'BD' not in name:
            if verbose == True: log.info('## Working on : '+name)

            idx   = [ii for ii in range(len(comb0)) if comb0[ii] == name]
            n_idx = len(idx)

            # Mod on 30/11/2017, 24/01/2018
            if '-long' in name:
                hdr0      = LS_hdr0.copy()
            if 'mos' in name:
                hdr0      = mos_hdr0.copy()

            datedir = Time(tab0['dateobs'][0]).datetime.strftime('%Y.%m%d')
            orig_dir = '/data/ccd/MMIRS/'+datedir+'/'
            print orig_dir
            hdr0['RAW_DIR'] = orig_dir

            # + on 25/01/2018
            if bright:
                log.info('## Bright source flag set!')
                hdr0['BRIGHT'] = 1
            else:
                log.info('## No bright source indicated!')
                hdr0['BRIGHT'] = 0

            # on 08/12/2017
            c_dict0 = get_calib_files(name, tab0)

            tell_dict0 = get_tellurics(tab0, idx, comb0) # Mod on 28/01/2018

            # Mod on 31/01/2018
            if w_dir == '': w_dir = rawdir + 'reduced/'+name+'/'

            if not exists(w_dir):
                commands.getoutput('mkdir -p '+w_dir)

            # Generate file containing science frames and darks to run
            # MMIRS IDL preprocessing
            # + on 12/02/2018
            sci_files  = tab0['filename'][idx].data

            # + on 12/02/2018
            t_comps   = c_dict0['comp_str'].split(',')
            sci_files = np.append(sci_files, t_comps)

            t_flats   = c_dict0['flat_str'].split(',')
            sci_files = np.append(sci_files, t_flats)

            # Get telluric files | + on 12/02/2018
            for tt in range(len(tell_dict0['name'])):
                t_tell    = tell_dict0['name'][tt].split(',')
                sci_files = np.append(sci_files, t_tell)

            idl_input_file = rawdir + 'IDL_input_'+name+'.lis'
            if silent == False: log.info('### Writing : '+idl_input_file)
            asc.write([sci_files], idl_input_file, format='no_header',
                      overwrite=True)

            # + on 11/12/2017, Mod on 28/01/2018
            temp1 = generate_taskfile(hdr0, rawdir, w_dir, name, c_dict0,
                                      tell_dict0, tab0, idx, dither=dither)

            # Write IDL script for main mmirs_pipeline | + on 17/02/2018
            main_script_outfile = rawdir + 'run_mmirs_pipeline_'+name+'.idl'
            if not exists(main_script_outfile):
                if silent == False: log.info('### Writing : '+main_script_outfile)

                f1 = open(main_script_outfile, 'w')
                str0 = [".run run_pipeline\n\n",
                        "run_pipeline, 'reduced/%s\n\n" % name, "exit\n"]
                f1.writelines(str0)
                f1.close()
            else:
                log.warn('### File exists! Will not overwrite : '+main_script_outfile)
        #endif
    #endfor

    # Write IDL script for non-linear pre-processing with
    # mmirs_pipeline_nonlin_script.pro | + on 17/02/2018
    script_outfile = rawdir + 'run_mmirs_pipeline_nonlin_script.idl'
    if not exists(script_outfile):
        if silent == False: log.info('### Writing : '+script_outfile)

        f0 = open(script_outfile, 'w')
        str0 = [".run mmirs_pipeline_nonlin_script\n\n",
                "mmirs_pipeline_nonlin_script, '%s', compress='.gz', /verbose\n\n" % rawdir,
                "exit\n"]
        f0.writelines(str0)
        f0.close()
    else:
        log.warn('### File exists! Will not overwrite : '+script_outfile)

    if silent == False: log.info('### End create : '+systime())
#enddef

