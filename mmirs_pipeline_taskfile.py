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
    '''

    n_files0 = len(files0)
    filename = [] #np.array(['']*n_files0)
    seqno    = []
    exptime  = np.zeros(n_files0)
    object0  = [] #np.array(['']*n_files0)
    imagetyp = [] #np.array(['']*n_files0)
    aptype   = [] #np.array(['']*n_files0)
    aperture = [] #np.array(['']*n_files0)
    filter0  = [] #np.array(['']*n_files0)
    disperse = [] #np.array(['']*n_files0)
    
    for ii in range(n_files0):
        hdr = fits.getheader(files0[ii], ext=1)

        exptime[ii] = hdr['EXPTIME']

        t_filename = hdr['FILENAME'].split('/')[-1]
        seqno.append(t_filename.split('.')[-1])
        filename.append(t_filename)
        object0.append(hdr['OBJECT'])
        imagetyp.append(hdr['IMAGETYP'])
        aptype.append(hdr['APTYPE'])
        aperture.append(hdr['APERTURE'])
        filter0.append(hdr['FILTER'])
        disperse.append(hdr['DISPERSE'])
    #endfor

    arr0   = [filename, seqno, exptime, object0, imagetyp, aptype, aperture,
              filter0, disperse]
    names0 = ('filename','seqno','exptime','object','imagetype','aptype',
              'aperture','filter','disperse')
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

    return keyword, text0
#enddef

def generate_taskfile(keyword, text0, tab0):
    '''
    Modify the default task file templates for each case

    Parameters
    ----------

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 29 November 2017
    '''

    col1 = ['RAW_DIR', 'R_DIR', 'W_DIR', 'RAWEXT', 'SLIT', 'GRISM', 'FILTER',
            'BRIGHT', 'SCI', 'SCI2', 'DITHPOS', 'DITHPOS2', 'DARKSCI', 'ARC',
            'DARKARC', 'FLAT', 'DARKFLAT']

    common0 = ['STAR', 'DARKST', 'STTYPE']
    for ss in range(1,6):
        col1 += [t0 + ('%02i' % ss) for t0 in common0]

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

    # Create task files | + on 30/11/2017
    for name in obj_comb0:
        if 'HD' not in name and 'HIP' not in name and 'BD' not in name:
            if verbose == True: log.info('## Working on : '+name)

            idx   = [ii for ii in range(len(comb0)) if comb0[ii] == name]
            n_idx = len(idx)

            longslit, mos = False, False
            if '-long' in name: longslit = True
            if 'mos' in name: mos = True

            keyword, text0 = read_template(longslit=longslit, mos=mos)


    if silent == False: log.info('### End create : '+systime())
#enddef

