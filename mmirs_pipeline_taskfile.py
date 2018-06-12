"""
mmirs_pipeline_taskfile
====

Python code to create task files for MMIRS IDL pipeline:
http://tdc-www.harvard.edu/software/mmirs_pipeline.html

Operates in a given path containing raw files, find common files, and create
task files to execute with MMIRS IDL data reduction pipeline, mmirs-pipeline

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

0. First, you will need python v2.7 (This code will not work python 3.xx).
   I recommend installing through anaconda:
     https://www.anaconda.com/download/

   Next you will need the Astropy package.  This code has been tested to work
   with v1.3 and v2.0.2 of Astropy. Install via the conda command:
     conda install astropy

   In addition, you will need astroquery package. Install via pip command:
     pip install astroquery


1. After python, astropy, and astroquery are successfully installed, you
   need to clone Chun Ly's MMTtools package to your local machine:
     git clone https://github.com/astrochun/MMTtools.git


2. Within ipython (or python) import this code:
     from MMTtools import mmirs_pipeline_taskfile

   Note: If this does not work, MMTtools is not within your PYTHONPATH
         environment


3. Remove FITS files from your raw path that you do not want this code to detect.
   This code does a file search for '*.????.fits*'.  If there are bad FITS
   files or those that are saturated, relocate them to another folder or
   delete them.


4. Call code for specified path in ipython (or python):
     rawdir = '/path/to/raw/files/' <- Should end with forward slash
     mmirs_pipeline_taskfile.create(rawdir, w_dir='', dither='ABApBp',
                                    bright=True)

   Notes:
    1. w_dir can be changed. Default is to create a 'reduced' folder in [rawdir]
    2. dither: If NOT specified, code will determine dither pattern based on
       FITS header
    3. Set bright to True of False if there is a bright object in slit or MOS
    4. Note that mmirs-pipeline will look for a 'calib_MMIRS' folder. This
       is needed in the pre-processing (e.g., applying non-linearity correction).
       This Python code will prompt you to provide a path such that a symbolic
       link is created.

       For example, if mmirs-pipeline is installed in /codes/idl, then you would
       type '/codes/idl', and it would create a symbolic link as follow:
         ln -s /codes/idl/mmirs-pipeline/pipeline/calib_MMIRS [rawdir]/calib_MMIRS

   If your [rawdir] contains multiple targets, mmirs_pipeline_taskfile _should_
   separate out the targets in a respective manner.


5. Next install IDL and set-up it up to have an appropriate license.
   Then have Igor Chilingarian's mmirs-pipeline on your computer and included
   in your IDL_PATH environment:
     git clone https://bitbucket.org/chil_sai/mmirs-pipeline

   Also make sure that you have the IDL Astrolib installed and it is in your
   IDL_PATH environment:
     git clone https://github.com/wlandsman/IDLAstro

   Note: legend.pro has been deprecated (part of IDL as al_legend.pro),
   which will cause mmirs-pipeline to crash. Either change mmirs-pipeline to
   use al_legend or include this legend.pro somewhere in your IDL_PATH:
     https://idlastro.gsfc.nasa.gov/ftp/obsolete/legend.pro


6. Run the IDL script run_mmirs_pipeline_nonlin_script.idl that is
   automatically generated from step 3 in the [rawdir] path:
     idl run_mmirs_pipeline_nonlin_script.idl

   Note: All the pre-processed files will be placed in the 'preproc' folder
   within [rawdir].


7. After creating pre-processed files, you can now run the MMIRS pipeline via
   the IDL scripts (run_mmirs_pipeline_[name].idl) that are automatically
   generated from step 3 in the [rawdir] path

   If your [rawdir] contains multiple targets, mmirs_pipeline_taskfile _should_
   separate out the targets in a respective manner.  Thus, there should be
   multiple run_mmirs_pipeline.idl scripts



TIPS:
   This code has a log file that is created in rawdir called
   'mmirs_pipeline_taskfile.log'. It logs everything that is written to
   stdout

   If a bug is encountered please submit an issue ticket here:
     https://github.com/astrochun/MMTtools/issues

   Also, please email the creator, Chun Ly, at chunly [at] mmto.org
   your mmirs_pipeline_taskfile.log, any error messages on the ipython
   or python screen, and your 'obs_summary.tbl' file

"""

__version__ = '0.1' # Set on 16/02/2018

import sys, os

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import glob

from astropy.table import Table
from astropy import log

from astropy.time import Time # + on 31/01/2018

from socket import gethostname

# Mod on 01/05/2018
use_simbad = 0
if use_simbad:
    from astroquery.simbad import Simbad
    if not any('sptype' in sfield for sfield in Simbad.get_votable_fields()):
        Simbad.add_votable_fields('sptype')

import collections

import logging
formatter = logging.Formatter('%(asctime)s - %(module)12s.%(funcName)20s - %(levelname)s: %(message)s')
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
sh.setFormatter(formatter)

# + on 30/11/2017
co_filename = __file__
co_path     = os.path.dirname(co_filename) + '/'
# Fix if code is placed locally | Mod on 02/03/2018
if co_path == '/': co_path = ''

class mlog:
    '''
    Main class to log information to stdout and ASCII file

    To execute:
    mylog = mlog(rawdir)._get_logger()

    Parameters
    ----------
    rawdir : str
      Full path for where raw files are

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 19 February 2018
    '''

    def __init__(self,rawdir):
        self.LOG_FILENAME = rawdir + 'mmirs_pipeline_taskfile.log'
        self._log = self._get_logger()

    def _get_logger(self):
        loglevel = logging.INFO
        log = logging.getLogger(self.LOG_FILENAME) # + Mod on 14/12/2017
        if not getattr(log, 'handler_set', None):
            log.setLevel(logging.INFO)
            sh = logging.StreamHandler()
            sh.setFormatter(formatter)
            log.addHandler(sh)

            fh = logging.FileHandler(self.LOG_FILENAME)
            fh.setLevel(logging.INFO)
            fh.setFormatter(formatter)
            log.addHandler(fh)

            log.setLevel(loglevel)
            log.handler_set = True
        return log
#enddef


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
    Modified by Chun Ly, 17 February 2018
     - Bug fix: Handle FITS files that are not multi-extension FITS
       (i.e., only one read)
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
        zhdr = fits.getheader(files0[ii], ext=0)

        extend = 'EXTEND' in zhdr.keys()
        if not extend:
            hdr = zhdr.copy()
        else:
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
    Extract FITS header comments. This function is to resolve issue with FITS
    block size, which prevents fits.Header.tofile() from working

    NOTE: This code is now obsolete

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
     - Moved up to work with read_template()
    Modified by Chun Ly, 22 February 2018
     - Updated documentation
    '''

    split0 = [str0.split(' / ')[-1].replace('\n','') for str0 in f0]

    comm0 = [split0[xx] if split0[xx] != f0[xx].replace('\n','') else ''
             for xx in range(len(split0))]
    return comm0
#enddef

def remove_padding(outfile, mylog=None):
    '''
    Removes space padding for FITS keyword value string. The padding is done
    automatically by astropy.io.fits.  This breaks mmirs-pipeline as the
    IDL code does not remove space paddings

    Parameters
    ----------
    c_hdr0 : astropy.io.fits.header.Header
      FITS header

    outfile : str
      Full path of filename to write updated FITS header file

    Returns
    -------
    None

    Notes
    -----
    Created by Chun Ly, 21 February 2018
     - Add mylog keyword input; Implement stdout and ASCII logging with mlog()
     - Bug fix: indentation mistake
     - Bug fix: CONTINUE cards are lost. Major modification: Need to modify
                ASCII file directly instead of FITS header object
    '''

    if type(mylog) == type(None): mylog = log

    f1 = open(outfile, 'r')
    str_hdr = f1.readlines()
    f1.close()

    i_fix = [xx for xx in range(len(str_hdr)) if '=' in str_hdr[xx]]

    for tt in range(len(i_fix)):
        split0 = str_hdr[i_fix[tt]].split('= ')
        if len(split0) > 1:
            right0 = split0[1]
            if (right0[:3] != "' '") and ('OGIP' not in right0) and \
               ('SIMPLE' not in split0[0]):
                lloc = right0.find("'")
                rloc = right0.rfind("'")
                if lloc != -1 and rloc != -1:
                    t_val   = right0[lloc:rloc+1]
                    r_t_val = t_val.replace(' ','')
                    right0  = r_t_val.ljust(len(t_val), ' ')
                    str_hdr[i_fix[tt]] = split0[0].ljust(8)+'= '+right0+'\n'
                else:
                    t_split = right0.split('/')
                    t_val   = t_split[0]
                    r_t_val = t_val.replace(' ','')
                    if '/' in right0:
                        right0  = r_t_val+' /'+t_split[1].replace("\n",'')
                    else:
                        right0  = r_t_val
                    str_hdr[i_fix[tt]] = split0[0].ljust(8)+'= '+right0+'\n'
    #endfor

    mylog.info('Updating : '+os.path.basename(outfile))
    f1 = open(outfile, 'w')
    f1.writelines(str_hdr)
    f1.close()
#enddef

def read_template(longslit=False, mos=False, mylog=None):
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

    Modified by Chun Ly, 19 February 2018
     - Add mylog keyword input; Implement stdout and ASCII logging with mlog()

    Modified by Chun Ly, 22 February 2018
     - Remove extraneous commented out code
    '''

    if type(mylog) == type(None): mylog = log

    if longslit == False and mos == False:
        mylog.warn('Must specify longslit or mos keyword!!!')
        mylog.warn('Exiting!')
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

    return temp_hdr

#enddef

def get_calib_files(name, tab0, mylog=None):
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
    Modified by Chun Ly, 20 February 2018
     - Add mylog keyword input; Implement stdout and ASCII logging with mlog()
    Modified by Chun Ly, 1 May 2018
     - Handle mos case when object is not the same as filename
       (i.e., mask vs mask_filter)
     - Handle case when comps are not available
     - Handle case when flats are not available
     - Handle mylog calls and calib_dict0 when comps/flats are not available
    Modified by Chun Ly, 25 May 2018
     - Change from underscore to colon separator
    '''

    if type(mylog) == type(None): mylog = log # + on 20/02/2018

    len0 = len(tab0)

    itype0 = tab0['imagetype']
    aper0  = tab0['aperture']
    filt0  = tab0['filter']
    disp0  = tab0['disperse']

    # + on 11/12/2017
    pi     = tab0['PI']
    propid = tab0['PropID']

    t_str  = name.split(':')
    t_obj  = t_str[0]
    t_ap   = t_str[1]
    t_filt = t_str[2]
    t_disp = t_str[3]

    i_obj = [ii for ii in range(len0) if t_obj in tab0['object'][ii]]

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
        t_txt = ",".join(tab0['filename'][i_dark])
        dark_str0.append(t_txt)
        mylog.info("List of science dark files for %.3fs : %s" % (etime, t_txt))

    ## COMPS
    # Mod on 11/12/2017
    i_comp = [ii for ii in range(len0) if
              (itype0[ii] == 'comp' and aper0[ii] == t_ap and
               filt0[ii] == t_filt and disp0[ii] == t_disp and \
               pi[ii] == t_pi and propid[ii] == t_propid)]

    if len(i_comp) != 0: # Mod 01/05/2018
        if len(i_comp) > 2:
            mylog.warn('Too many comps!!!') # Mod on 20/02/2018

        comp_str0  = ",".join(tab0['filename'][i_comp])
        comp_itime = tab0['exptime'][i_comp[0]]

        # Darks for comps
        id_comp = [ii for ii in range(len0) if
                   (itype0[ii] == 'dark' and tab0['exptime'][ii] == comp_itime)]
        comp_dark = ",".join(tab0['filename'][id_comp])
    else:
        mylog.warn('No comps found!!')

    ## FLATS
    # Mod on 11/12/2017
    i_flat = [ii for ii in range(len0) if
              (itype0[ii] == 'flat' and aper0[ii] == t_ap and \
               filt0[ii] == t_filt and disp0[ii] == t_disp and \
               pi[ii] == t_pi and propid[ii] == t_propid)]

    if len(i_flat) != 0: # Mod 01/05/2018
        if len(i_flat) > 2:
            mylog.warn('Too many flats!!!') # Mod on 20/02/2018

        flat_str0 = ",".join(tab0['filename'][i_flat])
        flat_itime = tab0['exptime'][i_flat[0]]

        # Darks for flats
        id_flat = [ii for ii in range(len0) if
                   (itype0[ii] == 'dark' and tab0['exptime'][ii] == flat_itime)]
        flat_dark = ",".join(tab0['filename'][id_flat])
    else:
        mylog.warn('No flats found!!')

    # Mod on 20/02/2018
    if len(i_comp) != 0:
        mylog.info("List of comp files : "+comp_str0)
        mylog.info("List of comp dark files : "+comp_dark)

    # Mod on 20/02/2018
    if len(i_flat) != 0:
        mylog.info("List of flat files : "+flat_str0)
        mylog.info("List of flat dark files : "+flat_dark)

    calib_dict0 = collections.OrderedDict()

    # Mod on 01/05/2018
    if len(i_comp) != 0:
        calib_dict0['comp_str']  = comp_str0
        calib_dict0['comp_dark'] = comp_dark # + on 11/12/2017
    else:
        calib_dict0['comp_str']  = ''
        calib_dict0['comp_dark'] = ''

    if len(i_flat) != 0:
        calib_dict0['flat_str']  = flat_str0
        calib_dict0['flat_dark'] = flat_dark # + on 11/12/2017
    else:
        calib_dict0['flat_str']  = ''
        calib_dict0['flat_dark'] = ''

    calib_dict0['dark_etime'] = dark_etime
    calib_dict0['dark_str']   = dark_str0

    return calib_dict0
#enddef

def handle_tellurics(tab0, object0, PropID, i_tell, obj_etime, tell_comb0,
                     idx, target_setup, mmirs_setup0, inter=False, mylog=None):
    '''
    Handle when multiple telluric datasets are available

    Parameters
    ----------
    tab0: astropy.table.table
      Astropy Table containing FITS header

    object0 : list
      List of strings of the source name (handle MOS case)

    PropID: str
      Proposal ID full name (for science exposures)

    i_tell : list or np.array
      Index of entries that are telluric stars

    obj_etime : list
      Telluric name and exposure time

    tell_comb0 : list
      List of strings that combines telluric name and exposure time

    idx : list or np.array
      Index of entries for a given target

    target_setup : str
      String that combines aperture, filter, and disperse for specific
      target

    mmirs_setup0 : list
      List of strings that combines aperture, filter, and disperse.
      This gets proper setup for telluric star identification

    inter : boolean
      For interactive telluric star selection.  Default: False

    Returns
    -------
    rev_tell_comb0 : list
      Improved telluric name and exposure time

    Notes
    -----
    Created by Chun Ly, 5 March 2018

    Modified by Chun Ly, 6 March 2018
     - Simplification improvements
     - Add idx input
     - Add mylog keyword input; Implement stdout and ASCII logging with mlog()
     - Simplification improvements (cont'd)
     - Add pass_score to tracking
     - Check if sci data is bracketed with telluric data
     - Return improved list of telluric data

    Modified by Chun Ly, 6 March 2018
     - Bug fix: log -> mylog changes

    Modified by Chun Ly, 1 May 2018
     - Fix to work with python3
     - Handle no telluric case for rev_tell_comb0
     - Fix typo for rev_tell_comb0

    Modified by Chun Ly, 4 June 2018
     - List all tellurics if more than one

    Modified by Chun Ly, 8 June 2018
     - Include inter keyword option
     - Add user prompts to identify telluric star
     - Bug fix: indexing issue

    Modified by Chun Ly, 9 June 2018
     - Add mmirs_setup0 and target_setup inputs
     - Require target_setup in telluric selection

    Modified by Chun Ly, 10 June 2018
     - Require target_setup in telluric selection (cont'd)
     - mylog.info spec setup info for telluric star list based on PropID
     - Clean up code
     - Include PropID info for interactive case
     - Bug fix: Correct index from user inputs
    '''

    if type(mylog) == type(None): mylog = log # + on 06/03/2018

    obj_etime = np.array(obj_etime) # Mod on 06/03/2018
    mmirs_setup0 = np.array(mmirs_setup0) # + on 09/06/2018

    pass_score = 0 # + on 06/03/2018

    # First distinguish by PropID
    i_prop = [xx for xx in range(len(tab0)) if tab0['PropID'][xx] == PropID]
    i_pid  = np.array(list(set(i_prop) & set(i_tell)))
    if len(i_pid) > 0:
        obj_etime_pid = list(set(obj_etime[i_pid]))
        setup_pid     = list(set(mmirs_setup0[i_pid])) # + on 10/06/2018
        if len(obj_etime_pid) == 1:
            mylog.info('Only one telluric dataset found using PropID !!!')
            mylog.info('## '+obj_etime_pid+' '+setup_pid) # + on 10/06/2018
            pass_score = 1 # + on 06/03/2018
            rev_tell_comb0 = list(obj_etime_pid) # + on 06/03/2018
        else:
            mylog.warn('More than one telluric dataset found using PropID : '+\
                     str(len(obj_etime_pid)))
            for tell_name_test in obj_etime_pid: # + on 04/06/2018
                tell_idx = [xx for xx in range(len(i_pid)) if
                            obj_etime[i_pid][xx] == tell_name_test]
                # + on 10/06/2018
                mylog.info('## '+tell_name_test+' '+mmirs_setup0[i_pid[tell_idx[0]]])
    else:
        mylog.warn('No tellurics found using PropID !!!')
        rev_tell_comb0 = list(tell_comb0) # + on 06/03/2018

    # Determine if telluric data bracket sci data | + on 06/03/2018
    if not pass_score:
        i_obj = [ii for ii in range(len(tab0)) if
                 (tab0['imagetype'][ii] == 'object' and tab0['aptype'][ii] != 'open')]
        tab0_nocalib = tab0[i_obj]

        obj_etime_nocalib = obj_etime[i_obj]

        mmirs_setup_nocalib = mmirs_setup0[i_obj] # + on 09/06/2018

        tell_idx_min = np.zeros(len(tell_comb0))
        tell_idx_max = np.zeros(len(tell_comb0))

        for tt in range(len(tell_comb0)):
            # Mod on 09/06/2018
            t_idx = [xx for xx in range(len(tab0_nocalib)) if
                     (obj_etime_nocalib[xx] == tell_comb0[tt] and
                      mmirs_setup_nocalib[xx] == target_setup)]
            #t_idx = np.array(i_obj)[t_idx]
            tell_idx_min[tt], tell_idx_max[tt] = min(t_idx), max(t_idx)

        idx_nocalib = [ii for ii in range(len(tab0_nocalib)) if
                       obj_etime_nocalib[ii] == obj_etime[idx[0]]]
        #idx_nocalib = np.array(i_obj)[idx_nocalib]
        sci_idx_min, sci_idx_max = min(idx_nocalib), max(idx_nocalib)

        tmp_tell_comb0 = [] # + on 06/03/2018

        # Mod on 07/06/2018
        if inter == False:
            # Check before
            bef0 = np.where(tell_idx_max - sci_idx_min == -1)[0]
            if len(bef0) == 1:
                mylog.info('Telluric data found before science data : '+\
                           tell_comb0[bef0[0]])
                tmp_tell_comb0.append(tell_comb0[bef0[0]])
                pass_score += 1
            else:
                mylog.info('NO telluric data before science data')

            # Check after
            aft0 = np.where(tell_idx_min - sci_idx_max == 1)[0]
            if len(aft0) == 1:
                mylog.info('Telluric data found after science data : '+tell_comb0[aft0[0]])
                tmp_tell_comb0.append(tell_comb0[aft0[0]]) # + on 06/03/2018
                pass_score += 1
            else:
                mylog.info('NO telluric data after science data')

            # + on 06/03/2018
            if len(bef0) == 1 or len(aft0) == 1:
                rev_tell_comb0 = tmp_tell_comb0
            if len(bef0) == 0 or len(aft0) == 0:
                rev_tell_comb0 = []
        else:
            print tab0[idx]

            # Check before
            bef0 = np.where(tell_idx_max - sci_idx_min < 0)[0]
            if len(bef0) > 0:
                mylog.info('Select telluric star BEFORE science observations : ')
                for bb in range(len(bef0)):
                    bb_bef = bef0[bb]
                    t_idx = [xx for xx in range(len(tab0_nocalib)) if
                             (obj_etime_nocalib[xx] == tell_comb0[bb_bef] and
                              mmirs_setup_nocalib[xx] == target_setup)]
                    tmpt    = tab0_nocalib[t_idx][0]
                    tmp_num = tab0_nocalib['seqno'][t_idx]
                    tell_str  = '(%i) %i-%i ' % (bb, min(tmp_num), max(tmp_num))
                    tell_str += '%s %s %s+%s %s' % (tell_comb0[bb_bef], tmpt['aperture'],
                                                    tmpt['filter'], tmpt['disperse'],
                                                    tmpt['PropID'])
                    mylog.info(tell_str)

                raw_bef = raw_input("Select from above telluric star to use : ")
                tmp_tell_comb0.append(tell_comb0[bef0[np.int(raw_bef)]])
                mylog.info("User selected : (%s) " % raw_bef)

            # Check after
            aft0 = np.where(tell_idx_min - sci_idx_max > 0)[0]
            if len(aft0) > 0:
                mylog.info('Select telluric star AFTER science observations : ')
                for bb in range(len(aft0)):
                    bb_aft = aft0[bb]
                    t_idx = [xx for xx in range(len(tab0_nocalib)) if
                             (obj_etime_nocalib[xx] == tell_comb0[bb_aft] and
                              mmirs_setup_nocalib[xx] == target_setup)]
                    tmpt    = tab0_nocalib[t_idx][0]
                    tmp_num = tab0_nocalib['seqno'][t_idx]
                    tell_str  = '(%i) %i-%i ' % (bb, min(tmp_num), max(tmp_num))
                    tell_str += '%s %s %s+%s %s' % (tell_comb0[bb_aft], tmpt['aperture'],
                                                    tmpt['filter'], tmpt['disperse'],
                                                    tmpt['PropID'])
                    mylog.info(tell_str)

                raw_aft = raw_input("Select from above telluric star to use : ")
                tmp_tell_comb0.append(tell_comb0[aft0[np.int(raw_aft)]])
                mylog.info("User selected : (%s) " % raw_aft)

            rev_tell_comb0 = tmp_tell_comb0
    #endif

    return rev_tell_comb0
#enddef

def get_tellurics(tab0, idx, comb0, object0, mmirs_setup0, inter=False, mylog=None):
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

    object0 : list
      List of strings that contains name of source. This is from
      organize_targets() and handles mask observations

    mmirs_setup0 : list
      List of strings that combines aperture, filter, and disperse.
      This handle MOS cases to get proper setup for telluric star identification

    inter : boolean
      For interactive telluric star selection.  Default: False

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
    Modified by Chun Ly, 18 February 2018
     - Handle mos case for tellurics with object0 input
     - Bug fix : obj -> object0
    Modified by Chun Ly, 20 February 2018
     - Add mylog keyword input; Implement stdout and ASCII logging with mlog()
    Modified by Chun Ly,  5 March 2018
     - Call handle_tellurics() function
    Modified by Chun Ly,  6 March 2018
     - Pass idx to handle_tellurics()
     - Get improved telluric list from handle_tellurics()
    Modified by Chun Ly,  7 March 2018
     - Pass mylog to handle_tellurics()
     - str_tell, tell_etime, str_dark and tell_stype definitions after
       handle_tellurics() call
    Modified by Chun Ly, 12 March 2018
     - Add BD telluric star in telluric list comprehension
     - Handle BD telluric star names for call to Simbad
    Modified by Chun Ly, 17 March 2018
     - mmirs-pipeline only accept a0v or g0v stars, changing this to work with
       templates
    Modified by Chun Ly, 1 May 2018
     - Turn off Simbad query if use_simbad=0
     - Add mmirs_setup0 input
     - Bug fix: mmirs_setup -> mmirs_setup0
    Modified by Chun Ly, 8 June 2018
     - Include inter keyword option
     - Pass inter keyword to handle_tellurics()
    Modified by Chun Ly, 9 June 2018
     - Pass target_setup and mmirs_setup0 to handle_tellurics()
     - Handle multi-spec setups for same telluric + integration time
    '''

    if type(mylog) == type(None): mylog = log # + on 20/02/2018

    etime = tab0['exptime'] # + on 28/01/2018

    target_setup = mmirs_setup0[idx[0]]

    # Mod on 12/03/2018
    i_tell = [xx for xx in range(len(object0)) if
              (('HD' in object0[xx] or 'HIP' in object0[xx] or 'BD_' in object0[xx]) and
               (mmirs_setup0[xx] == target_setup))]

    # Include exptime should data with multiple exptime for same target is taken
    # Mod on 28/01/2018
    obj_etime  = [a+'_'+str(b) for a,b in zip(object0, etime)]
    tell_comb0 = list(set(np.array(obj_etime)[i_tell]))

    n_tell = len(tell_comb0)

    if n_tell == 0:
        mylog.warn('No telluric data found!!!') # Mod on 20/02/2018

    if n_tell == 1:
        mylog.info('Only one telluric star is found!!!') # Mod on 20/02/2018
        mylog.info(object0[i_tell[0]]+' '+str(etime[i_tell[0]]))

    # + on 05/03/2018, Mod on 06/03/2018, 07/03/2018
    if n_tell > 1:
        PropID = tab0['PropID'][idx[0]]
        tell_comb0 = handle_tellurics(tab0, object0, PropID, i_tell, obj_etime,
                                      tell_comb0, idx, target_setup,
                                      mmirs_setup0, inter=inter, mylog=mylog)
        n_tell = len(tell_comb0)

    # Moved lower on 07/03/2018
    str_tell   = ['']  * n_tell
    tell_etime = [0.0] * n_tell
    str_dark   = ['']  * n_tell # + on 28/01/2018
    tell_stype = ['']  * n_tell # + on 29/01/2018

    if n_tell >= 1:
        for tt in range(n_tell):
            tmp = [xx for xx in range(len(object0)) if
                   (obj_etime[xx] == tell_comb0[tt] and
                    mmirs_setup0[xx] == target_setup)]

            # tell_time[tt] = tab0['dateobs'][tmp[0]]
            str_tell[tt]   = ",".join(tab0['filename'][tmp])
            tell_etime[tt] = etime[tmp[0]] # + on 28/01/2018

            # + on 28/01/2018
            i_dark = [xx for xx in range(len(object0)) if
                      (tab0['imagetype'][xx] == 'dark' and
                       etime[xx] == tell_etime[tt])]
            str_dark[tt] = ",".join(tab0['filename'][i_dark])

            # + on 29/01/2018; Mod on 12/03/2018
            tell_name = object0[tmp[0]]
            if tell_name[:2] == 'BD':
                tmp = tell_name.replace('BD_','BD+')
                tell_name = tmp[:5]+' '+tmp[5:]

            # Mod on 01/05/2018
            if use_simbad:
                t_simbad = Simbad.query_object(tell_name) #object0[tmp[0]])

                # + on 17/03/2018
                sp_type = t_simbad['SP_TYPE'][0].lower()
                if sp_type[0] == 'a' or sp_type[0] == 'g':
                    if sp_type[1:] != '0v':
                        mylog.warn(tell_comb0[tt]+' is '+sp_type+', not A0V or G0V!')
                        tell_stype[tt] = sp_type[0]+'0v'
                        mylog.warn('Setting : '+tell_comb0[tt]+' as '+tell_stype[tt])
                    else:
                        tell_stype[tt] = sp_type
                else:
                    mylog.warn(tell_comb0[tt]+' is '+sp_type+', NEITHER A or G star!')
                    tell_stype[tt] = ''
            else:
                tell_stype[tt] = 'a0v'
        #endfor

    # Mod on 28/01/2018, 29/01/2018
    tell_dict0 = {'name': str_tell, 'etime': tell_etime, 'dark': str_dark,
                  'stype': tell_stype}
    return tell_dict0
#enddef

def get_diff_images(tab0, idx, dither=None, mylog=None):
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
    Modified by Chun Ly, 20 February 2018
     - Add mylog keyword input; Implement stdout and ASCII logging with mlog()
    Modified by Chun Ly, 2 March 2018
     - Bug fix: Add option for only two exposures
     - Bug fix: Correctly handle two-exposure case
    Modified by Chun Ly, 6 April 2018
     - Round dither offsets to nearest integer pixel (0.2")
    Modified by Chun Ly, 9 April 2018
     - mylog call for when dither offsets are not integer pixels
    Modified by Chun Ly, 29 April 2018
     - Handle case when there is NO dithering
    '''

    if type(mylog) == type(None): mylog = log # + on 20/02/2018

    tab0 = tab0[idx]
    n_files = len(tab0)

    instel = tab0['instel']
    offset_val = list(set(instel))

    if dither == None:
        mylog.info('Determining dither sequence...') # Mod on 20/02/2018

        if len(offset_val) == 1: dither='None' # + on 29/04/2018

        if len(offset_val) == 4: dither="ABApBp"

        # Mod on 08/03/2018
        if len(offset_val) == 2:
            if len(instel) == 2:
                dither="ABAB"
            else:
                if instel[1] == instel[2]: dither="ABBA"
                if instel[1] != instel[2]: dither="ABAB"

    #Mod on 29/04/2018
    if dither != 'None':
        mylog.info('Dither sequence is : '+dither) # Mod on 20/02/2018

        i_off = [1, -1] * (n_files/2)
        if dither != 'ABBA':
            if n_files % 2 == 1: i_off.append(-1) # Odd number correction
        i_sky = np.arange(n_files)+np.array(i_off)

        im_dict = collections.OrderedDict()
        im_dict['sci']      = tab0['filename']
        im_dict['sci2']     = tab0['filename'][i_sky] # This is the sky frame

        # Mod on 06/04/2018, 09/04/2018
        base = 0.2 # pixel scale in arcsec for MMIRS
        t_instel = np.round_(base * np.round(np.float_(instel)/base), decimals=1)
        dither_fix = np.where(t_instel != instel)[0]
        if len(dither_fix) > 0:
            mylog.warn('Dither sequence is NOT integer pixels!')
            mylog.warn('Fixing for mmirs-pipeline!!!')
            instel = t_instel.copy()

        im_dict['dithpos']  = instel
        im_dict['dithpos2'] = instel[i_sky]
    else:
        mylog.warn('Dither sequence is : '+dither+' !!!')
        im_dict = dither
    return im_dict
#enddef

def generate_taskfile(hdr0, rawdir, w_dir, name, c_dict0, tell_dict0, tab0,
                      idx, dither=None, mylog=None):
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

    Modified by Chun Ly, 20 February 2018
     - Add mylog keyword input; Implement stdout and ASCII logging with mlog()
     - Pass mylog into get_diff_images()

    Modified by Chun Ly, 21 February 2018
     - Call remove_padding() to fix extraneous padding for FITS keyword value string
     - Simplify logging to exclude full outfile path
     - Change remove_padding() input (no longer providing FITS header)
    Modified by Chun Ly, 11 April 2018
     - Handle Kspec data: 'Kspec' -> 'K' for task files
    Modified by Chun Ly, 29 April 2018
     - Handle case when there is NO dithering
    Modified by Chun Ly, 5 May 2018
     - Telluric star processing once (with first taskfile)
    Modified by Chun Ly, 22 May 2018
     - Handle HK3 data: 'HK3' -> 'HK' for task files
    Modified by Chun Ly, 25 May 2018
     - Change from underscore to colon separator
     - Add uscore_name since name separator is colon
    Modified by Chun Ly, 27 May 2018
     - Handle case when comps or flats are not available (set to '')
    Modified by Chun Ly, 6 May 2018
     - Add inter keyword option
     - Remove inter keyword option
    Modified by Chun Ly, 11 June 2018
     - mmirs-pipeline taskfiles are for pairs, so only need half of the files
    '''

    if type(mylog) == type(None): mylog = log # + on 20/02/2018

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

    # + on 11/12/2017, Mod on 25/05/2018
    t_str = name.split(':')
    t_obj, t_ap, t_filt, t_disp  = t_str[0], t_str[1], t_str[2], t_str[3]

    uscore_name = name.replace(':','_') # + on 25/05/2018

    # + on 11/04/2018
    if t_filt == 'Kspec': t_filt = 'K'

    # + on 22/05/2018
    if t_filt == 'HK3': t_filt = 'HK'

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
        if val0[vv] != '': # Mod on 27/05/2018
            c_hdr0[col1[vv]] = val0[vv]
        else:
            mylog.info("Deleting : "+col1[vv])
            del c_hdr0[col1[vv]]

    col2 = ['SCI', 'SCI2', 'DITHPOS', 'DITHPOS2']

    # + on 24/01/2018
    im_dict = get_diff_images(tab0, idx, dither=dither, mylog=mylog)

    # Only write files if dithering is done | Mod on 29/04/2018
    if type(im_dict) != str:
        # Write ASCII taskfiles | + on 24/01/2018
        keys1 = c_hdr0.keys()

        for ii in range(0,len(im_dict['sci']),2):
            mylog.info('Writing taskfile for : '+im_dict['sci'][ii]) # Mod on 20/02/2018
            for t_key in col2:
                c_hdr0[t_key] = im_dict[t_key.lower()][ii]

            # Only perform telluric star processing once | + on 05/05/2018
            if ii == 0:
                c_hdr0['S08PROC'] = 1
            else:
                c_hdr0['S08PROC'] = 0

            c_hdr0['W_DIR'] = w_dir + format(ii/2+1, '02') + '/'

            outfile = w_dir+uscore_name+'_'+format(ii/2+1, '02')+'.txt' # Mod on 31/01/2018
            mylog.info('Writing : '+os.path.basename(outfile)) # Mod on 20/02/2018

            # Mod on 16/02/2018
            c_hdr0.tofile(outfile, sep='\n', padding=False, overwrite=True)
            # False padding to avoid IOError: Header size (5725) is not a multiple
            # of block size (2880)

            # Fix extraneous padding for FITS keyword value string | + on 21/02/2018
            remove_padding(outfile, mylog=mylog)
        #endfor
    else:
        log.warn('Taskfiles are not generated, NO dithering - not supported!')

    return c_hdr0
#enddef

def organize_targets(tab0, mylog=None):
    '''
    Use FITS header information to organize targets based on name, aperture,
    filter, and disperse

    Parameters
    ----------
    tab0: astropy.table.table
      Astropy Table containing FITS header info

    Returns
    -------
    comb0 : list
      List of strings that combines name, aperture, filter, and disperse

    obj_comb0 : list
      A unique list of strings that combines name, aperture, filter, and disperse

    object0 : list
      List of strings of the source name (handle MOS case)

    mmirs_setup0 : list
      List of strings that combines aperture, filter, and disperse.
      This handle MOS cases to get proper setup for telluric star identification

    Notes
    -----
    Created by Chun Ly, 28 November 2017
     - Get unique lists of combinations
    Modified by Chun Ly, 17 February 2018
     - Ignore imaging data
     - Bug fix: Typo with obj list comprehension
     - Handle mos data and telluric data for mos
    Modified by Chun Ly, 18 February 2018
     - Return object name as object0 list
    Modified by Chun Ly, 20 February 2018
     - Add mylog keyword input; Implement stdout and ASCII logging with mlog()
    Modified by Chun Ly, 22 February 2018
     - Change print statement for combination to mylog call
    Modified by Chun Ly, 1 May 2018
     - Determine number of spec for each combination, ignore those with single spec
     - Define and return mmirs_setup0
    Modified by Chun Ly, 25 May 2018
     - Change from underscore to colon separator
     - log aesthetics: use underscore name instead of colon-separated name
    '''

    if type(mylog) == type(None): mylog = log # + on 20/02/2018

    len0 = len(tab0)

    itype  = tab0['imagetype']
    aptype = tab0['aptype'] # + on 17/02/2018

    # Mod on 17/02/2018
    obj     = [ii for ii in range(len0) if
               (itype[ii] == 'object' and aptype[ii] != 'open')]
    comp    = [ii for ii in range(len0) if itype[ii] == 'comp']
    flat    = [ii for ii in range(len0) if itype[ii] == 'flat']

    comb0        = ['N/A'] * len0
    object0      = ['N/A'] * len0 # + on 18/02/2018
    mmirs_setup0 = ['N/A'] * len0 # + on 01/05/2018

    for kk in [obj, comp, flat]:
        for ii in kk:
            tab0_o = tab0[ii]

            # Mod on 17/02/2018
            if tab0_o['aptype'] == 'mos':
                t_name = tab0_o['filename'].split('_')[0]
            else:
                t_name = tab0_o['object']

            object0[ii] = t_name # + on 18/02/2018
            comb0[ii] = t_name + ':' + tab0_o['aperture'] + ':' + \
                        tab0_o['filter'] + ':' + tab0_o['disperse']
            # + on 01/05/2018
            mmirs_setup0[ii] = tab0_o['aperture'] + ':' + \
                               tab0_o['filter'] + ':' + tab0_o['disperse']

    obj_comb0 = list(set(np.array(comb0)[obj]))

    n_obj_comb0 = len(obj_comb0)
    mylog.info('Total number of combinations found : '+str(n_obj_comb0))

    # Mod on 01/05/2018
    cnt_comb0 = np.zeros(n_obj_comb0)
    for oo in range(n_obj_comb0):
        oo_idx = [xx for xx in range(len(obj)) if
                  (np.array(comb0)[obj[xx]] == obj_comb0[oo])]
        cnt_comb0[oo] = len(oo_idx)

        mylog.info('## '+obj_comb0[oo].replace(':','_')+' N=%i' % cnt_comb0[oo])
    #endfor

    # Exclude those with single spectra | + on 01/05/2018
    single_data = np.where(cnt_comb0 == 1)[0]
    if len(single_data) > 0:
        for jj in single_data:
            mylog.info('The following will be excluded : '+\
                       obj_comb0[jj].replace(':','_'))

    obj_comb0 = np.delete(obj_comb0, single_data).tolist()

    return comb0, obj_comb0, object0, mmirs_setup0 # Mod on 18/02/2018, 01/05/2018
#enddef

def create(rawdir, w_dir='', dither=None, bright=False, extract=False,
           inter=False, silent=False, verbose=True, debug=False):

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

    extract: boolean
      Indicate whether to turn on 1-D extraction. Default: False

    inter : boolean
      For interactive telluric star selection.  Default: False

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
     - Write tab0 when no ASCII file is present
     - Bug fix: With using LS or MOS template
     - Bug fix: With using LS or MOS template (cont'd)
    Modified by Chun Ly, 18 February 2018
     - Bug fix: taskfiles were placed in incorrect directory when dealing with
                multiple targets. Fixed w_dir definition
     - Handle path when w_dir keyword input is given
     - Get proper source name to handle mask observations, object0
    Modified by Chun Ly, 19 February 2018
     - Implement stdout and ASCII logging with mlog()
    Modified by Chun Ly, 20 February 2018
     - Pass mylog into organize_targets, read_template, get_calib_files, get_tellurics
     - Rename sci_files np.array to idl_files for simplicity
     - Add all necessary darks to idl_files for IDL_input.lis
     - Bug fix with c_dict0 use
    Modified by Chun Ly, 21 February 2018
     - Bug fix: Appears that mmirs_pipeline does use RAW_DIR input
    Modified by Chun Ly, 22 February 2018
     - Remove call to chun_codes.systime()
    Modified by Chun Ly, 23 February 2018
     - Bug fix: missing single quote
     - Add check for calib_MMIRS and symlink command
    Modified by Chun Ly, 1 March 2018
     - Bug fix: typo in if statement
    Modified by Chun Ly, 7 March 2018
     - Add debug keyword option for code testing without having FITS files
    Modified by Chun Ly, 1 May 2018
     - Get mmirs_setup0 from organize_targets(), pass to get_tellurics()
     - Modify IDL script file for w_dir
     - Minor bug fix for IDL script file
    Modified by Chun Ly, 4 May 2018
     - Add extract keyword input, Update S07PROC keyword for hdr0
    Modified by Chun Ly, 21 May 2018
     - Handle cases without comps, flats, and darks
    Modified by Chun Ly, 23 May 2018
     - Change location of .lis and .idl files (rawdir to w_dir)
     - Change preproc location to w_dir in generate_taskfile() call
     - Include outdir in call to mmirs_pipeline_nonlin_script.pro
    Modified by Chun Ly, 24 May 2018
     - Change outdir to w_dir in call to mmirs_pipeline_nonlin_script.pro
     - Change calib_MMIRS symlink to w_dir path
     - Save mmirs_pipeline_taskfile log to w_dir
     - Save obs_summary.tbl to w_dir
    Modified by Chun Ly, 25 May 2018
     - Add uscore_name since name separator is colon
     - Use name instead of uscore_name for generate_taskfile
     - log aesthetics: use underscore name instead of colon-separated name
    Modified by Chun Ly, 26 May 2018
     - Copy msk files over to preproc folder (fields case)
    Modified by Chun Ly,  5 June 2018
     - Add on_error call for idl scripts for mmirs_pipeline
    Modified by Chun Ly,  7 June 2018
     - Add inter keyword option
    Modified by Chun Ly, 8 June 2018
     - Pass inter keyword to get_tellurics()
    '''

    mylog = mlog(rawdir)._get_logger() # + on 19/02/2018

    mylog.info('Begin create ! ') # Mod on 19/02/2018

    if rawdir[-1] != '/': rawdir = rawdir + '/'

    files0 = glob.glob(rawdir+'*.????.fits*')
    files0 = [t_file for t_file in files0 if 'align_' not in t_file]
    files0 = [t_file for t_file in files0 if '_stack' not in t_file]

    n_files0 = len(files0)

    # Mod on 19/02/2018
    mylog.info('Number of FITS files found : '+str(n_files0))

    tab0_outfile = rawdir + 'obs_summary.tbl' # Moved up on 15/02/2018

    # + on 07/03/2018
    if n_files0 == 0 and debug:
        tmp      = asc.read(tab0_outfile, format='fixed_width_two_line')
        files0   = [files+'.fits' for files in tmp['filename']]
        n_files0 = len(files0)

    # Use ASCII catalog if available and size matches | + on 15/02/2018
    if exists(tab0_outfile):
        mylog.info('Reading : '+tab0_outfile) # Mod on 19/02/2018
        tab0 = asc.read(tab0_outfile, format='fixed_width_two_line')
        if n_files0 != len(tab0):
            # Get header information
            tab0 = get_header_info(files0)

            # + on 30/11/2017. Mod on 19/02/2018
            if exists(tab0_outfile):
                mylog.info('Overwriting : '+tab0_outfile)
            else:
                mylog.info('Writing : '+tab0_outfile)

            tab0.write(tab0_outfile, format='ascii.fixed_width_two_line',
                       overwrite=True)
        else:
            mylog.info('Using existing tab0_outfile') # Mod on 19/02/2018
    else: # + on 16/02/2018
        tab0 = get_header_info(files0)

        tab0.write(tab0_outfile, format='ascii.fixed_width_two_line',
                   overwrite=True)

    tab0.pprint(max_lines=-1, max_width=-1)

    # Mod on 18/02/2018, 20/02/2018
    comb0, obj_comb0, object0, mmirs_setup0 = organize_targets(tab0, mylog=mylog)

    # Get default FITS template headers
    # + on 30/11/2017, Mod on 18/12/2017, 24/01/2018, 17/02/2018
    LS_hdr0  = read_template(longslit=True, mylog=mylog)
    mos_hdr0 = read_template(mos=True, mylog=mylog)

    # Create task files | + on 30/11/2017
    for name in obj_comb0:
        if 'HD' not in name and 'HIP' not in name and 'BD' not in name:
            uscore_name = name.replace(':','_')

            mylog.info('Working on : '+uscore_name) # Mod on 19/02/2018

            idx   = [ii for ii in range(len(comb0)) if comb0[ii] == name]
            n_idx = len(idx)

            # Mod on 30/11/2017, 24/01/2018, 17/02/2018
            aptype = tab0['aptype'][idx[0]]
            if 'longslit' in aptype:
                hdr0  = LS_hdr0.copy()
            if 'mos' in aptype:
                hdr0  = mos_hdr0.copy()

            # Mod on 21/02/2018
            hdr0['RAW_DIR'] = rawdir #orig_dir

            # + on 25/01/2018. Mod on 19/02/2018
            if bright:
                mylog.info('Bright source flag set!')
                hdr0['BRIGHT'] = 1
            else:
                mylog.info('No bright source indicated!')
                hdr0['BRIGHT'] = 0

            # + on 04/05/2018
            if extract:
                mylog.info('Extract flag set!')
                hdr0['S07PROC'] = 1
            else:
                mylog.info('Will not perform 1-D extraction!')
                hdr0['S07PROC'] = 0

            # on 08/12/2017
            c_dict0 = get_calib_files(name, tab0, mylog=mylog)

            # Mod on 28/01/2018, 18/02/2018, 20/02/2018, 08/06/2018
            tell_dict0 = get_tellurics(tab0, idx, comb0, object0, mmirs_setup0,
                                       inter=inter, mylog=mylog)

            # Mod on 31/01/2018, 18/02/2018
            if w_dir == '':
                w_dir_tmp = rawdir + 'reduced/'+uscore_name+'/'
            else:
                if w_dir[-1] != '/': w_dir = w_dir + '/'
                w_dir_tmp = w_dir + uscore_name + '/'

            if not exists(w_dir_tmp): # Mod on 18/02/2018
                commands.getoutput('mkdir -p '+w_dir_tmp)

            # Generate file containing science frames and darks to run
            # MMIRS IDL preprocessing
            # + on 12/02/2018
            idl_files  = tab0['filename'][idx].data

            # + on 12/02/2018, Mod on 21/05/2018
            if c_dict0['comp_str'] != '':
                t_comps   = c_dict0['comp_str'].split(',')
                idl_files = np.append(idl_files, t_comps)

            if c_dict0['flat_str'] != '':
                t_flats   = c_dict0['flat_str'].split(',')
                idl_files = np.append(idl_files, t_flats)

            # Add darks | + on 20/02/2018, Mod on 21/05/2018
            if c_dict0['dark_str'][0] != '':
                idl_files = np.append(idl_files, c_dict0['dark_str'][0].split(','))
            if c_dict0['comp_dark'] != '':
                idl_files = np.append(idl_files, c_dict0['comp_dark'].split(','))
            if c_dict0['flat_dark'] != '':
                idl_files = np.append(idl_files, c_dict0['flat_dark'].split(','))

            # Get telluric files | + on 12/02/2018, Mod on 21/05/2018
            for tt in range(len(tell_dict0['name'])):
                t_tell    = tell_dict0['name'][tt].split(',')
                t_dark = tell_dict0['dark'][tt].split(',') # + on 20/02/2018
                idl_files = np.append(idl_files, t_tell)
                idl_files = np.append(idl_files, t_dark) # + on 20/02/2018

            idl_input_file = w_dir + 'IDL_input_'+uscore_name+'.lis' # Mod on 23/05/2018
            mylog.info('Writing : '+idl_input_file) # Mod on 19/02/2018
            asc.write([idl_files], idl_input_file, format='no_header',
                      overwrite=True)

            # + on 11/12/2017, Mod on 28/01/2018, 18/02/2018
            temp1 = generate_taskfile(hdr0, w_dir, w_dir_tmp, name,
                                      c_dict0, tell_dict0, tab0, idx,
                                      dither=dither, mylog=mylog)

            # Write IDL script for main mmirs_pipeline | + on 17/02/2018
            main_script_outfile = w_dir + 'run_mmirs_pipeline_'+uscore_name+'.idl' # Mod on 23/05/2018
            if not exists(main_script_outfile):
                mylog.info('Writing : '+main_script_outfile)

                f1 = open(main_script_outfile, 'w')

                # Mod on 01/05/2018
                if w_dir == '':
                    str0 = [".run run_pipeline\n\n", "on_error, 0\n"
                            "run_pipeline, 'reduced/%s'\n\n" % uscore_name, "exit\n"] # Mod on 23/02/2018
                else:
                    str0 = [".run run_pipeline\n\n",
                            "run_pipeline, '%s'\n\n" % w_dir_tmp, "exit\n"] # Mod on 23/02/2018

                f1.writelines(str0)
                f1.close()
            else:
                mylog.warn('File exists! Will not overwrite : '+main_script_outfile)
        #endif
    #endfor

    # Write IDL script for non-linear pre-processing with
    # mmirs_pipeline_nonlin_script.pro | + on 17/02/2018
    script_outfile = w_dir + 'run_mmirs_pipeline_nonlin_script.idl' # Mod on 23/05/2018
    if not exists(script_outfile):
        mylog.info('Writing : '+script_outfile) # Mod on 19/02/2018

        f0 = open(script_outfile, 'w')
        str0 = [".run mmirs_pipeline_nonlin_script\n\n",
                "mmirs_pipeline_nonlin_script, '%s', w_dir='%s', compress='.gz', /verbose\n\n"% \
                (rawdir, w_dir), "exit\n"]
        f0.writelines(str0)
        f0.close()
    else:
        mylog.warn('File exists! Will not overwrite : '+script_outfile)

    # Check if calib_MMIRS exists
    dir_calib = w_dir+'calib_MMIRS'
    if not exists(dir_calib):
        mylog.warn('calib_MMIRS does NOT exists in w_dir !!!')

        m_path = raw_input('Enter path to main folder of IDL mmirs-pipeline : ')
        if m_path[-1] != '/': m_path += '/'
        c_path = m_path + 'pipeline/calib_MMIRS'
        cmd0   = 'ln -s '+c_path+' '+dir_calib
        mylog.info('symlink command : '+cmd0)
        commands.getoutput(cmd0)
        if not exists(dir_calib):
            mylog.warn('Path is not correct for calib_MMIRS !!!')
            mylog.warn('User need to fix manually !!!')

    # Copy msk files over to preproc folder | + on 26/05/2018
    pproc_dir = w_dir+'preproc'
    mylog.info("Copying files over to : " + pproc_dir)
    use_dir = ''
    if gethostname() == 'fields': # Check if it's fields
        m_files = glob.glob(rawdir+'*.msk')
        if len(m_files) > 0:
            use_dir = rawdir
        else:
            ccd_dir = rawdir.replace('/data/crunch','/data/ccd')
            mylog.info("Trying : "+ccd_dir)
            m_files = glob.glob(ccd_dir+'*.msk')
            if len(m_files) > 0:
                use_dir = ccd_dir
            else:
                archive_dir = rawdir.replace('/data/crunch','/data/archive')
                mylog.info("Trying : "+archive_dir)
                m_files = glob.glob(archive_dir+'*.msk')
                if len(m_files) > 0:
                    use_dir = archive_dir
            #endelse
        #endelse

        if use_dir != '':
            if not exists(pproc_dir):
                commands.getoutput('mkdir -p '+pproc_dir)

            cmd1a = 'cp -a %s*.msk %s' % (use_dir, pproc_dir)
            mylog.info(cmd1a)
            commands.getoutput(cmd1a)

            cmd1b = 'cp -a %s*_msk.fits %s' % (use_dir, pproc_dir)
            mylog.info(cmd1b)
            commands.getoutput(cmd1b)
        #endif
    #endif

    mylog.info('End create ! ') # Mod on 19/02/2018

    # Save mmirs_pipeline_taskfile log to w_dir | + on 24/05/2018
    cmd1 = 'cp -a %s %s ' % (rawdir+'mmirs_pipeline_taskfile.log', w_dir)
    commands.getoutput(cmd1)

    # Save obs_summary.tbl to w_dir | + on 24/05/2018
    cmd2 = 'cp -a %s %s ' % (tab0_outfile, w_dir)
    commands.getoutput(cmd2)

#enddef

