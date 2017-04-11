"""
mmirs
=====

A set of codes for MMIRS observing
"""

import sys, os

#from os.path import exists
#from astropy.io import ascii as asc
#from astropy.io import fits

import numpy as np

# import matplotlib.pyplot as plt
# import glob

#from astropy.table import Table
from astropy import log
import astropy.coordinates as coords
import astropy.units as u
from astropy.time import Time

def systime():
  import time
  return time.strftime("%d_%b_%Y_%H:%M:%S", time.localtime())
#enddef

def compute_mmirs_offsets(ra, dec, pm_ra, pm_dec, date='', epoch=0,
                          silent=False, verbose=True):

    '''
    This code computes offsets between two coordinates for a given epoch.
    It factors in proper motion. The coordinate should be entered as:
     (offset, target) for ra, dec, pm_ra, pm_dec, epoch

    Parameters
    ----------
    ra : list or array
     J2000 RA in hours. Can either be a string or a decimal value

    dec : list or array
     J2000 Dec in degrees. Can either be a string or a decimal value

    pm_ra : list or array
     RA proper motion in units of mas/yr. The convention is:
     pm_ra = mu_RA * cos(Dec)

    pm_dec : list or array
     Dec proper motion in units of mas/yr.

    date : string
     Date of observations in ISO format 'YYYY-MM-DD'

    epoch : float
     The epoch to compute offsets for. Optional. Will use 

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
    Offset values

    Notes
    -----
    Created by Chun Ly, 10 April 2017
    '''
    
    if silent == False: log.info('### Begin compute_mmirs_offsets : '+systime())

    if date != '':
        epoch = Time(date).decimalyear
    dtime = epoch - 2000.0

    # J2000 coordinates
    c0 = coords.SkyCoord(ra[0], dec[0], frame='icrs',
                         unit=(u.hr, u.deg)) # offset star
    ct = coords.SkyCoord(ra[1], dec[1], frame='icrs',
                         unit=(u.hr, u.deg)) # target

    log.info('J2000 coord, offset star : '+ra[0]+' '+dec[0])
    log.info('J2000 coord, target      : '+ra[1]+' '+dec[1])

    # offset star
    pmRA = pm_ra[0] * dtime/np.cos(np.radians(c0.dec.value))
    ra_off  = coords.Angle(pmRA, unit=u.marcsec)
    dec_off = coords.Angle(pm_dec[0] * dtime, unit=u.marcsec)
    new_c0  = coords.SkyCoord(c0.ra + ra_off, c0.dec + dec_off, 'icrs')
    
    # target
    pmRA    = pm_ra[1] * dtime/np.cos(np.radians(ct.dec.value))
    ra_off  = coords.Angle(pmRA, unit=u.marcsec)
    dec_off = coords.Angle(pm_dec[1] * dtime, unit=u.marcsec)
    new_ct  = coords.SkyCoord(ct.ra + ra_off, ct.dec + dec_off, 'icrs')

    new_c0_str = new_c0.to_string('hmsdms')
    new_ct_str = new_ct.to_string('hmsdms')
    log.info('Epoch='+str(epoch)+' coord, offset star : '+new_c0_str)
    log.info('Epoch='+str(epoch)+' coord, target      : '+new_ct_str)
    
    # For J2000
    dra0  = (ct.ra.deg - c0.ra.deg) * 3600.0 # * np.cos(new_ct.dec.radian)
    ddec0 = (ct.dec.deg - c0.dec.deg) * 3600.0

    PA1  = c0.position_angle(ct).degree - 180.0 # + => East of North

    print ''
    log.info('IF PROPER MOTION WAS NOT CONSIDERED!!!')
    log.info('This offset in RA  for alignboxmmirs : %.3f' % dra0)
    log.info('This offset in Dec for alignboxmmirs : %.3f' % ddec0)
    log.info('PA for longslit alignment : %.3f, %.3f, or %.3f' % (PA1, PA1-180,
                                                                  PA1+180))

    dra0  = (new_ct.ra.deg - new_c0.ra.deg) * 3600.0 # * np.cos(new_ct.dec.radian)
    ddec0 = (new_ct.dec.deg - new_c0.dec.deg) * 3600.0

    PA  = new_c0.position_angle(new_ct).degree - 180.0 # + => East of North

    print ''; print ''
    log.info('IF PROPER MOTION IS CONSIDERED!!!')
    log.info('This offset in RA  for alignboxmmirs : %.3f' % dra0)
    log.info('This offset in Dec for alignboxmmirs : %.3f' % ddec0)
    log.info('PA for longslit alignment : %.3f, %.3f, or %.3f' % (PA, PA-180, PA+180))
    
    if silent == False: log.info('### End compute_mmirs_offsets : '+systime())
#enddef

