"""
create_mmirs_obs_summary_table
====

Creates a table that summarizes observations from MMIRS. To be executed in a
folder that is organized by UTC date.
"""

import sys, os

from chun_codes import systime

from os.path import exists
from astropy.io import ascii as asc

import numpy as np

from glob import glob

from astropy.table import Table
from astropy import log

def get_telluric_info(QA_tab, t_idx, t_names):
    ## + on 25/04/2017

    if len(t_idx) == 1:
        telstar   = QA_tab['object'][t_idx][0]
        t_exptime = QA_tab['exptime'][t_idx][0]
        telset    = str(len(t_idx))+'x'+str(t_exptime)+'s'

        AM0   = QA_tab['airmass'][t_idx]
        telAM = '%.3f-%.3f' % (np.min(AM0),np.max(AM0))
    else:
        telset = []
        telAM  = []
        for nn in range(len(t_names)):
            idx2      = [xx for xx in range(len(QA_tab)) if
                         (t_names[nn] in QA_tab['object'][xx])]
            nidx      = list(set(t_idx) & set(idx2))
            t_exptime = QA_tab['exptime'][nidx][0]
            telset.append(str(len(nidx))+'x'+str(t_exptime)+'s')

            AM0   = QA_tab['airmass'][nidx]
            telAM.append('%.3f-%.3f' % (np.min(AM0),np.max(AM0)))

        telstar = ','.join(t_names)
        telset  = ','.join(telset)
        telAM   = ','.join(telAM)
    return telstar, telset, telAM
#enddef

def main(path0, outfile=None, silent=False, verbose=True):
    '''
    Generate ASCII file summarizing observations

    Parameters
    ----------
    path0 : str
      Parent path for all files (above UTC date directory)

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 24 January 2019
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    Targets, ObsDate, ObsSet = [], [], []
    TotalTime, gratwave, Airmass = [], [], []
    TellStar, TellSet, TellAM = [], [], []

    dirs = glob(path0+'????.????')
    
    for dir0 in dirs:
        cnt = 0
        obs_file = dir0+'/obs_summary.tbl'
        if not exists(obs_file):
            log.warn('## File not found! '+obs_file)

            ObsDate.append('N/A')
            gratwave.append('N/A')
            
            ObsSet.append('N/A')
            TotalTime.append('N/A')
            Airmass.append('N/A')
            
            TellStar.append('N/A')
            TellSet.append('N/A')
            TellAM.append('N/A')
        else:
            obs_tab  = asc.read(obs_file, format='fixed_width_two_line')

            obj0 = obs_tab['object'].data
            # All science targets

            idx = [ii for ii in range(len(obs_tab)) if
                   (obs_tab['imagetype'][ii] == 'object' and
                    ('HIP' not in obj0[ii] and 'HD' not in obj0[ii] and
                     'BD_' not in obj0[ii]))]
            idx = np.array(idx)

            target_names = list(set(obj0[idx]))
            print("Targets : "+", ".join(target_names))

            for target in target_names:
                s_idx = np.where(obj0[idx] == target)[0]
                s_idx = idx[s_idx]

                tab_ref = obs_tab[s_idx][0]
                t_date  = tab_ref['dateobs'].split('T')[0]
                exptime = tab_ref['exptime']

                Targets.append(target)
                ObsDate.append(t_date)
                gratwave.append(tab_ref['filter']+'_'+tab_ref['disperse'])

                ObsSet.append(str(len(s_idx))+'x'+str(exptime)+'s')
                TotalTime.append('%.2f' % (len(s_idx)*exptime/60.0))

                AM0 = obs_tab['airmass'][s_idx]
                Airmass.append('%.3f-%.3f' % (np.min(AM0),np.max(AM0)))

                #t_idx = [ii for ii in range(len(obs_tab)) if
                #         (obs_tab['imagetype'][ii] == 'object' and
                #          ('HIP' in obj0[ii] or 'HD' not in obj0[ii]))]
                #t_names = list(set(obj0[t_idx]))

            #telstar, telset, telAM = get_telluric_info(obs_tab, t_idx, t_names)
            #TellStar.append(telstar)
            #TellSet.append(telset)
            #TellAM.append(telAM)
            cnt += 1
        #endelse
    #endfor

    arr0   = [Targets, ObsDate, ObsSet, TotalTime, gratwave, Airmass]
    # TellStar, TellSet, TellAM]
    names0 = ('Name', 'UT_Date', 'Sequence', 'Int_Time', 'Grating_Wave',
              'Airmass', 'Telluric_Star', 'Telluric_Seq', 'Telluric_AM')
    tab0 = Table(arr0, names=names0)

    print tab0

    if outfile == None: outfile = path0+'obs_summary.txt'

    # Mod on 06/05/2017
    if silent == False:
        stat0 = 'Overwriting : ' if exists(outfile) else 'Writing : '
        log.info(stat0+outfile)
    asc.write(tab0, output=outfile, format='fixed_width_two_line', overwrite=True)

    if silent == False: log.info('### End main : '+systime())
#enddef
