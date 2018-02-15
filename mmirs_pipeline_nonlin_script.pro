PRO mmirs_pipeline_nonlin_script, rawdir, first=first, linear=linear, $
                                  keepfirst=keepfirst, verbose=verbose, $
                                  debug=debug, biasframe=biasframe, $
                                  badamp=badamp, crosstalk=crosstalk, $
                                  compress=compress, tmpdir=tmpdir, clean=clean
;+
; NAME:
;       MMIRS_PIPELINE_NONLIN_SCRIPT
;
; PURPOSE:
;       Runs mmirs_pipeline's mmfixen_nonlin
;
; CALLING SEQUENCE:
;       mmirs_pipeline_nonlin_script, rawdir
;
; INPUTS:
;       rawdir - Directory path to raw files. Must end with '/'
;
; OPTIONAL KEYWORD INPUT:
;       None.
;
; OUTPUTS:
;
; OPTIONAL OUTPUT KEYWORD:
;       None.
;
; PROCEDURES USED:
;
; NOTES:
;
; REVISON HISTORY:
;       Created by Chun Ly, 14 February 2018
;
;-

  suffix='.fix.fits'

  IDL_infile = rawdir + 'IDL_input.lis'
  print, '### Reading : '+IDL_infile + ' | '+systime()
  READCOL, IDL_infile, files0, format='A,X'

  outdir = rawdir+'preproc/'
  if not file_test(outdir) then begin
     spawn, 'mkdir '+outdir
  endif

  for ii=0L,N_elements(files0)-1 do begin
     outfile = outdir + files0[ii] + suffix
     mmfixen_nonlin, rawdir+files0[ii]+'.fits', outfile, first=first, $
                     linear=linear, keepfirst=keepfirst, verbose=verbose, $
                     debug=debug, biasframe=biasframe, badamp=badamp, $
                     crosstalk=crosstalk, compress=compress, tmpdir=tmpdir, $
                     clean=clean  
  endfor

END
