PRO mmirs_pipeline_nonlin_script, rawdir, w_dir=w_dir, first=first, $
                                  linear=linear, keepfirst=keepfirst, $
                                  verbose=verbose, debug=debug, $
                                  biasframe=biasframe, badamp=badamp, $
                                  crosstalk=crosstalk, compress=compress, $
                                  tmpdir=tmpdir, clean=clean
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
;       Modified by Chun Ly, 17 February 2018
;        - Get all IDL_input.lis files and run mmfixend_nonlin for
;          each set
;
;       Modified by Chun Ly, 20 February 2018
;        - Bug fix: 'IDL_infiles' to 'IDL_files'
;        - Bug fix: READCOL format change (one column)
;
;       Modified by Chun Ly, 21 February 2018
;        - Check if file exists before running mmfixend_nonlin
;
;       Modified by Chun Ly, 23 May 2018
;        - Allow outdir keyword input option
;
;       Modified by Chun Ly, 24 May 2018
;        - Change outdir keyword input to w_dir
;-

  ; Moved up on 17/02/2018
  if not keyword_set(w_dir) then begin
     outdir = rawdir+'preproc/'
  endif else $
     outdir = w_dir+'preproc/'

  if not file_test(outdir) then $
     spawn, 'mkdir '+outdir

  suffix='.fix.fits'

  ; Get all IDL_input*lis files
  if not keyword_set(w_dir) then begin
     spawn, 'ls '+rawdir+'IDL_input*.lis', IDL_files
  endif else $
     spawn, 'ls '+w_dir+'IDL_input*.lis', IDL_files

  n_IDL_files= n_elements(IDL_files)

  for ff=0,n_IDL_files-1 do begin
     print, '### Reading : '+IDL_files[ff] + ' | '+systime()
     READCOL, IDL_files[ff], files0, format='A'

     for ii=0L,N_elements(files0)-1 do begin
        outfile = outdir + files0[ii] + suffix
        if not file_test(outfile + compress) then begin
           mmfixen_nonlin, rawdir+files0[ii]+'.fits', outfile, first=first, $
                           linear=linear, keepfirst=keepfirst, verbose=verbose, $
                           debug=debug, biasframe=biasframe, badamp=badamp, $
                           crosstalk=crosstalk, compress=compress, tmpdir=tmpdir, $
                           clean=clean
        endif else print, '## File exists! | ' + outfile + compress + ' '+systime()
     endfor
  endfor

END
