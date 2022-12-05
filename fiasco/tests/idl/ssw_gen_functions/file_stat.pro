function file_stat, files, exist=exist, size=size, fstat_all=fstat_all
;+
; Project     :	SOHO - SSW
;
; Name        :	FILE_STAT()
;
; Purpose     :	Vector version of FSTAT
;
; Category    :	Utility, Operating_system
;
; Explanation :	Vector version of FSTAT
;
; Syntax      :	Result = FILE_STAT( FILES )
;
; Examples    :	
;
; Inputs      :	FILES	= List of files to check.
;
; Opt. Inputs :	None.
;
; Outputs     :	None.
;
; Opt. Outputs:	None.
;
; Keywords    :	EXIST	= If set, then the returned values are whether the
;			  files exist or not.  This is the default behavior.
;		SIZE	= If set, then the returned values are the sizes of the
;			  files.
;               FSTAT_ALL = If set, return vector of 'fstat' structures.  This
;                         keyword is only compatible with files less than about
;                         2 GB each in size.
;
; Calls       :	DATA_CHK
;
; Common      :	None.
;
; Restrictions:	None.
;
; Side effects:	None.
;
; Prev. Hist. :	11-Mar-1994 (SLF) Wanted faster file_exist function
;
; History     :	Version 1, 11-Mar-1994, S. Freeland
;               Version 1.1 05-Jun-1998, J. Newmark changed loop to long
;               Version 1.2 15-Feb-2000, S.L.Freeland - work around RSI
;                 Version 5.3 non-backwardly compatible change....
;               Version 1.3 10-Mar-2000, S.L.Freeland - fix 5.3 /SIZE bug
;               Version 1.4 14-Jun-2002, S.L.Freeland - add FSTAT_ALL
;               Version 2, 1-Mar-2005, W.T. Thompson - Add support for files
;                                                      greater than 2 GB
;  
; Contact     :	SFREELAND
; Restrictions:
;   file size returned for directories under +5.3 is not valid
;-
;
if keyword_set(fstat_all) then template=fstat(-1) else begin
    if !version.os ge '5.2' then sizeval = long64(-1) else sizeval = -1L
    template = {name: '', size: sizeval}
endelse
; initialize some parameters
template.name=''
template.size=-1

if data_chk(files,/string) then begin
   template=replicate(template,n_elements(files))
   for i=0l,n_elements(files)-1 do begin
on_ioerror,err
      openr,lun,/get_lun,files(i)
      fstat_lun = fstat(lun)
      if keyword_set(fstat_all) then begin
          if n_elements(files) eq 1 then template = fstat_lun else $
            template(i) = fstat_lun
      end else begin
          template(i).name = fstat_lun.name
          template(i).size = fstat_lun.size
      endelse
      free_lun,lun
      goto,ok
err:       
   on_ioerror,null
ok:
   endfor
endif else begin
   message,/info,'Expect file or file array...'
   return,-1
endelse

case 1 of 
   keyword_set(exist): retval=template.name ne ''
   keyword_set(size):  retval=template.size
   keyword_set(fstat_all): retval=temporary(template)
   else: retval=template.name ne ''		; default to /exist
endcase

if since_version('5.3') and keyword_set(exist) then begin
  ssno=where(1-retval,ssncnt)
  if ssncnt gt 0 then retval(ssno)=is_dir(files(ssno))
endif  

return,retval

end
