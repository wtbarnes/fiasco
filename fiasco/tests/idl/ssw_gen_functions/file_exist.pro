;+
;   Name: file_exist
;
;   Purpose: returns true(1) if any files match pattern=file
;	     false(0) if no files found
;
;   Input Parameters:
;      file - file, pathname or file filter to check
;
;   Optional Keyword Parameters
;      direct - set if want full check (slower) - use if might be an
;		empty directory you are looking for.
;      old - use pre-2022 version that uses findfile, otherwise new version
;               uses file_test and file_search (new version ignores direct keyword)
;
;   History: written by slf, 11/30/91
;	4-Sep-92 (MDM) - Modified to handle an array of input files
;      23-sep-93 (slf) - make sure null file returns false
;			 (findfile count=1! when null requested)
;      17-Mar-94 (SLF) - use file_stat(/exist) for non wild card cases
;      April Fools'94 (DMZ) -- added check for :: in file name
;                              (file_stat doesn't work across DECNET)
;      31-May-94 (DMZ)      -- added check for VMS directory name input
;                              (file_stat doesn't seem to work on VMS
;                               dir names like '[ys.atest.soft]')
;      16-Aug-94 (DMZ)      -- added another VMS patch to handle
;                              case when user enters  a VMS subdirectory name
;                              without a file name.
;      6-Aug-97 (Zarro/GSFC) -- added check for nonstring input
;      30-may-99 (rdb)      -- stop use of file_stat with WINDOWS
;                              cause access violation - for wildcard search
;      8-Jan-1999 - S.L.Freeland - put in IDL version check since 5.3/IRIX
;                   (at least) not backwardly compatible (must use findfile)
;     10-Jan-1999 - S.L.Freeland - extended 8-jan mod to all UNIX 5.3
;     15-Feb-2000 - S.L.Freeland - removed 5.3 changes (moved fix to file_stat)
;     10-Feb-2005 - Kim Tolbert - changed '/' to get_delim()
;     24-oct-2021 - Zarro (ADNET) - updated to use more robust FILE_SEARCH
;                   (old version uses findfile, which fails in certain cases in IDL 8.8.1)
;     27-Feb-2022 - Zarro (ADNET) - updated to expand environment variable using more robust FILE_TEST2
;                   (intrinsic FILE_TEST doesn't expand on Windows)
;     
;-
function file_exist, file, direct=direct, verbose=verbose,old=old

  new=~keyword_set(old)
  direct=keyword_set(direct)
  if new then begin
   if n_elements(file) eq 0 then return,0b
   nf=n_elements(file)
   out=bytarr(nf)
   if is_blank(file) then if nf eq 1 then return,0b else return,out
 
;-- search for wildcard filenames

   chk=stregex(file,'\*',/bool)
   wc=where(chk,wcount,complement=nwc,ncomplement=nwcount)
   if (wcount gt 0) then begin
    for i=0,wcount-1 do begin
     test=file_search(file[wc[i]],count=count,/expand_env)
     if count gt 0 then out[wc[i]]= 1b
    endfor
   endif

;-- search for remaining filenames

   if (nwcount gt 0) then out[nwc]=file_test2(file[nwc],/expand_env)

   if nf eq 1 then out=out[0]
   return,out
  endif
  
;----------------------------------------------------------------------------------------------------------

  if datatype(file) ne 'STR' then return,0
  verbose=keyword_set(verbose)
  direct=keyword_set(direct)

  file_sav=file            ;-- protect input

  n = n_elements(file)
  out = bytarr(n)

  ;-- check for VMS directory name only

  if !version.os eq 'vms' then begin
    break_file,file,dsk,dir,name,ext  ;-- is this a directory name only?
    dcp=where( (strtrim(name,2) eq '') and (strtrim(dir,2) ne ''),dcount)
    if (dcount gt 0) then file(dcp)=file(dcp)+'*.*'  ;-- need to add wild chars
  endif

  wc=  where( (strpos(file,'*') ne -1) or (strpos(file,'::') ne -1),wccnt)
  nowc=where((strpos(file,'*') eq -1) and (strpos(file,'::') eq -1),nowccnt)

  case 1 of
    strlowcase(!version.os_family) eq 'windows': begin
      if verbose then box_message,'Using findfile method'
      wc=indgen(n_elements(file))
      wccnt = n_elements(wc)
    endcase
    else: if nowccnt gt 0 then out(nowc)=file_stat(file(nowc),/exist)
  endcase

  ; now do wild card (file filter) logic
  for i=0, wccnt-1 do begin
    case keyword_set(direct) of
      0:test=findfile(file(wc(i)),count=count)
      ;
      ;  if direct keyword set, can detect empty directories
      ;
      1:begin
        filedelim=str_lastpos(file(i),get_delim())
        name=strmid(file(wc(i)),filedelim+1,strlen(file(wc(i)))-filedelim)
        path=strmid(file(wc(i)),0,filedelim)
;        help,name,path
        files=findfile(path)
        found=where(name eq files)
        count = found(0) ge 0
      end
    endcase
    out(wc(i)) = count ne 0
  end
  ;
  out = out and (file ne '')
  if (n eq 1) then out = out(0)	;turn it into a scalar

  file=file_sav
  return, out
end

