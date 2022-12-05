function wc_where, inarray, inpattern, mcount, case_ignore=case_ignore
;
;+
;   Name: wc_where
; 
;   Purpose: return subscripts of input array where a pattern match is
;	     found - allows use of multiple wild card characters (*)
;
;   Input Paramters:
;      inarray - string array to search
;      inpattern - string (scaler) to match - may include wild cards (*)
;
;   Output:
;      function returns subscripts of inarray where match found (-1 if none)
;      mcount - number of matches found 
;
;   Calling Examples:
;      ss=wc_where(files_arr,'*9201*',mcount)	
;      ss=wc_where(files_arr,'sfr*1230*')	
;      ss=wc_where(routine_arr,'*time*.pro',mcount)	
;
;   History:
;      slf, 8-Jan-1993
;      slf,13-Jan-1993 - patched last segment logic
;      slf,15-Jan-1993 - added case_ignore keyword
;      slf,12-apr-1993 - dont clobber pattern via case_ignore
;      slf,17-feb-1994 - fix bug in last segment logic
;      slf,30-jun-1994 - call wc_whereq if pattern has embedded "?" character
;      acs,23-jul-2004 - checks for idl version tu use strmatch if avail
;-
;
;

if since_version( 5.3 ) then return, where( strmatch( inarray, inpattern, fold_case = case_ignore ), mcount )
 
if (strpos(inpattern,'?'))(0) ne -1 then $
   ss=wc_whereq(inarray,inpattern,mcount,case_ignore=case_ignore) $
else begin
;  initialize search array and output vector
   newarray = inarray
   pattern  = inpattern
   ss=lindgen(n_elements(newarray))	
   remss=ss
; 
; slf, 15-jan - ignore case on request
   if keyword_set(case_ignore) then begin
      newarray=strlowcase(newarray)
      pattern=strlowcase(pattern)
  endif
;
;  break pattern into search segments
   pparts=str2arr(pattern,'*')		; break at wild card (*)
   nparts=n_elements(pparts)
   case nparts of
      1: pparts=[pparts,pparts]		; single segment
      2:					; 1 wildcard
      else: midparts=pparts(1:nparts-2)	; central segments
   endcase         
;
   firstpart=pparts(0)			; null if leading wild card
   lastpart=pparts(n_elements(pparts)-1)	; null if trailing wild card
   nmparts=n_elements(midparts)		; central segments
; -----------------------------------------------------------------------
; ------------------------------------------------------------------------
;  check for leading wild card - if none, process first segment
   remcount=1
   if firstpart ne '' then begin
;     no leading wc, match pattern to 1st character position
      ss=where(strpos(newarray,firstpart) eq 0,mcount)
      if mcount gt 0 then newarray(ss)=strrempat(newarray(ss),firstpart, $
	 patss=remss, remcount, /trunc)
   endif
; ------------------------------------------------------------------------

; ------------------------------------------------------------------------
;  do middle segments 
   midcount=0
   while ss(0) ne -1 and midcount lt (nmparts) do begin
      remarr=strrempat(newarray(ss),midparts(midcount),/trunc,remcount, $
	patss=remss)
      if remcount gt 0 then begin
         newarray(ss(remss))=remarr(remss)
         ss=ss(remss)
      endif else ss=remss
      midcount=midcount+1
   endwhile
;
; ------------------------------------------------------------------------
;    do last segment
;    (got to be a better way...)

   if lastpart ne '' and ss(0) ne -1 then begin
      vcount=0
      if firstpart ne lastpart then begin		; should already be null 
         valid=where(newarray(ss) ne '', mcount)   
         if mcount gt 0 then begin
            ss=ss(valid)
            pos=str_lastpos(newarray(ss),lastpart)
            valid=where((pos ne -1) and $
   	      ( (pos + strlen (lastpart)) eq strlen(newarray(ss))),vcount)
         endif 
      
         if vcount gt 0 then ss=ss(valid) else ss=-1
      endif else ss = where(newarray eq '')
   endif
; --------------------------------------------------------------------------
;
   if ss(0) ne -1 then mcount=n_elements(ss) else mcount=0
endelse

;
return,ss
end



