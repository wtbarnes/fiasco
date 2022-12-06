function get_logenv, logenvs, count=count, status=status, $
        environ=environ, home=home, curdir=curdir, outenv=outenv, $
	full_trans=full_trans, uniq=uniq, case_ignore=case_ignore
;+
;   Name: get_logenv
;
;   Purpose: provide system independent front end to get_env and trn_log
;
;   Input Parameters:
;      logenvs - string or string array of vms logicals or unix environments
;		 (may contain wild cards)
;      
;   Output Parameters:
;	function returns equivilence string/string array
;
;   Optional Keyword Parameters:
;      outenv - actual environmentals (output) - may differ from logenvs
;		if wild cards were used (null if associated env not defined)
;
;   Calling Sequence:
;      envs=get_logenv(logenvs [,count=count, outenv=outenv, 
;
;   Calling Examples:
;      sfddirs=get_logenv('DIR_SFD*')		  ; sfd directories
;      gendirs=get_logenv('dir_gen*',/case_ignor) ; all DIR_GENxxxxx
;      home   =get_logenv(/home)		  ; home directory
;      curdir =get_logenv(/curdirq)		  ; current [same as curdir()]
;      ysenvs =get_logenv('ys_*',/environ)	  ; environmentals form ys_...
;      ysflags=get_logenv('ys_*')		  ; same, only translated value
;
;   Optional Keyword Parameters
;	count  (output) - number of non-null elements returned
; 	status (output) - boolean success vector (1:match, 0:nomatch)
;       outenv (output) - string array - actual environmentals found
;	environ(input)  - switch - return logical/environ, not translation
;			  (could be useful if logenvs contains wild card)
;	home   (input)  - switch, if set, return home directory translation
;       curdir (input)  - if set, return current directory
;	case_ignore (input) if set, case INsensitive (only functional in unix)
;     
;   History:
;      slf - 5-apr-1993
;      SLF - 13-Jan-1993 (use /noshell on printenv for speed)
;      SLF - 16-Jan-1993 (work around for multiple line envs)
;      SLF - 30-Mar-1994 fix VMS bug
;      DMZ - 20-May-1994 added check for '$' in first character
;      SLF - 21-May-1994 added a semicolon to comment (4 man-months effort)
;      SLF -  2-Jun-1994 fix (define allenv if wild cards used)
;      JSN - 03-Jun-1998 changed loops to long integer
;
;   Restrictions:
;      if logenvs is an array with one or more elements containg a wild
;      card AND more than one env/log satisfies the wildcard, then only
;      one match is returned (maintains 1 to 1 correspond between 
;      input <logenvs> and <output>  
;      Does not yet handle multiple line environmentals (1st line returned)
;-
env=keyword_set(environ)

sav_logenvs=logenvs  ;-- save inputs to avoid clobbering

;-- remove '$' in first character so that UNIX getenv/printenv 
;    will work correctly

if !version.os ne 'vms' then begin  
 for i=0l,n_elements(logenvs)-1 do begin
  ilog=logenvs(i) & doll=strpos(ilog,'$')
  if doll eq 0 then logenvs(i)=strmid(ilog,1,strlen(ilog))
 endfor
endif
 

status=1 & count=0 			; success vector (found/not found)
outenv=''
case_ignore=keyword_set(case_ignore)	; if set, case INsensitive

case 1 of 
;  process keywords first (special cases)
   keyword_set(home):   retval=getenv('HOME')	; home directory
   keyword_set(curdir): cd,current=retval	; same as curdir()
;
   else: begin
; initialize
      status=0 & count = 0
      retval=strarr(n_elements(logenvs))	; return value
      outenv=strarr(n_elements(logenvs))	; ouput environmentals
      status=intarr(n_elements(logenvs))	; boolean status
      wc_check=where(strpos(logenvs,'*') ne -1,wc_count) ; wild cards??
      wcs=wc_count gt 0
      mwcs=wcs and n_elements(logenvs) gt 1	; multiple wild card entries
      case strlowcase(!version.os) of
         'vms': begin
	    if wcs then begin
;	       if not mcs then spawn,'show log ' + logenvs(0),wcout
               message,/info,'vms wild cards not implemented'
	    endif else begin
	       for i=0l, n_elements(logenvs)-1 do begin
		  temp=''
                  stat=trnlog(logenvs(i), temp)
                  if stat then begin
                     status(i)=stat
		     retval(i)=temp
                  endif
                endfor
             endelse
         endcase
;
         else: begin
            if wcs then begin
	       spawn,'printenv',allenv,/noshell
	       lead_dollar=where(strpos(allenv,'$') eq 1,dollcnt)
               if dollcnt gt 0 then allenv(lead_dollar)= $
		  strmid(allenv(lead_dolloar),1,max(strlen(allenv)))
	       equals=strpos(allenv,'=')
	       equalw=where(equals ne -1, ecount)
	       allenv=allenv(equalw)
	       equals=equals(equalw)
;	       gotta be a better way (strmid for vectors)
	       environ='' & trans=''
	       for i=0l,n_elements(equalw)-1 do begin
	          environ=[environ,strmid(allenv((i)), 0, equals(i))]
		  trans  =[trans, strmid(allenv((i)), equals(i)+1, $
		     strlen(allenv((i))))]
	       endfor
	       environ=environ(1:*) & trans=trans(1:*)
               for i=0l, n_elements(logenvs)-1 do begin
	          match=wc_where(environ, logenvs(i), mcount, $
	             case_ignore=case_ignore)
		  if mcount gt 0 then begin
                     if mwcs then begin
		        retval(i)=trans(match(0)) 
			outenv(i)=environ(match(0))
		     endif else begin
		        retval=trans(match)		
			outenv=environ(match)
		     endelse
		  endif
	       endfor
               miss=where( (strmid(outenv,0,1) ne '') and $
                           (strmid(outenv,0,1) ne '$'),mcnt)
               if mcnt gt 0 then outenv(miss)='$' + outenv(miss)

	       if env then retval=outenv

	    endif else begin
	       for i=0l,n_elements(logenvs)-1 do begin
                  retval(i)=getenv(logenvs(i))
		  if case_ignore and retval(i) eq '' then begin
                    retval(i)=getenv(strupcase(logenvs(i)))	; try upper
		    if retval(i) eq '' then $
 		       retval(i)=getenv(strlowcase(logenvs(i)))  ; try lower
		  endif
                  status(i)=retval(i) ne ''
               endfor
               ss=where(status,cnt)
               if cnt gt 0 then outenv(ss)=logenvs(ss)
               miss=where( (strmid(outenv,0,1) ne '') and $
                           (strmid(outenv,0,1) ne '$'),mcnt)
               if mcnt gt 0 then outenv(miss)='$' + outenv(miss)
               if keyword_set(environ) then retval=outenv
	    endelse
	 endcase
      endcase
   endcase
endcase        

if keyword_set(uniq) and n_elements(retval) gt 1 then begin
   uretval=uniq(retval,sort(retval))
   retval=retval(uretval)
   outenv=outenv(uretval)
endif

if n_elements(retval) eq 1 then retval=retval(0)	; make scaler

if n_elements(outenv) gt 1 then count=n_elements(outenv)

status=outenv ne ''

logenvs=sav_logenvs

return,retval

end
