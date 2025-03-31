;+
; Project     : SSW     
;                   
; Name        : ARR2STR()
;               
; Purpose     : Convert an array to a delimited string.
;               
; Explanation : 
;               
; Use         : IDL> s = arr2str(array,'-')
;                    s = arr2str(array,delim='-')
;    
; Inputs      : arr  -  input data array which is to be converted
;                       to a simple string.
;               
; Opt. Inputs : delim - an alternative positional parameter to specify the
;                       delimiter.
;               
; Outputs     : Function returns a simple string using the delimiter to 
;               separate the elements of the original array.
;               
; Opt. Outputs: 
;               
; Keywords    : delimiter  -  specify the delimiter to be used, default 
;                             delimiter is ','
;		trim_str   -  If set, call STRTRIM when converting to string
;		compress   -  If set, call STRCOMPRESS after converting
;               no_duplicate  If set, inhibit having string with consecutive
;                             delimiters such as // 
;
; Category    : Util, string
;               
; Prev. Hist. : Sam Freeland 11/19/91 
;               (Various Slf,MDM,DP mods)
;
; Written     : Sam Freeland 
;               
; Modified    : Version 2, William Thompson, GSFC, 15 June 1995
;			Added /TRIM keyword to be compatible with Yohkoh
;			version.  Added /COMPRESS keyword
;               Version 2.1, Sam Freeland, SSW merge
;               Version 3, Zarro (SAC/GSFC) - added /NO_DUPLICATE &
;                       renamed TRIM keyword to TRIM_STR to avoid
;                       name conflict with TRIM function
;               Modified, 17-Mar-03, Zarro (EER/GSFC) - changed loop variable
;                       to long, change to use [], added error checks,
;                       changed 'string' variable name to 'ostring' in case
;                       of reserved name conflicts
;               Modified, 20-April-12, Zarro (ADNET) 
;                       Vectorized and sped up with strjoin.
;-            

function arr2str, starray, delim, delimiter=delimiter, trim_str=trim_str,$
	compress=compress,no_duplicate=no_duplicate,_extra=extra

if n_elements(starray) eq 0 then return,''

;
;  delimiter specified as positional parameter
;

if n_params() eq 2 then delimiter = delim
if (n_elements(delimiter) eq 0) then delimiter=','

if keyword_set(trim_str) then $
 ostring=strjoin(strtrim(starray,2),delimiter,/single,_extra=extra) else $
  ostring=strjoin(starray,delimiter,/single,_extra=extra)

if keyword_set(no_duplicate) then $
 ostring=str_replace(ostring,delimiter+delimiter,delimiter)

if keyword_set(compress) then ostring=strcompress(ostring,_extra=extra)

return,ostring

;------------------------------------------------------------------
;-- pre-STRJOIN

sz=size(starray)
dtype=sz[n_elements(sz)-2]
if dtype eq 7 then strings=starray else strings=string(starray)
if (keyword_set(trim_str)) then strings=(strtrim(strings,2))

;
;  concatenate elements with required delimiter
;

ostring=strings[0]
no_dup=keyword_set(no_duplicate)
np=n_elements(starray)
for i=1l,np-1 do begin
 temp_limiter=delimiter
 if no_dup then begin
  first_char=strmid(strtrim(strings[i],2),0,1)
  if first_char eq delimiter then temp_limiter=''
 endif
 ostring = ostring + temp_limiter + strings[i]
endfor
if keyword_set(compress) then ostring=strcompress(ostring)

return,ostring

end
