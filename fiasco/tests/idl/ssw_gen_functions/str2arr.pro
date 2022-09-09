function str2arr, instring,  delim,  array=array, delimit=delimit, $
      nomult=nomult, list=list
;+
;
; PROJECT: SSW
;
; NAME:
; 	STR2ARR	
;
; PURPOSE:
;	Convert delimited string into string array
;
; INPUT PARAMETERS:
;	instring  - delimited string to be split into components
;       delim     - delimiter to use (default=comma)
;
; OUTPUT:
;      function output is string array, n_elements=number delimiters+1
; 
; KEYWORD PARAMETERS:
;      delimit - delimiter - equiv to positional DELIM; for backward compat.
;      numult  - switch, if set, dont return nulls for consecutive delimiters
;      list    - if set, display array to terminal (via more)
;      array  (output) - string array, number elements=number of delimeters+1
;                        (same as function output)
;
; CALLING SEQUENCE:
;	array=STR2ARR(string [, delimiter, delimit=delimit, /nomult, /list] 
;
; CALLING EXAMPLES:
;       IDL> more,str2arr('this,is,a,test')                    ; default delim
;       this                                                   ; (ie: comma)
;       is
;       a
;       test
;
;       IDL> more,str2arr('this$$$is$$$another$$$test','$$$') ; delim='$$$'
;       this
;       is
;       a
;       test
;
;
; COMMON BLOCKS;
;	NONE
;
; MODIFICATION HISTORY:
;	Version 0 - Sam Freeland (Yohkoh)
;	slf - feature correction (occur is array, then prob.)
;	slf - 25-feb-92 - added positional delim parameter
;       slf -  2-feb-93 - changed recursion to loop for memory problems
;       slf - 19-mar-93 - optimize case where delimiter is 1 character (comma)
;       slf - 20-mar-93 - fixed a minor bug with major implications
;       slf - 16-Jan-97 - merge Dave Pike (RAL) changes for SSW GEN 
;                         (LIST and NOMULT keyword and function)
;       slf -  8-may-97 - added to header
;       Zarro (SM&A/GSFC), 14-Sep-99, threw in couple of temporary's for memory
;-
;
list=keyword_set(list)
nomult=keyword_set(nomult)

if n_params() eq 2 then delimit=delim		; slf, 25-feb-92
if not keyword_set(delimit) then delimit=','
delim_len=strlen(delimit)
if n_elements(array) eq 0 then array=''
text=instring(0)					 ; expect/make scaler
maxlen=strlen(text)

; slf, optimize case where delimiter is 1 character long      - BYTE operations
if delim_len eq 1 then begin
   bdelim=byte(delimit)
   btext=byte(text)
   wdelim=where(btext eq bdelim(0),dcount)
   if dcount gt 0 then begin
      wdelim=[0,wdelim+1]      
      sizewd=deriv_arr(wdelim)-1
      array=strarr(dcount+1)
      for i=0,dcount-1 do begin
         array(i)=strmid(text,wdelim(i),sizewd(i))
      endfor
      array(i)=strmid(text,wdelim(i),strlen(text)-wdelim(i) )      
   endif else array=text
endif else begin
   occur=strpos(text,delimit)
;
   while occur(0) ne -1 do begin
      substring=strmid(text,0,occur)
      array=[temporary(array),substring]
      text=strmid(text,occur+delim_len,maxlen)
      occur=strpos(text,delimit)
   endwhile
   array=[temporary(array),text]
   array=array(1:*)
endelse


if nomult then begin
   nnulls=where(array ne '',nncnt)
   if nncnt gt 0 then array=array(nnulls)   
endif
if list then more,array

;
return,array
end 
