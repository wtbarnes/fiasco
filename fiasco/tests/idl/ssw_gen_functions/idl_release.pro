;+
; Project     : SOHO - CDS     
;                   
; Name        : IDL_RELEASE
;               
; Purpose     : check if IDL release version within specified range
;               
; Category    : system
;               
; Explanation : 
;               
; Syntax      : IDL> a=idl_release(lower=lower,upper=upper)
;    
; Examples    :
;
; Inputs      : None
;               
; Opt. Inputs : 
;               
; Outputs     : 1/0 if IDL version is within specified range
;
; Opt. Outputs: None
;               
; Keywords    : LOWER = lower version to check
;               UPPER = upper version to check
;               INCLUSIVE = make check inclusive
;               VERS = IDL version
;
; Common      : None
;               
; Restrictions: None
;               
; Side effects: None.
;               
; History     : Version 1,  27-Feb-1997,  D M Zarro.  Written
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-            


function idl_release,upper=upper,lower=lower,inclusive=inclusive,vers=vers

if not exist(upper) then upper=100
if not exist(lower) then lower=0

vers=float(strmid(!version.release,0,3))

if keyword_set(inclusive) then begin
 ok=(vers ge lower) and (vers le upper)
endif else begin
 ok=(vers gt lower) and (vers lt upper)
endelse

return,ok
end


