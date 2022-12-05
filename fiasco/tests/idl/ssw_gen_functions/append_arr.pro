;+
; Project     : SOHO - CDS
;
; Name        : APPEND_ARR
;
; Category    : Utility
;
; Purpose     : Manageable method for concatanating arrays
;
; Syntax      : IDL> result=append_arr(input,append)
;
; Inputs      : INPUT = array (or scalar) that requires appending
;               APPEND = array (or scalar) to append
;
; Outputs     : concatanated arrays
;
; Keywords    : NO_COPY = set to not create internal copy of
;               INPUT & APPEND (it will be destroyed)
;
; History     : Written: 1-Oct-1998, Zarro (SM&A/GSFC)
;               Modified: 26-March-2000, Zarro - sped up with SIZE and CATCH
;               Modified: 25-Feb-2003, Zarro - improved error checking
;               26-Dec-2018, Zarro (ADNET) - modernized
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-

function append_arr,input,append,no_copy=no_copy,_extra=extra

error=0
catch,error
if error ne 0 then begin
 mprint,err_state()
 catch,/cancel
 return,-1
endif

if exist(input) && exist(append) then begin
 if keyword_set(no_copy) then return,[temporary(input),append] else $
  return,[input,append]
endif

if exist(input) && ~exist(append) then return,input
if ~exist(input) && exist(append) then return,append

return,-1

end

