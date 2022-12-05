;+
; Project     : HESSI
;
; Name        : ERR_STATE
;
; Purpose     : return !error_state.msg or !err_string
;
; Category    : utility help
;
; Syntax      : IDL> print,err_state()
;
; Inputs      : None
;
; Outputs     : !error_state.msg (if supported), else !err_string
;               sys_msg = system error message
;
; Keywords    : None
;
; History     : 6-Jan-2003, D. Zarro (EER/GSFC) - written
;               27 May-2017, Zarro (ADNET) - added sys_msg
;
; Contact     : dzarro@solar.stanford.edu
;-

function err_state,sys_msg

err=''
sys_msg=''

defsysv,'!error_state',exists=exists
if exists then begin
 s=execute('err=!error_state.msg')
 s=execute('sys_msg=!error_state.sys_msg')
 return,err
endif

defsysv,'!err_string',exists=exists
if exists then s=execute('err=!err_string')

return,err

end

