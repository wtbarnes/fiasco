;+
; Project     : VSO
;
; Name        : IS_WINDOWS
;
; Purpose     : Return true if Windows OS
;
; Category    : utility system
;
; Syntax      : IDL> chk=is_windows()
;
; Inputs      : None
;
; Outputs     : CHK = true if Windows OS
;
; Keywords    : None
;
; History     : 23-Feb-2022, Zarro (ADNET)
;-

function is_windows

return,os_family(/lower) eq 'windows'

end

