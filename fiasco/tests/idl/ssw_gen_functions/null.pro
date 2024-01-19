;+
; Project     : VSO
;
; Name        : NULL
;
; Purpose     : Return !NULL 
;
; Category    : Utility
;
; Inputs      : None
;
; Outputs     : !NULL = '' if !NULL not defined
;
; Keywords    : None
;
; History     : 9-Dec-2015, Zarro (ADNET) - written
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-

function null

null=''
defsysv,'!null',exists=i
if i eq 0 then defsysv,'!null',''
return,!null
end
