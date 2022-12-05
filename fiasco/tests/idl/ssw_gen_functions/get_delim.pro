;+
; NAME:
;	get_delim
; PURPOSE:
;	returns file delimiter that is appropriate to VMS or UNIX
; CALLING SEQUENCE:
;	delim=get_delim()
; INPUTS:
;	none
; OUTPUTS:
;       delim=':' if VMS, '/' otherwise
; PROCEDURE:
;	checks !version.os system variable
; MODIFICATION HISTORY:
;       Written DMZ (ARC) May 1992
;       Modified DMZ(SAC) Sept 1997 - added Windows

function get_delim

on_error,1

os=strlowcase(os_family())

case os of
 'vms'    : delim=':'
 'windows': delim='\'
 else     : delim='/'
endcase

return,delim & end
