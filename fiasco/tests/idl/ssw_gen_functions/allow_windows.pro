;+
; Project     : HESSI
;                  
; Name        : ALLOW_WINDOWS
;               
; Purpose     : platform/OS independent check if current device
;               allows windows 
;                             
; Category    : system utility
;               
; Syntax      : IDL> a=allow_windows()
;                                        
; Outputs     : 1/0 if yes/no
;
; Keyword     : REMOTE - check if remote window can be opened
;                   
; History     : Version 1,  4-Nov-1999, Zarro (SM&A/GSFC)
;               13 Dec 2001, Zarro (EITI/GSFC) - added DEVICE call
;               22 Oct 2002, Zarro (EER/GSFC) - added FSTAT check
;               25 Feb 2012, Zarro (ADNET) - added /REMOTE
;               10 Mar 2022, Zarro (ADNET) - cleaned up
;
; Contact     : dzarro@solar.stanford.edu
;-    

function allow_windows,verbose=verbose,err=err,_extra=extra,remote=remote

;-- catch any open errors (vers > 3)

err=''
werr='Error opening widget window.'
error=0
catch,error
if error ne 0 then begin
 err=werr
 if keyword_set(verbose) then mprint,err,/info
 catch,/cancel
 message,/reset
 return,0b
endif

wdisp=widget_info(/display)
if ~wdisp then begin
 err=werr & return,0b
endif   

;-- check if windows or widgets are even supported

if ~have_windows() || ~have_widgets() then begin
 err=werr
 return,0b
endif

if keyword_set(remote) then begin
 check=remote_window(err=err,verbose=verbose)
 return,check
endif

device,get_screen_size=s

return,1b

end
