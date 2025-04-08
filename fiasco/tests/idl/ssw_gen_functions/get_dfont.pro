;+
; Project     : SOHO - CDS
;
; Name        :
;	GET_DFONT()
; Purpose     :
;       Return widget font with size compatible with current device
; Explanation :
;	Useful for selecting fonts to fit into  widgets
; Use         :
;	Result = GET_DFONT(UFONT)
; Inputs      :
;       UFONT = user's optional input fonts (string array)
; Opt. Inputs :
;	None.
; Outputs     :
;            DFONT = returned fonts (string array)
; Prev. Hist. :
;       Written Elaine Einfalt (HSTX) May 1992.
; Modified    :
;       Version 1, Dominic Zarro, GSFC, 1 August 1994.
;               Corrected some bugs
;	Version 2, 23-Oct-1997, William Thompson, GSFC,
;		Only select X device if on Unix or VMS.
;		Call SELECT_WINDOWS
;       Modified, 18-Jun-01, Zarro (EITI/GSFC) - added call to IS_WOPEN
;       Modified, 10-Jul-02, Kim Tolbert - added branches for Windows
;       3-Feb-2017, Zarro (ADNET)
;        - deprecated as some device fonts were incompatible with some
;          window managers.
;        27-Feb-2017, Zarro (ADNET)
;        - only deprecate for Mac OS
;       9-March-2018, Zarro (ADNET)
;       - added CATCH for window managers that don't support certain fonts
;
;-

function get_dfont,ufont

if stregex(!version.os_name,'^mac',/bool,/fold) then begin
 if is_string(ufont) then return,ufont[0] else return,''
endif

;
; If no window currently exists then create a window,
; or else the DEVICE,FONT command will create one, but won't delete it later.
;

if ~allow_windows(/quiet) then return,''

if ~is_wopen(!d.window) then begin
 window,/pix,/free,xsize=1,ysize=1
 wpix=!d.window
endif

if datatype(ufont) ne 'STR' then begin
   case os_family() of
     'unix': ufont = ["*menu*-r-*--12-*", 		$
             "*menu*-r-*--13-*", 			$
	         "*helvetica*bold*-r-*--12-*", 	$
	         "*time*-r-*--12-*"]
     'Windows': ufont=["MS Sans Serif", "Arial", "Times New Roman"]
      else: ufont="helvetica"
   endcase
endif

floop 	 = n_elements(ufont)			; number of different searches
got_one = 0
counter = 0

;
; Search the hardware fonts for a font that matches the user's input.
;

error=0
catch,error
if error ne 0 then begin
 err=err_state()
 mprint,err
 message,/reset
 catch,/cancel
 if is_string(ufont) then return,ufont[0] else return,''
endif

while (counter lt floop) && (got_one eq 0) do begin
   ;
   ; 1) see if any system hardware fonts fit the font designation.
   ; 2) return those matches in string or string array DEFNAMES,
   ; 3) and, the number of matches are returned in in DEF_NUM.
   ;

 device, font=ufont[counter], get_fontname=defnames, get_fontnum=def_num

 if def_num gt 0 then got_one = 1 	; there is at least one font match
 counter = counter + 1
endwhile

if def_num ne 0 then begin
   ; on windows, loses any attributes like size, bold, etc, so return original name passed in
   if os_family() eq 'Windows' then dfont=ufont[counter-1] else dfont=defnames
endif else dfont=''

if n_elements(dfont) eq 1 then dfont=dfont[0]

wdel,wpix
return,dfont
end


