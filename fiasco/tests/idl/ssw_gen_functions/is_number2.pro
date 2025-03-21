function is_number2, inarray
;+
;NAME:
;	is_number
;PURPOSE:
;	Test string (or array of strings) to see if it is a number
;SAMPLE CALLING SEQUENCE:
;	out = is_number(inarray)
;	print, is_number(v)
;	print, is_number('xx')
;	out = is_number('66.6e')
;INPUT:
;	inarray - The string(s) to test
;OUTPUT:
;	out	- Boolean array (0 means not a string, 1 means it is)
;METHOD:
;	Use READS and trap on any errors with (ON_IOERROR earlier) now with catch       
;HISTORY:
;	Written 29-Oct-97 by M.Morrison
;       Modified 29-Jun-99, Zarro (SM&A/GSFC) - added check for undefined input
;       Modified 29-Sep-00, Zarro (EIT/GSFC) - added check for invalid inputs
;       Modified 21-jun-05, Csillaghy (UAS Switzerland) -- make sure it does
;                           not crash when 'information' or 'nan...' are
;                           passed. Also, move from on_ioerror to catch.
;       Modified 18-feb-15, Etesi (UAS Switzerland) - resetting error caused by reads
;-
;

; acs factor out n_elements
n = n_elements( inarray )
if n eq 0 then return, 0b

; acs change datatype to size /type
;or (datatype(inarray,2) gt 7) then return,0b
arrtype = size( inarray, /type )
out =  n gt 1 ?  bytarr(n) : 0B

if arrtype lt 7 or arrtype eq 9 or arrtype gt 11 then return, out + 1B 
if arrtype ne 7 then return, out
 
; acs see above
;n = n_elements(inarray)

for i=0,n-1 do begin

    err_no = 0
    catch, err_no
    if err_no ne 0 then begin 
; if the string is not a valid number here it will generate an error and go to
; the next iteration
        catch, /cancel
        message, /reset
        continue
    endif

    f = 0.
    reads, inarray[i], f

    ; now correct in case inarray[i] starts with 'inf' or 'nan'
    if not finite( inarray[i] ) and strlen( inarray[i] ) gt 3 then begin
        out[i] = 0b
    endif else begin 
; here no problems 
        out[i] = 1b
    endelse

end
;
if (n_elements(out) eq 1) then out = out[0]

return, out

end
