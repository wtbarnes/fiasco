function is_number, inarray
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
;	Use READS and trap on any errors with ON_IOERROR
;HISTORY:
;	Written 29-Oct-97 by M.Morrison
;       Modified 29-Jun-99, Zarro (SM&A/GSFC) - added check for undefined input
;       Modified 29-Sep-00, Zarro (EIT/GSFC) - added check for invalid inputs
;       Modified 21-jun-05, Csillaghy (UAS Switzerland)- added call to is_number2
;       Modified 18-feb-15, Etesi (UAS Switzerland) - resetting error caused by reads
;-
;


if n_elements(inarray) eq 0 then return, 0b

if since_version( '5.4' ) then return, is_number2(inarray)

if (datatype(inarray,2) gt 7) then return, 0b

n = n_elements(inarray)
out = bytarr(n) + 1b
;
for i=0,n-1 do begin
    on_ioerror, err
    f = 0.
    reads, inarray[i], f
  goto, no_err
  err:
    catch, /cancel
    message, /reset
    out[i] = 0b
  no_err:
end
;
if (n_elements(out) eq 1) then out = out[0]
return, out
end
