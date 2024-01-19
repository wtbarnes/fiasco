;+
; Project     : SOHO - CDS
;
; Name        : STRPAD
;
; Purpose     : Pads a string with blanks (or whatever) to specified width
;
;
; Explanation : Spaces are added in front of the input string to
;		make the returned string have LENGTH characters.
;		Use /AFTER to add spaces at the end.
;
; Use         : str = strpad(value, length, /after, fill=fill])
;
; Inputs      : VALUE:  A string..
;		LENGTH: The desired length of the result in characters
;
; Opt. Inputs : None.
;
; Outputs     : Returns the padded string.
;
; Opt. Outputs: None.
;
; Keywords    : AFTER : Set to add spaces at the end.
;               FILL  : The character with which to pad out the string.
;                        Default is the space character
;
; Restrictions: Value must be a scalar string
;
; Side effects: If the input string is longer than the desired
;		width, it is returned without truncation
;
; Category    : Utilities, Strings
;
; Written     : Stein Vidar Hagfors Haugan, 27 September 1993
;
; Modified    : Corrected typo which stopped /after working.  CDP, 28-Sep-94
;               Increased possible length used.  CDP, 22-Dec-94
;               Handle arrays of strings.  CDP, 16-Mar-95
;               Add FILL keyword and stop overwriting input. CDP, 30-Apr-95
;		Vectorized, richard.schwartz@gsfc.nasa.gov, 23-jan-2003
;               Vectorized even better, Zarro (EER/GSFC), 24-Jan-2003
;               Fixed degenerate dimension bug, Zarro (EER/GSFC), 29-Mar-2003
;-


function strpad, in_txt, length, after=after, fill=fill, no_copy=no_copy

;
;  check basic input
;

sz=size(in_txt)
if sz[n_elements(sz)-2] ne 7 then return,''
if exist(fill) then fill=fill else fill = ' '
if not exist(length) then return,in_txt

;-- only process strings with length less than new length

tlen=strlen(in_txt)
mlen=max(tlen)
process=where(tlen lt length,lcount)
if lcount eq 0 then return,in_txt else begin
 chk=where(tlen ge length,pcount)
 if pcount gt 0 then keep=in_txt[chk]
 mlen=max(tlen[process])
endelse
mlen=mlen > 1

;-- convert to byte array

if keyword_set(no_copy) then byte_in=temporary(in_txt) else byte_in=in_txt

;-- remove degenerate dimensions

if sz[0] gt 0 then byte_in=reform(byte_in)

byte_in = transpose(byte(strmid(temporary(byte_in),0,mlen)))

;-- pad with 32b blanks

blank=where(byte_in eq 0b,count)
if count gt 0 then byte_in[blank]=(byte(' '))[0]

;-- create output byte array 

byte_out = bytarr(n_elements(byte_in[*,0]),length) + (byte(fill))[0]

if keyword_set(after) then byte_out[0,0]=temporary(byte_in) else $
 byte_out[0,length-mlen]=temporary(byte_in)

;-- convert back to string

byte_out=string(transpose(temporary(byte_out)))

;-- insert back unprocessed strings

if pcount gt 0 then byte_out[chk]=temporary(keep)

return,byte_out

end
