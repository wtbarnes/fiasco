function DERIV_arr, INARR
;
;+
;NAME:
;	deriv_arr
;PURPOSE:
;	Procedure to return the rate of change of the
;	input array
;INPUT:
;	inarr	- the vector to find the derivative of
;OUTPUT:
;	returns	- If the input is "n" elements, the
;		  output is "n-1" elements.
;HISTORY:
;	Written 1988 by M.Morrison
;	17-Apr-92 (MDM) - Changed to be a "lindgen" 
;	 8-Sep-93 (MDM) - Return value as scalar if only one value
;-
;
N=N_ELEMENTS(INARR)
ss = lindgen(n-1)
outarr = inarr(ss+1) - inarr(ss)
;
if (n_elements(outarr) eq 1) then outarr = outarr(0)
RETURN, outarr
END
