;+
; Project     :	SDAC
;
; Name        :	EXIST
;
; Purpose     :	To See if variable Exists
;
; Explanation :	So obvious, that explaining it will take more
;               lines than the code.
;
; Use         :	A=EXIST(VAR)
;
; Inputs      :	VAR = variable name
;
; Opt. Inputs : None.
;
; Outputs     :	1/0 if variable exists or not
;
; Opt. Outputs:	None.
;
; Keywords    :	None.
;
; Calls       :	None.
;
; Common      :	None.
;
; Restrictions:	None.
;
; Side effects:	None.
;
; Category    :	Useful stuff
;
; Prev. Hist. :	None.
;
; Written     :	Dominic Zarro (ARC)
;
; Version     :	Version 1.0, 18 September 1993
;-

function exist,var

return,n_elements(var) ne 0

end

