;+
; Project     : SOHO - CDS     
;                   
; Name        : DEFAULT
;               
; Purpose     : Supply default values for variables
;               
; Explanation : If the first parameter is not defined, it is
;		set to the value of the second parameter.
;               
; Use         : DEFAULT,VARIABLE,DEFAULT_VALUE
;    
; Inputs      : VARIABLE : The variable that could take on the default value
;
;		DEFAULT_VALUE : The default value.
;               
; Opt. Inputs : None.
;               
; Outputs     : None.
;               
; Opt. Outputs: None.
;               
; Keywords    : None.
;
; Calls       : None.
;
; Common      : None.
;               
; Restrictions: None.
;               
; Side effects: None.
;               
; Category    : Utility, Misc.
;               
; Prev. Hist. : Taken from my private library.
;
; Written     : Stein Vidar Hagfors Haugan
;               
; Modified    : Never
;
; Version     : 1, 4-Sept-1995
;-            

PRO DEFAULT,VAR,VAL

If N_params() lt 2 then message,"Use: DEFAULT,VARIABLE,DEFAULT_VALUE"

IF N_ELEMENTS(VAR) EQ 0 THEN VAR=VAL

END
