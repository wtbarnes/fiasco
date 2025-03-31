;+
; Project     :	HESSI
;
; Name        :	is_struct
;
; Purpose     :	returns, 1/0 if valid/invalid input structure
;
; Category    :	Structure handling
;
; Syntax      : IDL> output=is_struct(input)
;
; Inputs      :	INPUT = input structure array
;
; Outputs     :	OUTPUT = 1/0
;
; Written     : Zarro (EITI/GSFC), 17 Sept 2001
;
; Contact     : dzarro@solar.stanford.edu
;-


function is_struct,input

sz=size(input)
return,sz[n_elements(sz)-2] eq 8

end
