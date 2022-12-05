;+
; Project     : VSO
;
; Name        : FILE_TEST2
;
; Purpose     : Wrapper around FILE_TEST that supports /EXPAND_EVIRONMENT
;               which matches FILE_SEARCH. Note this keyword is the
;               default on Mac/Unix-OS but not Windows. 
;
; Category    : utility strings
;
; Syntax      : IDL> output=file_test2(input)
;
; Inputs      : INPUT = string array to test
;
; Outputs     : OUTPUT = boolean results
;
; Keywords    : EXPAND_ENVIRONMENT = expand environment variable
;
; History     : 21-Feb-2022, Zarro (ADNET)
;-

function file_test2,file,expand_environment=expand_environment,_extra=extra

;-- input sanity checks

if n_elements(file) eq 0 then return,0b
nf=n_elements(file)
if is_blank(file) then if nf eq 1 then return,0b else return,bytarr(nf)

expand=keyword_set(expand_environment)

if expand && is_windows() then return,file_test(str_env(file),_extra=extra)

return,file_test(file,_extra=extra)

end
