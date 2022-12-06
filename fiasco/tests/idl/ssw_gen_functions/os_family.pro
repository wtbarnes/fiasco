function OS_FAMILY,LOWER=LOWER,NAME=NAME
;+
; Project     :	SOHO - CDS
;
; Name        :	OS_FAMILY()
;
; Purpose     :	Return current operating system as in !VERSION.OS_FAMILY
;
; Category    :	Utilities, Operating_system
;
; Explanation :	Return the current operating system as in !VERSION.OS_FAMILY
;
;	OS_FAMILY is assumed to be 'unix' if !VERSION.OS is not 'windows',
;		'MacOS' or 'vms'
;
;	To make procedures from IDL V4.0 and later compatibile with earlier
;	versions of IDL, replace calls to !VERSION.OS_FAMILY with OS_FAMILY().
;
;	Can also be used to replace calls to !VERSION.OS if care is taken with
;	the change of case between 'Windows', which is what this routine
;	returns, and 'windows' which is what !VERSION.OS returned in versions
;	of IDL prior to 4.0.
;
; Syntax      :	Result = OS_FAMILY()
;
; Examples    :	IF OS_FAMILY() EQ 'vms' THEN ...
;
; Inputs      :	None.
;
; Outputs     :	The result of the function is a scalar string containing one of
;		the four values 'Windows','MacOS','vms' or 'unix'
;
; Opt. Outputs:	None.
;
; Keywords    :	LOWER - set to return lowercase strings
;               NAME - return OS_NAME (if available)
;
; Prev. Hist. :	Written,  W. Landsman
;
; History     :	Version 1, 29-Aug-1995, William Thompson, GSFC
;			Incorporated into CDS library
;               Version 2, 15 May 2000, Zarro (SM&A/GSFC) - added /LOWER
;         
;               Version 3, 22-Dec-2002, Zarro (EER/GSFC) - saved result in
;                        common block for speed
;               Modifed, 2-June-2006, Zarro (L-3Com/GSFC) - removed
;                        COMMON and added NAME keyword
;
; Contact     :	WTHOMPSON
;-

 lower=keyword_set(lower)
 name=''
 if have_tag(!version,'os_name',index) then begin
  name=!version.(index)
  if lower then name=strlowcase(name)
 endif

 if have_tag(!version,'os_family',index) then begin
  os=!version.(index)
  if lower then return,strlowcase(os) else return,os
 endif

;-- pre version 5

 os='unix'
 if have_tag(!version,'os',index) then begin
  os=!version.(index)
  if os eq 'windows' then os='Windows'
 endif

 if lower then return,strlowcase(os) else return,os

 end
