;+
; PROJECT:
;       SOHO - CDS
;
; NAME:
;       GREP
;
; PURPOSE:
;       Search for string through a string array (cf grep in Perl)
;
; CATEGORY: String utility
;
; CALLING SEQUENCE:
;       Result = GREP(str, strarray)
;
; INPUTS:
;       STR      - Target string to be searched for
;       STRARRAY - String array over which the string is seached
;
; OPTIONAL INPUTS:
;       None.
;
; OUTPUTS:
;       Result   - String array of elements from STRARRAY which contains
;                  the searched string STR. A null string will be returned if
;                  no elements in STRARRAY match STR.
;
; OPTIONAL OUTPUTS:
;       INDEX    - Index of elements in STRARRAY which match the target string
;
; KEYWORD PARAMETERS:
;       SENSITIVE - Make the search pattern case sensitive, if set
;       EXACT     - Requires exact match (instead of substring match)
;       EXCLUDE   - Set this keyword to return string array that does 
;                   not contain string STR
;       START     - match at start of string
;
; PREVIOUS HISTORY:
;       Written December 15, 1994, Liyun Wang, GSFC/ARC
;
; MODIFICATION HISTORY:
;       Version 2, Liyun Wang, GSFC/ARC, January 20, 1995
;          Added the INDEX keyword
;       Version 3, July 26, 1995, Liyun Wang, GSFC/ARC
;          Made it return a null string if error occurs (instead of
;             stopping the program.
;       Version 4, April 7, 1996, Liyun Wang, GSFC/ARC
;          Used intrinsic function WHERE; thanks to Richard Schwartz
;       Version 5, March 12, 1997, Liyun Wang, NASA/GSFC
;          Added the EXCLUDE keyword
;       Version 6, August 22, 2000, Zarro (EIT/GSFC)
;          Removed DATATYPE calls
;       Version 7, June 19, 2002, Zarro (LAC/GSFC)
;          Added /START
;-
;
FUNCTION GREP, str, strarray, sensitive=sensitive, exact=exact, index=index,$
               exclude=exclude, start=start

   sz=size(str)
   dtype=sz(n_elements(sz)-2)
   sz=size(strarray)
   stype=sz(n_elements(sz)-2)

   IF (dtype ne 7) OR (stype ne 7) THEN BEGIN
      MESSAGE, 'Input parameter must be of string type.', /cont
      RETURN, ''
   ENDIF

   aa = strtrim(strarray,2)
   a = strtrim(str,2)
   start=KEYWORD_SET(start)
   exact=KEYWORD_SET(exact)

   IF NOT KEYWORD_SET(sensitive) THEN BEGIN
    aa = STRUPCASE(aa)
    a = STRUPCASE(a)
   ENDIF

   IF KEYWORD_SET(exclude) THEN BEGIN 
    IF exact THEN index = WHERE(aa NE a) ELSE $
     index = WHERE(STRPOS(aa, a) EQ -1)
   ENDIF ELSE BEGIN 
    IF exact THEN index = WHERE(aa EQ a) ELSE BEGIN
     IF start then index=WHERE(STRPOS(aa,a) eq 0) ELSE $
      index = WHERE(STRPOS(aa, a) NE -1)
    ENDELSE
   ENDELSE

;---------------------------------------------------------------------------
;        Make sure the returned values are of scalar type
;---------------------------------------------------------------------------

   IF index(0) GE 0 THEN BEGIN
    IF N_ELEMENTS(index) EQ 1 THEN BEGIN
     index = index(0)
     RETURN, (strarray(index))(0)
    ENDIF ELSE RETURN, strarray(index)
   ENDIF ELSE BEGIN
    index = -1
    RETURN, ''
   ENDELSE

   END