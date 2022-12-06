function str_lastpos, source, substring
;+
; Project     :	SOHO - CDS
;
; Name        :	STR_LASTPOS()
;
; Purpose     :	Find last occurence of a substring in the source string
;
; Category    :	String
;
; Explanation :	
;
; Syntax      :	Result = STR_LASTPOS( SOURCE, SUBSTRING )
;
; Examples    :	
;
; Inputs      :	SOURCE	  = String or string array to search
;		SUBSTRING = String to search for
;
; Opt. Inputs :	None.
;
; Outputs     :	Return value is position in string or -1 if not present.
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
; Prev. Hist. :	slf, 11/1/91
;		modified, 11/19/91 to allow string arrays for source
;
; History     :	Version 1, 15-Jan-1996, William Thompson, GSFC
;			Incorporated into SSW tree.
;
; Contact     :	WTHOMPSON
;-
;
;
rev_bytes=reverse(byte(source)) 
;
; since byte operation will generate nulls for array padding,
; terminators must be purged before converting back to string
;
terminators=where(rev_bytes eq 0)
if terminators(0) ge 0 then rev_bytes(terminators) = 32
;
rev_string=strtrim(string(rev_bytes),2)
rev_substring=string(reverse(byte(substring)))
;
; use standard strpos on now reversed operands
backpos=strpos(rev_string, rev_substring)
found=where(backpos ge 0)
if found(0) ge 0 then backpos(found) = $
   strlen(source(found)) - backpos(found) - strlen(substring(found))
;
return,backpos
end
